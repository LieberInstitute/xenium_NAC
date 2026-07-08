#Goal: Choose the neighDist parameter for CRAWDAD by visualizing cell-type
#   proportion histograms at several candidate distances. Run this BEFORE
#   the full 04_CRAWDAD.R analysis to pick the right neighDist.
#   Ref: https://github.com/JEFworks-Lab/CRAWDAD/blob/main/docs/choosing_neighDist.md
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(patchwork)
library(tidyverse)
library(HDF5Array)
library(jsonlite)
library(crawdad)
library(here)

################################################################################
#   Setup — load one sample to explore
################################################################################

message(Sys.time(), " | Loading SPE")

spe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
))

res <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD_colData.Rds"
))
stopifnot(identical(rownames(res), colnames(spe)))
colData(spe) <- res

weight_cols <- grep("^rctd_weight_", names(colData(spe)), value = TRUE)
weight_mat  <- as.matrix(colData(spe)[, weight_cols])
spe$rctd_max_type <- gsub(
  "^rctd_weight_", "",
  weight_cols[max.col(weight_mat, ties.method = "first")]
)
spe <- spe[, !is.na(spe$rctd_max_type)]

################################################################################
#   Get current sample from SLURM array
################################################################################
message(Sys.time(), " | Setting up")
samples <- unique(spe$sample_id)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message(Sys.time(), " | Processing sample: ", current_sample,
        " (", task_id, " of ", length(samples), ")")

spe_sub <- spe[, spe$sample_id == current_sample]
message("  Number of cells: ", ncol(spe_sub))

plot_dir <- here("plots", "HD_Full_Analysis", "CRAWDAD", "neighDist_exploration")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#   Candidate neighDist values in microns. We want to find the value where the
#   proportion histogram is NOT piled up at 0 (too small) or 100 (too large),
#   but has a spread distribution.
candidate_um <- c(10, 25, 50, 75, 100, 150, 200)

#   Cell types to visualize neighborhoods spatially with vizClusters().
#   Adjust these to match your exact RCTD label names.
ref_types <- c("DRD1_MSN_A",
               "DRD2_MSN_A")

#   Subsample size — 50k cells is plenty for proportion histograms
max_cells <- 50000

##############################################################################
#   Get mpp from scalefactors
##############################################################################

scalefactors_path <- here(
  "processed-data", "01_spaceranger", current_sample, "outs",
  "binned_outputs", "square_002um", "spatial", "scalefactors_json.json"
)

if (!file.exists(scalefactors_path)) {
  scalefactors_path <- here(
    "processed-data", "01_spaceranger", current_sample, "outs",
    "binned_outputs", "square_008um", "spatial", "scalefactors_json.json"
  )
}

scalefactors <- fromJSON(scalefactors_path)
mpp <- scalefactors$microns_per_pixel
message("  microns_per_pixel: ", mpp)

##############################################################################
#   Subsample for speed
##############################################################################
message(Sys.time(), " | Subsampling")
coords <- spatialCoords(spe_sub)
pos <- data.frame(
  x        = coords[, "pxl_col_in_fullres"],
  y        = coords[, "pxl_row_in_fullres"],
  celltype = spe_sub$rctd_max_type
)

if (nrow(pos) > max_cells) {
  message("  Subsampling to ", max_cells, " cells for exploration")
  set.seed(1554)
  pos <- pos[sample(nrow(pos), max_cells), ]
}

message("  Cell type distribution:")
print(table(pos$celltype))

##############################################################################
#   Convert to sf
##############################################################################
message(Sys.time(), " | Converting to SF")
cells <- crawdad::toSF(
  pos       = pos[, c("x", "y")],
  cellTypes = pos$celltype
)

##############################################################################
#   Test neighDist candidates
##############################################################################

message("  Testing neighDist candidates: ",
        paste(candidate_um, "um", collapse = ", "))

plots <- list()

for (nd_um in candidate_um) {
  nd_px <- round(nd_um / mpp)
  message("    neighDist = ", nd_um, " µm (", nd_px, " px)")
  
  p <- vizCelltypeProportions(cells, neighDist = nd_px) +
    labs(title = paste0("neighDist = ", nd_um, " um (", nd_px, " px)")) +
    theme(plot.title = element_text(size = 10, hjust = 0.5))
  
  plots[[as.character(nd_um)]] <- p
  
  ggsave(
    file.path(plot_dir, paste0(
      current_sample, "_neighDist_", nd_um, "um.png"
    )),
    p, width = 10, height = 6, dpi = 150
  )
}

#   Combined comparison panel
p_combined <- wrap_plots(plots, ncol = 2) +
  plot_annotation(
    title    = paste0(current_sample, " — neighDist exploration"),
    subtitle = "Choose the value with a well-spread histogram (not piled at 0 or 100)"
  )

ggsave(
  file.path(plot_dir, paste0(current_sample, "_neighDist_comparison.png")),
  p_combined, width = 16, height = 4 * ceiling(length(candidate_um) / 2),
  dpi = 200
)

##############################################################################
#   Visualize neighborhoods in tissue space
##############################################################################

for (ref in ref_types) {
  if (!(ref %in% unique(pos$celltype))) {
    message("  Skipping vizClusters for '", ref, "' — not found in cell types")
    next
  }
  
  for (nd_um in c(25, 50, 100)) {
    nd_px <- round(nd_um / mpp)
    p <- vizClusters(cells, ref = ref, neighDist = nd_px, lineWidth = 0.5) +
      labs(title = paste0(
        current_sample, " — ", ref, " neighborhood — neighDist = ", nd_um, " µm"
      ))
    
    ggsave(
      file.path(plot_dir, paste0(
        current_sample, "_vizClusters_", ref, "_", nd_um, "um.png"
      )),
      p, width = 10, height = 8, dpi = 200
    )
  }
}

message(Sys.time(), " | Done with sample: ", current_sample)

#### Reproducibility
sessioninfo::session_info()
