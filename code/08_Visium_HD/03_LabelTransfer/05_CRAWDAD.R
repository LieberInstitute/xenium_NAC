#Goal: Run CRAWDAD multi-scale spatial colocalization on cell-level Visium-HD data
#   Uses RCTD cell type labels from 02_addRCTD / 03_plot pipeline
#   Runs per-sample as a SLURM array job
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(HDF5Array)
library(tidyverse)
library(jsonlite)
library(crawdad)
library(here)

################################################################################
#   Setup
################################################################################

message(Sys.time(), " | Loading SPE")

spe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)
)

#   Load RCTD colData and attach (same approach as 03_plot.R)
res <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD_colData.Rds"
))
stopifnot(identical(rownames(res), colnames(spe)))
colData(spe) <- res

#   Assign max RCTD cell type
weight_cols <- grep("^rctd_weight_", names(colData(spe)), value = TRUE)
weight_mat  <- as.matrix(colData(spe)[, weight_cols])
spe$rctd_max_type <- gsub(
  "^rctd_weight_", "",
  weight_cols[max.col(weight_mat, ties.method = "first")]
)

#   Remove the 51 NA cells (as in 03_plot.R)
spe <- spe[, !is.na(spe$rctd_max_type)]

################################################################################
#   Get current sample from SLURM array
################################################################################

samples   <- unique(spe$sample_id)
task_id   <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message(Sys.time(), " | Processing sample: ", current_sample,
        " (", task_id, " of ", length(samples), ")")

spe_sub <- spe[, spe$sample_id == current_sample]
message("  Number of cells: ", ncol(spe_sub))

################################################################################
#   Determine pixel-to-micron conversion from scalefactors
################################################################################

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

################################################################################
#   Define scales and neighDist (in pixel units to match spatial coords)
################################################################################

#   Scales for shuffling — biologically meaningful distances in the NAc:
#     ~50 µm  : local microenvironment / dendritic arbor scale
#    ~100 µm  : small cellular neighborhood
#    ~200 µm  : medium-range interactions
#    ~300 µm  : subregional transitions
#    ~400 µm  : patch/matrix compartment scale
#    ~500 µm  : NAc core-shell boundary zone
#    ~750 µm  : large subregion
#   ~1000 µm  : macro-scale tissue organization

distances_um <- c(50, 100, 200, 300, 400, 500, 750, 1000)
scales_px    <- round(distances_um / mpp)

message("  Scales in µm: ", paste(distances_um, collapse = ", "))
message("  Scales in px: ", paste(scales_px, collapse = ", "))

#   neighDist: radius defining the local neighborhood around each reference
#   cell for proportion counting. Should be tuned with vizCelltypeProportions()
#   so that the proportion histogram is not piled up at 0 or 100.
#   Start with ~50 µm as a reasonable default for cell-level Visium-HD data.
neighDist_px <- round(100 / mpp) #03 script showed 100 is best. 
message("  neighDist in px: ", neighDist_px, " (~50 µm)")

################################################################################
#   (Optional) Subsample if needed for memory
################################################################################

max_cells <- 300000
coords <- spatialCoords(spe_sub)

pos <- data.frame(
  x        = coords[, "pxl_col_in_fullres"],
  y        = coords[, "pxl_row_in_fullres"],
  celltype = spe_sub$rctd_max_type
)

if (nrow(pos) > max_cells) {
  message(Sys.time(), " | Subsampling from ", nrow(pos), " to ", max_cells, " cells")
  set.seed(42)
  keep_idx <- sample(nrow(pos), max_cells)
  pos <- pos[keep_idx, ]
}

message("  Cell type distribution:")
print(table(pos$celltype))

################################################################################
#   Convert to sf object
################################################################################

message(Sys.time(), " | Converting to sf object with toSF()")

cells <- crawdad::toSF(
  pos       = pos[, c("x", "y")],
  cellTypes = pos$celltype
)

################################################################################
#   Create shuffled null background
################################################################################

message(Sys.time(), " | Making shuffled cells for null background")

shuffle_list <- crawdad:::makeShuffledCells(
  cells,
  scales  = scales_px,
  perms   = 3,        
  ncores  = 1,
  seed    = 1643,
  verbose = TRUE
)

################################################################################
#   Run CRAWDAD — compute z-scores across scales
################################################################################

message(Sys.time(), " | Running findTrends()")

results <- crawdad::findTrends(
  cells,
  neighDist   = neighDist_px,
  shuffleList = shuffle_list,
  ncores      = 1,
  verbose     = TRUE,
  returnMeans = FALSE
)

################################################################################
#   Process results
################################################################################

message(Sys.time(), " | Processing results")

dat <- crawdad::meltResultsList(results, withPerms = TRUE)
dat$sample_id <- current_sample

#   Bonferroni-corrected z-score threshold
zsig <- crawdad::correctZBonferroni(dat)
message("  Bonferroni z-score threshold: ", round(zsig, 3))

################################################################################
#   Save results
################################################################################

outdir <- here("processed-data", "HD_Full_Analysis", "CRAWDAD")

saveRDS(
  list(dat = dat, zsig = zsig, results_raw = results,
       neighDist_px = neighDist_px, scales_px = scales_px, mpp = mpp),
  file.path(outdir, paste0("crawdad_results_", current_sample, ".Rds"))
)

################################################################################
#   Generate summary plots
################################################################################

message(Sys.time(), " | Generating plots")


plot_dir <- here("plots","HD_Full_Analysis","CRAWDAD")

#   Colocalization dotplot (main CRAWDAD summary)
p_dot <- vizColocDotplot(
  dat,
  zSigThresh  = zsig,
  zScoreLimit = 2 * zsig,
  dotSizes    = c(3, 15)
) +
  theme(axis.text.x = element_text(angle = 35, hjust = 0)) +
  labs(title = current_sample)

ggsave(
  file.path(plot_dir, paste0(current_sample, "_colocalization_dotplot.png")),
  p_dot, width = 16, height = 12, dpi = 200
)

message(Sys.time(), " | Done with sample: ", current_sample)

################################################################################
#   Reproducibility
################################################################################

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
