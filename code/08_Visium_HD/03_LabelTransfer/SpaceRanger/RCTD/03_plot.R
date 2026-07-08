#Goal: Add RCTD weights and predictions to the SPE
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(HDF5Array)
library(escheR)
library(here)
library(data.table)

## Load the full SPE
message("Loading SPE object - ", Sys.time())
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)
spe

#Load colData from 02 script  
res <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD_doublet_colData.Rds"
))

res

#replace colData
stopifnot(identical(rownames(res), colnames(spe)))
colData(spe) <- res

## ---- 1) Normalize RCTD weights and call max type on normalized signal ----

weight_cols <- grep("^rctd_weight_", names(colData(spe)), value = TRUE)
ct_names <- gsub("^rctd_weight_", "", weight_cols)

# Extract raw weights
weight_mat <- as.matrix(colData(spe)[, weight_cols])

# Normalize: each row sums to 1 (proportion of signal)
row_sums <- rowSums(weight_mat)
weight_mat_norm <- sweep(weight_mat, 1, row_sums, "/")

# Call max type based on NORMALIZED weights
spe$rctd_max_type <- ct_names[max.col(weight_mat_norm)]
spe$rctd_max_weight <- apply(weight_mat_norm, 1, max)

table(spe$rctd_max_type)

# Add normalized weights to colData with _prop suffix
colnames(weight_mat_norm) <- paste0(ct_names, "_prop")
colData(spe) <- cbind(colData(spe), weight_mat_norm)

# Rename raw weight columns (strip rctd_weight_ prefix)
names(colData(spe))[match(weight_cols, names(colData(spe)))] <- ct_names

# Remove cells with NA (filtered during RCTD)
message("Removing ", sum(is.na(spe$rctd_max_type)), " cells with NA RCTD results")
spe <- spe[, !is.na(spe$rctd_max_type)]
spe

## ---- 3) Save object with raw scores, weighted scores, cell id, and max type ----

raw_cols <- ct_names
prop_cols <- paste0(ct_names, "_prop")

rctd_results <- data.table(
  cell_id = colnames(spe),
  rctd_max_type = spe$rctd_max_type,
  rctd_max_weight = spe$rctd_max_weight
)

# Add raw scores
rctd_results <- cbind(rctd_results, as.data.table(as.matrix(colData(spe)[, raw_cols])))

# Add normalized scores
rctd_results <- cbind(rctd_results, as.data.table(as.matrix(colData(spe)[, prop_cols])))

saveRDS(rctd_results, here(
  "processed-data", "HD_Full_Analysis", "SPEs", "doublet_rctd_raw_and_normalized_weights.Rds"
))
message("Saved RCTD results object with ", nrow(rctd_results), " cells and ", ncol(rctd_results), " columns")

## ---- 2) Create plot directories dynamically, deleting if they exist ----

# Helper function to create/recreate a directory
make_plot_dir <- function(dir_path) {
  if (dir.exists(dir_path)) {
    unlink(dir_path, recursive = TRUE)
  }
  dir.create(dir_path, recursive = TRUE)
}

# Create per-sample dirs
for (sample in unique(spe$sample_id)) {
  make_plot_dir(here("plots", "HD_Full_Analysis", "LabelTransfer", "doublet_RCTD", sample))
}

# Create per-celltype dirs
for (ct in unique(spe$rctd_max_type)) {
  make_plot_dir(here("plots", "HD_Full_Analysis", "LabelTransfer", "doublet_RCTD", ct))
}

## ---- Plot normalized weights per cell type on tissue ----

for (sample in unique(spe$sample_id)) {
  message("Plotting normalized weights for sample: ", sample)
  dir_path <- here("plots", "HD_Full_Analysis", "LabelTransfer", "doublet_RCTD", sample)
  spe_sub <- spe[, spe$sample_id == sample]
  for (ct in unique(spe$rctd_max_type)) {
    ct_prop <- paste0(ct, "_prop")
    p <- make_escheR(spe_sub) |>
      add_fill(ct_prop) +
      scale_fill_gradientn(colors = c("lightgrey", "red", "black"))
    ggsave(filename = here(dir_path, paste0(ct, "_norm.png")),
           plot = p,
           width = 20, height = 14, dpi = 200)
  }
}

## ---- Plot max type (all cell types together) ----

load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda", verbose = TRUE)

for (sample in unique(spe$sample_id)) {
  message("Plotting max cell type for sample: ", sample)
  dir_path <- here("plots", "HD_Full_Analysis", "LabelTransfer", "doublet_RCTD", sample)
  spe_sub <- spe[, spe$sample_id == sample]
  p <- make_escheR(spe_sub) |>
    add_fill("rctd_max_type") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here(dir_path, paste0(sample, "_all_celltypes.png")),
         plot = p,
         width = 20, height = 14, dpi = 200)
}

## ---- Plot one predicted cell type at a time ----

for (ct in unique(spe$rctd_max_type)) {
  message("Plotting individual cell type: ", ct)
  dir_path <- here("plots", "HD_Full_Analysis", "LabelTransfer", "doublet_RCTD", ct)
  for (sample in unique(spe$sample_id)) {
    spe_sub <- spe[, spe$sample_id == sample]
    spe_sub$CellType <- ifelse(spe_sub$rctd_max_type == ct, ct, "Other")
    ct_colors <- c("lightgrey", cluster_cols[ct])
    names(ct_colors) <- c("Other", ct)
    p <- make_escheR(spe_sub) |>
      add_fill("CellType") +
      scale_fill_manual(values = ct_colors)
    ggsave(filename = here(dir_path, paste0(sample, ".png")),
           plot = p,
           width = 20, height = 14, dpi = 200)
  }
}

sessionInfo()
