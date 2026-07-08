#Goal: Add RCTD weights and predictions to the sfe
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(here)

## Load the full sfe
message("Loading SFE object - ", Sys.time())

#Read in the filtered sfe object
#Now save.
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

sfe

## Get cell type names from first available result
rctd_dir <- here("processed-data", "HD_Full_Analysis", "LabelTransfer","SpaceRanger","RCTD")

samples <- unique(sfe$sample_id)

## Get cell type names AND prediction column names from first available resultS
first_file <- file.path(rctd_dir, paste0("doublet_snRNA_VisiumHD_RCTD_", samples[1], ".Rds"))
first_results <- readRDS(first_file)
cell_types <- rownames(assay(first_results, "weights"))

## Pick the colData columns you want to carry over.
## Check colnames(colData(first_results)) to confirm what's there.
pred_cols <- c("spot_class", "first_type", "second_type")

## Initialize weight columns (numeric)
for (ct in cell_types) {
  sfe[[paste0("rctd_weight_", ct)]] <- NA_real_
}

## Initialize prediction columns (character)
for (pc in pred_cols) {
  sfe[[paste0("rctd_", pc)]] <- NA_character_
}

## Load and merge each sample's RCTD results
for (s in samples) {
  message(s, " - ", Sys.time())
  rctd_file <- file.path(rctd_dir, paste0("doublet_snRNA_VisiumHD_RCTD_", s, ".Rds"))
  
  if (!file.exists(rctd_file)) {
    message("Missing RCTD results for sample: ", s)
    next
  }
  
  message("Loading RCTD results for sample: ", s)
  res <- readRDS(rctd_file)
  wt  <- assay(res, "weights")
  cd  <- colData(res)
  
  ## Match barcodes
  shared <- intersect(colnames(sfe)[sfe$sample_id == s], colnames(res))
  idx <- match(shared, colnames(sfe))
  
  ## Transfer weights
  for (ct in cell_types) {
    sfe[[paste0("rctd_weight_", ct)]][idx] <- wt[ct, shared]
  }
  
  ## Transfer predictions
  for (pc in pred_cols) {
    sfe[[paste0("rctd_", pc)]][idx] <- as.character(cd[shared, pc])
  }
}

## Save updated sfe colData
message("Saving colData - ", Sys.time())
sfe <- as(sfe, "SpatialExperiment")
saveRDS(colData(sfe), here(
  "processed-data", "HD_Full_Analysis","sfe_cell_RCTD_doublet_colData.Rds"
))

## Save updated sfe
message("Saving sfe object - ", Sys.time())
saveHDF5SummarizedExperiment(sfe, here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_RCTD_doublet"
), replace = TRUE)

sessionInfo()
