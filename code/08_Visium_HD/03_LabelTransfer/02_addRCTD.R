#Goal: Add RCTD weights and predictions to the SPE
library(SpatialExperiment)
library(HDF5Array)
library(here)

## Load the full SPE
message("Loading SPE object - ", Sys.time())
spe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
))

## Get cell type names from first available result
rctd_dir <- here("processed-data", "HD_Full_Analysis", "LabelTransfer", "RCTD")
samples <- unique(spe$sample_id)

first_file <- file.path(rctd_dir, paste0("snRNA_VisiumHD_RCTD_", samples[1], ".Rds"))
first_results <- readRDS(first_file)
cell_types <- rownames(assay(first_results, "weights"))

## Initialize weight columns
for (ct in cell_types) {
  print(ct)
  spe[[paste0("rctd_weight_", ct)]] <- NA_real_
}

## Load and merge each sample's RCTD results
for (s in samples) {
  message(s," - ", Sys.time())
  rctd_file <- file.path(rctd_dir, paste0("snRNA_VisiumHD_RCTD_", s, ".Rds"))
  
  if (!file.exists(rctd_file)) {
    message("Missing RCTD results for sample: ", s)
    next
  }
  
  message("Loading RCTD results for sample: ", s)
  res <- readRDS(rctd_file)
  wt <- assay(res, "weights")
  
  ## Match barcodes
  shared <- intersect(colnames(spe)[spe$sample_id == s], colnames(res))
  idx <- match(shared, colnames(spe))
  
  ## Transfer weights
  for (ct in cell_types) {
    spe[[paste0("rctd_weight_", ct)]][idx] <- wt[ct, shared]
  }
}

## Save updated SPE colData
message("Saving colData - ", Sys.time())
spe <- as(spe, "SpatialExperiment")
saveRDS(colData(spe), here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD_colData.Rds"
))

## Save updated SPE
message("Saving SPE object - ", Sys.time())
saveHDF5SummarizedExperiment(spe, here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD"
), replace = TRUE)

sessionInfo()
