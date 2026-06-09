library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(here)
library(sf)

sample_info <- read.csv(
  here("processed-data", "visiumHD_sample_info_NAc.csv"),
  stringsAsFactors = FALSE
)
sample_ids <- sample_info$sample_id
 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "spaceranger_sfe")
 
## ---- load each sample's SFE ------------------------------------------------
sfe_list <- list()
for (sid in sample_ids) {
  f <- file.path(sfe_dir, paste0(sid, "_sfe.rds"))
  if (!file.exists(f)) {
    message("  skipping (no rds): ", sid)
    next
  }
  message(format(Sys.time()), " | reading ", sid)
  sfe_list[[sid]] <- readRDS(f)
}
stopifnot(length(sfe_list) > 0)
 
## ---- merge -----------------------------------------------------------------
#   All samples were run through the same Space Ranger reference, so the genes
#   (Ensembl IDs set in script 1) are identical and in the same order -- no
#   harmonization needed. cbind still enforces identical rownames, so it errors
#   loudly if a sample was ever run against a different reference.
merged <- do.call(cbind, sfe_list)
message(format(Sys.time()), " | merged: ",
        nrow(merged), " genes x ", ncol(merged), " cells across ",
        length(unique(merged$sample_id)), " samples; colGeometries: ",
        paste(colGeometryNames(merged), collapse = ", "))
 
out_path <- file.path(sfe_dir, "merged_sfe.rds")
saveRDS(merged, out_path)
message(format(Sys.time()), " | wrote ", out_path)
 
sessionInfo()
