#Goal: Run clusterBanksy on the subsampled BANKSY object that already has
#   PCA and UMAP computed from script 01.
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(Banksy)
library(here)

################################################################################
#   Setup
################################################################################

message(Sys.time(), " | Loading subsampled BANKSY object")

spe_joint <- readRDS(here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_bin_banksy_subsampled.RDS"
))

message(sprintf("Loaded: %d genes x %d bins", nrow(spe_joint), ncol(spe_joint)))

################################################################################
#   Run clustering
################################################################################

lambda <- 0.8
use_agf <- TRUE

message(paste0("Running clusterBanksy - ", Sys.time()))
spe_joint <- clusterBanksy(
  spe_joint,
  use_agf = use_agf,
  lambda = lambda,
  algo = "leiden",
  resolution = c(0.5, 0.8, 1.0),
  seed = 2040
)
message(paste0("Finished clusterBanksy - ", Sys.time()))

clust_cols <- grep("^clust_", colnames(colData(spe_joint)), value = TRUE)
message("Cluster columns created: ", paste(clust_cols, collapse = ", "))

for (clust_col in clust_cols) {
  message(sprintf("  %s: %d clusters", clust_col, nlevels(spe_joint[[clust_col]])))
}

################################################################################
#   Save
################################################################################

message(Sys.time(), " | Saving")

saveRDS(spe_joint, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_bin_banksy_subsampled_clustered.RDS"
))

message("Done.")

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
