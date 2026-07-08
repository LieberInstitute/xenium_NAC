#Goal: Run BANKSY leiden clustering at a given resolution and plot per sample
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(optparse)
library(escheR)
library(Banksy)
library(here)

# Parse resolution from job script
option_list <- list(
  make_option("--res", type = "double", help = "Resolution for clusterBanksy")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$res)) {
  stop("Error: --res must be provided")
}

res <- opt$res
message(sprintf("Running clusterBanksy with resolution=%.2f", res))

# Read in sfe
sfe <- loadHDF5SummarizedExperiment( here(
  "processed-data", "HD_Full_Analysis", 
  "sfe_banksy_0_cell_level_nonspatial"
))

sfe


samples <- unique(sfe$sample_id)

# Cluster
message(paste0("Running Banksy clustering - ", Sys.time()))
sfe <- clusterBanksy(sfe, use_agf = FALSE, lambda = 0,
                     algo = "louvain", resolution = res, seed = 1313)
message(paste0("Finished Banksy clustering - ", Sys.time()))

#Save cluster assignments
saveRDS(colData(sfe),
        here("processed-data", "HD_Full_Analysis","sr_banksy","NonSpatial_colData",
             paste0("sr_nonspatial_banksy_res_", res, ".Rds")))

# Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
