# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(DeconvoBuddies)
library(sessioninfo)
library(HDF5Array)
library(scran)
library(here)

# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")

sfe <- loadHDF5SummarizedExperiment(sfe_dir)
#Add spatial clusters and non-spatial clusters
spatial_clust <- read.csv(here("processed-data", "HD_Full_Analysis",
                               "sr_spatial_banksy_clusters_res_0.4.csv"))
rownames(spatial_clust) <- spatial_clust$V2
spatial_clust <- spatial_clust[colnames(sfe),]

#Add spatial data to sfe 
stopifnot(identical(spatial_clust$V2,colnames(sfe)))

sfe$spatial_0.4 <- spatial_clust$V1

set.seed(1008)

# Group cells by cluster within sample
group <- interaction(sfe$sample_id, sfe$spatial_0.4, drop = TRUE)

# From each sample × cluster group, keep 33% of cells
keep_idx <- unlist(lapply(split(seq_len(ncol(sfe)), group), function(idx) {
  n_keep <- round(length(idx) * 0.33)
  if (n_keep < 1) n_keep <- 1   # guarantee at least 1 cell per non-empty group
  sample(idx, n_keep)
}), use.names = FALSE)

keep_idx <- sort(keep_idx)

# Subset the SFE
sfe <- sfe[, keep_idx]

###############
#Run get_mean_ratio to go ahead and 
message("Starting get_mean_ratio() - ",Sys.time())
gmr_spatial <- get_mean_ratio(
  sfe,
  "spatial_0.4",
  assay_name = "logcounts"
)

saveRDS(gmr_spatial,
        file = here("processed-data","HD_Full_Analysis","spatial_0.4_gmr_results.Rds"))

###############
message("Starting pairwise DEG testing - ",Sys.time())
mod <- with(colData(sfe), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(sfe, 
                                groups=sfe$spatial_0.4,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

saveRDS(markers_pairwise,file = here("processed-data","HD_Full_Analysis","spatial_0.4_pairwse.Rds"))

###############
#findMarkers_1vALL 
message("Starting 1vALL DEG testing - ",Sys.time())
markers_1vALL_enrich <- findMarkers_1vAll(sfe, 
                                          assay_name = "logcounts", 
                                          cellType_col = "spatial_0.4", 
                                          mod = "~sample_id")

saveRDS(markers_1vALL_enrich,file = here("processed-data","HD_Full_Analysis","spatial_0.4_1vALL.Rds"))



sessionInfo()
