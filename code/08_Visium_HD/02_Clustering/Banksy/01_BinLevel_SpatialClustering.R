#Goal: Run Banksy spatial clustering with parameters equating to spatial clustering for domain finding
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#code modified from: 
#https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/03_clustering/01_spatial_domains_banksy.R and
#https://www.bioconductor.org/packages/release/bioc/vignettes/Banksy/inst/doc/multi-sample.html
library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(Banksy)
library(here)
library(scran)


#Read in the filtered spe object
spe_filtered_path <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_norm_binQC_filtered.Rds"
)

spe <- readRDS(spe_filtered_path)

#Remove any cells with 0 counts for all genes 
#Remove cells with 0 counts
spe <- spe[, colSums(counts(spe)) > 0]

#Sanity check 
table(colSums(counts(spe)) > 0)

spe


####HVGs 
message(Sys.time(), " | Selecting top HVGs")

n_hvgs <- 4000
gene_var <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(gene_var, n = n_hvgs)
message(sprintf("  Selected %d HVGs", length(top_hvgs)))

################################################################################
#   Subsample bins (stratified by sample)
################################################################################
#Running into memory errors. Subsampling bins will help. Can then use label transfer to transfer to surrounding bins. 

set.seed(1607)
subsample_pct <- 0.33  # 33% of bins per sample

samples <- unique(spe$sample_id)

subsample_idx <- c()
for (sid in samples) {
  sample_idx <- which(spe$sample_id == sid)
  n_keep <- round(length(sample_idx) * subsample_pct)
  subsample_idx <- c(subsample_idx, sample(sample_idx, n_keep))
}
subsample_idx <- sort(subsample_idx)

spe_sub <- spe[top_hvgs, subsample_idx]

message(dim(spe_sub), " | Dimensions of subsetted object")


message(sprintf(
  "  Subsampled: %d -> %d bins (%.1f%%)",
  ncol(spe), ncol(spe_sub), ncol(spe_sub) / ncol(spe) * 100
))

for (sid in samples) {
  message(sprintf(
    "    %s: %d -> %d bins",
    sid,
    sum(spe$sample_id == sid),
    sum(spe_sub$sample_id == sid)
  ))
}

#######Banksy parameters
lambda <- 0.8 #0.8 for spatial
compute_agf <- TRUE 
use_agf <- TRUE
k_geom <- 50  #Increase k value to find larger domains. 

#split spe by sample
samples <- unique(spe$sample_id)
spe_list <- list()
for(i in samples){
  spe_list[[i]] <- spe_sub[, spe_sub$sample_id == i]
}

#Banksy Step 1 
message(paste0("Running computeBanksy - ", Sys.time()))
spe_list <- lapply(spe_list, computeBanksy,
                   assay_name = "logcounts", compute_agf = compute_agf, k_geom = k_geom)
message(paste0("Finished computeBanksy - ", Sys.time()))

#Merge SPEs
spe_joint <- do.call(cbind, spe_list)
rm(spe_list)
gc()

#Now run PCA and UMAP on the Banksy matrix
message(paste0("Running PCA - ", Sys.time()))
spe_joint <- runBanksyPCA(spe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
message(paste0("Finished PCA - ", Sys.time()))

message(paste0("Running UMAP - ", Sys.time()))
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
message(paste0("Finished UMAP - ", Sys.time()))

################################################################################
#   Save subsampled BANKSY object
################################################################################

message(Sys.time(), " | Saving subsampled BANKSY object")

saveRDS(spe_joint, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_bin_banksy_subsampled.RDS"
))

message(Sys.time(), " | Saving subsampled spe object pre-BANSKY")
saveRDS(object = spe_sub,file = here("processed-data","HD_Full_Analysis","SPEs","subsample_8umbins.Rds"))

message(Sys.time(), " | Saving subsampled cellular IDs")
saveRDS(object = subsample_idx,file = here("processed-data","HD_Full_Analysis","Banksy","subsample_idx_8umbins.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()



