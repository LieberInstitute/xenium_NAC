#Goal: Run Banksy spatial clustering with parameters equating to spatial clustering for domain finding
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(Banksy)
library(scran)
library(here)


#Read in the filtered spe object
spe_filtered_path <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_path)

#Remove any cells with 0 counts for all genes 
#Remove cells with 0 counts
spe <- spe[, colSums(counts(spe)) > 0]

#Sanity check 
table(colSums(counts(spe)) > 0)

spe

####HVGs 
message(Sys.time(), " | Selecting top HVGs")

n_hvgs <- 2000
gene_var <- modelGeneVar(spe)
top_hvgs <- getTopHVGs(gene_var, n = n_hvgs)
message(sprintf("  Selected %d HVGs", length(top_hvgs)))


###Subset the spe
message(Sys.time(), " | Subsetting SPE for the top 2000 HVGs")
spe <- spe[top_hvgs,]

spe

#######Banksy parameters
lambda <- 0.8 #0.8 for spatial
compute_agf <- TRUE 
use_agf <- TRUE
k_geom <- c(25,50)  #Increase k value to find larger domains. 

#split spe by sample
samples <- unique(spe$sample_id)
spe_list <- list()
for(i in samples){
  spe_list[[i]] <- spe[, spe$sample_id == i]
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

saveHDF5SummarizedExperiment(spe_joint, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_banksy_0.8_cell_level"
))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
