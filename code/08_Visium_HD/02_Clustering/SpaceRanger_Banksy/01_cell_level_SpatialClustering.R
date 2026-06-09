#Goal: Run Banksy spatial clustering with parameters equating to spatial clustering for domain finding
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(Banksy)
library(scran)
library(here)


#Read in the filtered sfe object
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

#Remove any cells with 0 counts for all genes 
#Remove cells with 0 counts
sfe <- sfe[, colSums(counts(sfe)) > 0]

#Sanity check 
table(colSums(counts(sfe)) > 0)

sfe

####
#Calcualte log-normalized counts
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

####HVGs 
message(Sys.time(), " | Selecting top HVGs")

n_hvgs <- 2000
gene_var <- modelGeneVar(sfe)
top_hvgs <- getTopHVGs(gene_var, n = n_hvgs)
message(sprintf("  Selected %d HVGs", length(top_hvgs)))


###Subset the sfe
message(Sys.time(), " | Subsetting sfe for the top 2000 HVGs")
sfe <- sfe[top_hvgs,]

sfe

#######Banksy parameters
lambda <- 0.8 #0.8 for spatial
compute_agf <- TRUE 
use_agf <- TRUE
k_geom <- c(25,50)  #Increase k value to find larger domains. 

#split sfe by sample
samples <- unique(sfe$sample_id)
sfe_list <- list()
for(i in samples){
  sfe_list[[i]] <- sfe[, sfe$sample_id == i]
}

#Banksy Step 1 
message(paste0("Running computeBanksy - ", Sys.time()))
sfe_list <- lapply(sfe_list, computeBanksy,
                   assay_name = "logcounts", compute_agf = compute_agf, k_geom = k_geom)
message(paste0("Finished computeBanksy - ", Sys.time()))

#Merge sfes
sfe_joint <- do.call(cbind, sfe_list)
rm(sfe_list)
gc()

#Now run PCA and UMAP on the Banksy matrix
message(paste0("Running PCA - ", Sys.time()))
sfe_joint <- runBanksyPCA(sfe_joint, use_agf = use_agf, lambda = lambda, group = "sample_id", seed = 1000)
message(paste0("Finished PCA - ", Sys.time()))

message(paste0("Running UMAP - ", Sys.time()))
sfe_joint <- runBanksyUMAP(sfe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
message(paste0("Finished UMAP - ", Sys.time()))

################################################################################
#   Save subsampled BANKSY object
################################################################################

message(Sys.time(), " | Saving subsampled BANKSY object")

saveHDF5SummarizedExperiment(sfe_joint, 
                             here("processed-data", "HD_Full_Analysis","sfe_banksy_0.8_cell_level"),
                             replace = TRUE)

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
