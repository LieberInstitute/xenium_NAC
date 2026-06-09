#Run louvain clustering with a variety of resolution values. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(Banksy)
library(here)
library(scran)

#############################################
############### SET UP OBJECT ###############
#############################################
##Read in the filtered spe object
#Read in the filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

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
lambda <- 0 #0 for non-spatial
compute_agf <- FALSE #Run simpler version of Banksy 
use_agf <- FALSE
k_geom <- 30  #Banksy reference manual suggest that values from 15-30 work well. 

#split spe by sample
samples <- unique(spe$sample_id)
spe_list <- list()
for(i in samples){
  spe_list[[i]] <- spe[,spe$sample_id == i]
}


#Banksy Step 1 
message(paste0("Running computeBanksy - ", Sys.time()))
spe_list <- lapply(spe_list,computeBanksy,
                   assay_name = "logcounts",compute_agf = compute_agf,k_geom = k_geom)
message(paste0("Finished computeBanksy - ", Sys.time()))

message("Assays after computeBanksy: ", paste(assayNames(spe_list[[1]]), collapse = ", "))
#Merge SPEs
spe_joint <- do.call(cbind, spe_list)
rm(spe_list)
gc()

#Now run PCA and UMAP on the Banksy matrix
message(paste0("Running PCA - ", Sys.time()))
spe_joint <- runBanksyPCA(spe_joint,use_agf = use_agf, lambda = lambda, group = "sample_id",seed = 1000) #Vignette uses seed of 1000
message(paste0("Finished PCA - ", Sys.time()))

message(paste0("Running UMAP - ", Sys.time()))
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
message(paste0("Finished UMAP - ", Sys.time()))

#Save the object
saveHDF5SummarizedExperiment(spe_joint, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_banksy_0_cell_level_nonspatial"
),replace = TRUE)


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
