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

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe


#Keep only brain Br6660
spe <- spe[,spe$Donor == "Br6660"]

spe

unique(spe$Donor)

#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe

#Remove any cells with 0 counts for all genes in the panel
#Remove cells with 0 counts
spe <- spe[, colSums(counts(spe)) > 0]

#Sanity check 
table(colSums(counts(spe)) > 0)

spe


#######Banksy parameters
lambda <- 0.8 #0 for non-spatial
compute_agf <- TRUE #Run simpler version of Banksy 
use_agf <- TRUE
k_geom <- 30  #Banksy reference manual suggest that values from 15-30 work well. 

#split spe by sample
samples <- unique(spe$Sample)
spe_list <- list()
for(i in samples){
  spe_list[[i]] <- spe[,spe$Sample == i]
}


#Banksy Step 1 
message(paste0("Running computeBanksy - ", Sys.time()))
spe_list <- lapply(spe_list,computeBanksy,
                   assay_name = "nucleus_normcounts",compute_agf = compute_agf,k_geom = k_geom)
message(paste0("Finished computeBanksy - ", Sys.time()))


#Merge SPEs
spe_joint <- do.call(cbind, spe_list)
rm(spe_list)
gc()

#Now run PCA and UMAP on the Banksy matrix
message(paste0("Running PCA - ", Sys.time()))
spe_joint <- runBanksyPCA(spe_joint,use_agf = use_agf, lambda = lambda, group = "Sample",seed = 1000) #Vignette uses seed of 1000
message(paste0("Finished PCA - ", Sys.time()))

message(paste0("Running UMAP - ", Sys.time()))
spe_joint <- runBanksyUMAP(spe_joint, use_agf = use_agf, lambda = lambda, seed = 1000)
message(paste0("Finished UMAP - ", Sys.time()))

#Save the object
saveRDS(spe_joint, here("processed-data", "05_Clustering", "SPEs", "spe_joint_banksy_lambda_pt8.RDS"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
