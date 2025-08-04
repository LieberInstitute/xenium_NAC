#Goal: Run Banksy spatial clustering with default parameters
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


#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe

#######Banksy parameters
lambda <- 0.8 #Per Banksy reference manual, "0.8 incorporates more spatial neighborhood" and is good for spatial domains
res <- 0.5 
compute_agf <- FALSE #Run simpler version of Banksy 
use_agf <- FALSE 
k_geom <- 15 #Banksy reference manual suggest that values from 15-30 work well. 
attr <- sprintf("clust_lam%s_k10_res%s",lambda, res)

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

#Save the Banksy matrix, PCA, and UMAP
#Banksy mat
message(paste0("Saving matrix - ", Sys.time()))
saveRDS(object = assay(spe_joint,"H0"),file = here("processed-data","05_Clustering","Banksy_matrix.Rds"))
#Banksy PCA
message(paste0("Saving PCA embedding - ", Sys.time()))
saveRDS(object = reducedDim(spe_joint,"PCA_M0_lam0.8"),file = here("processed-data","05_Clustering","Banksy_PCA_Embedding.Rds"))
#Banksy UMAP
message(paste0("Saving UMAP embedding - ", Sys.time()))
saveRDS(object = reducedDim(spe_joint,"UMAP_M0_lam0.8"),file = here("processed-data","05_Clustering","Banksy_UMAP_Embedding.Rds"))

#Run Banksy leiden clustering
message(paste0("Running Banksy clustering - ", Sys.time()))
spe_joint <- clusterBanksy(spe_joint, use_agf = use_agf, lambda = lambda, resolution = res, seed = 1000)
message(paste0("Finished Banksy clustering - ", Sys.time()))


#Save CSV of cluster assigments
#Pull cluster name
cluster_name <- colnames(colData(spe_joint))[grep('^clust_', colnames(colData(spe_joint)))]
cluster_assign <- cbind(colData(spe_joint)[,cluster_name],rownames(colData(spe_joint)))
head(cluster_assign)
write.csv(cluster_assign,here("processed-data","05_Clustering","Banksy_clusters.csv"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
