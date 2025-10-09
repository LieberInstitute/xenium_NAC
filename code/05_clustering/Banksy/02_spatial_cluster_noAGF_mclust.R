#Run louvain clustering with a variety of resolution values. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
library(Banksy)
library(here)

#Read in the object
spe_joint <- readRDS(here("processed-data", "05_Clustering", "SPEs", "spe_joint_banksy_lambda_pt8_noAGF.RDS"))

spe_joint

#Perform banksy clustering. 
#Run Banksy louvain clustering
message(paste0("Running Banksy clustering - ", Sys.time()))
spe_joint <- clusterBanksy(spe_joint, use_agf = FALSE, lambda = 0.8,mclust.G = 8, algo =  "mclust",seed = 1000)
message(paste0("Finished Banksy clustering - ", Sys.time()))

#Save CSV of cluster assigments
#Pull cluster name
cluster_name <- colnames(colData(spe_joint))[grep('^clust_', colnames(colData(spe_joint)))]
cluster_assign <- cbind(colData(spe_joint)[,cluster_name],rownames(colData(spe_joint)))
head(cluster_assign)
write.csv(cluster_assign,here("processed-data","05_Clustering","Banksy_Results","NonSpatial","mclust_Lambdapt8.csv"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


