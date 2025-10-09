#Run louvain clustering with a variety of resolution values. 
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
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


#Create a list that contains all of the clustering results from 02_nonspatial
files <- list.files(path = here("processed-data","05_Clustering",
                                "Banksy_Results","NonSpatial"),
	 	    pattern = "*.csv",
                    full.names = TRUE)

print(files)

file_list <- lapply(X = files,FUN = function(x) read.csv(x))

names(file_list) <- as.character(lapply(strsplit(x = basename(files),split = "_"),"[",2))

lapply(X = file_list,FUN = function(x) table(x$V1))



#Is everything in the same order? 
lapply(X = file_list,FUN = function(x){
  identical(colnames(spe),x$V2)
})



#Add all of the clusters to the object
spe$clust_res0.25 <- as.factor(file_list[["res0.25"]]$V1)
spe$clust_res0.50 <- as.factor(file_list[["res0.5"]]$V1)
spe$clust_res0.75 <- as.factor(file_list[["res0.75"]]$V1)
spe$clust_res1    <- as.factor(file_list[["res1"]]$V1)
spe$clust_res1.25 <- as.factor(file_list[["res1.25"]]$V1)
spe$clust_res1.50 <- as.factor(file_list[["res1.25"]]$V1)

#Now create escheR plots for each clustering
library(escheR)
for(i in c("clust_res0.25","clust_res0.50","clust_res0.75",
           "clust_res1","clust_res1.25","clust_res1.50")){
  new_dir_path <- here("plots","05_clustering","Banksy","Pre_Annotation","NonSpatial",i)
  dir.create(new_dir_path)
  
  #Generate cluster_cols
  cluster_cols <- Polychrome::createPalette(length(unique(colData(spe)[,i])),
                                            c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
  names(cluster_cols) <- as.character(unique(colData(spe)[,i]))
  saveRDS(cluster_cols,file = here("processed-data","05_Clustering","Banksy_Results","NonSpatial",paste0(i,"_nonspatial_colors.RDS")))
  
  for (l in unique(spe$Sample)){
    print(l)
    sub_spe <- spe[, spe$Sample == l]
    
    p <-make_escheR(sub_spe) |>
      add_fill(i)+
      scale_fill_manual(values = cluster_cols) +
      ggtitle(i) +
      theme(element_text(hjust = 0.5))
    
    png(filename = here(new_dir_path,paste0(l,"_",i,".png")),
        units = "in",height = 20,width = 20,res = 300)
    print(p)
    dev.off()}
  
}


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
