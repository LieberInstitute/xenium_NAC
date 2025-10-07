#Goal: Load Banksy results and explore
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe

#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe

#Add Banksy clusters to the object. 
Banksy_clusters <- read.csv(here("processed-data","05_Clustering","Banksy_clusters.csv"))
head(Banksy_clusters)

#X is just rownames, remove it. 
Banksy_clusters <- Banksy_clusters[,-1]

#Rename the columns
colnames(Banksy_clusters) <- c("Banksy_Cluster","key")

identical(rownames(colData(spe)),Banksy_clusters$key)
#[1] TRUE

#Just cbind the clusters into the colData
colData(spe)$Banksy_Cluster <- as.factor(Banksy_clusters$Banksy_Cluster)


table(spe$Sample,spe$Banksy_Cluster)

#Build a dataframe of annotations. 
ann_df <- data.frame(Banksy_Cluster = c(1:15),
                     Annotation = c("WM_Astro_A","WM_Astro_B","MSN_A",
                                    "WM_A","Endothelial","Excitatory",
                                    "WM_B","Microglia","D1_Island_A",
                                    "Inh_NPY_PNOC","Inh_CHAT","D1_Island_B",
                                    "Ependymal","WM_C","Immune_Other"))

#Add annotation
spe$Annotation <- ann_df$Annotation[match(spe$Banksy_Cluster,
                                          ann_df$Banksy_Cluster)]

table(spe$Annotation)

table(spe$Banksy_Cluster)


###Plot Annotations on the tissue sections
#Create a color palette that makes sense. 
#make some new brain colors
cluster_cols <- Polychrome::createPalette(length(unique(spe$Annotation)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(cluster_cols) <- unique(spe$Annotation)

#Plot with escheR 
for (i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[, spe$Sample == i]
  
  p <-make_escheR(sub_spe) |>
    add_fill("Annotation")+
    ggtitle(i) +
    scale_fill_manual(values = cluster_cols) +
    theme(element_text(hjust = 0.5))
  
  png(filename = here("plots","05_clustering","Banksy","Post_Annotation",paste0(i,".png")),
      units = "in",height = 20,width = 20,res = 300)
  print(p)
  dev.off()
  
}


#Save the colors
saveRDS(cluster_cols,file = here("processed-data","05_Clustering","Banksy_cluster_cols.Rds"))

#Save the object
saveRDS(spe,file = here("processed-data","02_build_spe","SPEs","spe_Annotated.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
