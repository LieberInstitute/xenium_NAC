# Use metaneighbor to compare human, NHP, rat MSNs
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(scuttle)
library(here)

#############################################
###########     Human   #####################
#############################################
###### Visium-HD
# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")

sfe <- loadHDF5SummarizedExperiment(sfe_dir)

#Add spatial clusters and non-spatial clusters
#spatial_clust <- read.csv(here("processed-data", "HD_Full_Analysis",
#                               "sr_spatial_banksy_clusters_res_0.4.csv"))
#rownames(spatial_clust) <- spatial_clust$V2
#spatial_clust <- spatial_clust[colnames(sfe),]

#Add spatial data to sfe 
#stopifnot(identical(spatial_clust$V2,colnames(sfe)))

#sfe$spatial_0.4 <- spatial_clust$V1

#Annotate
#anno_df <- read.csv(here("processed-data","spatial_0.4_Annotations.csv"))
#anno_df[11,"Annotation"] <- "Hypo"
#sfe$Spatial_Domain <- anno_df[match(sfe$spatial_0.4,anno_df$spatial_0.4),"Annotation"]

## Pseudobulk by spatial domain
sfe_pseudo <- aggregateAcrossCells(
  sfe,
  DataFrame(
    snRNA_label = colData(sfe)$snRNA_label,
    sample_id = colData(sfe)$sample_id
  ))
dim(sfe_pseudo)


#Now save. 
sfe_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_snRNA_label_pseudo"
)
saveHDF5SummarizedExperiment(sfe_pseudo, dir = sfe_dir,as.sparse = TRUE,replace = TRUE)

sessionInfo()
