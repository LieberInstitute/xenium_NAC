# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#Goal: Explore resolutions of clustering. 
#All clusters are at lambda of 0 (non-spatial/cell typing mode)

#load libraries
library(SingleCellExperiment)
library(escheR)
library(dplyr)
library(here)

#How many files? 
files <- list.files(path = here("processed-data","05_Clustering","Banksy_Results","NonSpatial"),
           pattern = "Lambda0_*",full.names = TRUE)

#Read in the files
files_list <- lapply(X = files,FUN = read.csv)
names(files_list) <- basename(files)
names_cols <-  c("row","cluster","id_cell")
files_list <- lapply(files_list, setNames, names_cols)

#How many clusters in each set
p <- lapply(files_list,FUN = function(x){
  length(unique(x$cluster)) 
}) %>% as.data.frame() %>% 
  t() %>% 
  as.data.frame() %>% 
  select("Number_of_clusters" = 1) %>% 
  tibble::rownames_to_column("Settings") %>% 
  mutate(Settings = as.character(lapply(strsplit(Settings,split = ".csv"),"[",1))) %>%
  ggplot(aes(x = Settings, y = Number_of_clusters)) +
  geom_bar(stat = "identity") +
  ggtitle("Number of Clusters Per Resolution Value\nNon-spatial clustering") +
  labs(y = "Number of Clusters") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5))

ggsave(plot = p,
       file = here("plots","05_clustering","Banksy","Non_spatial_Number_clusters.pdf"))


#make a dataframe of clusters
#First make the names of the lists the second column
for(i in names(files_list)){
  print(i)
  colnames(files_list[[i]])[2] <- as.character(lapply(strsplit(i,split = ".csv"),"[",1))
}

#Are they all in the same order
all(sapply(list(files_list[[2]]$id_cell, files_list[[3]]$id_cell,
                files_list[[4]]$id_cell, files_list[[5]]$id_cell,
                files_list[[6]]$id_cell), FUN = identical, files_list[[1]]$id_cell))
#[1] TRUE

#Now create a dataframe that is the cell_id and all of the clustering options? 
clustering_options <- data.frame(id_cell = files_list[[1]]$id_cell,
                                 Lambda0_res0.25_louvain  = files_list[[1]]$Lambda0_res0.25_louvain,
                                 Lambda0_res0.5_louvain   = files_list[[2]]$Lambda0_res0.5_louvain,
                                 Lambda0_res0.75_louvain  = files_list[[3]]$Lambda0_res0.75_louvain,
                                 Lambda0_res1_louvain     = files_list[[4]]$Lambda0_res1_louvain,
                                 Lambda0_res1.25_louvain  = files_list[[5]]$Lambda0_res1.25_louvain,
                                 Lambda0_res1.50_louvain  = files_list[[6]]$Lambda0_res1.5_louvain)

saveRDS(object = clustering_options,
        file = here("processed-data","05_Clustering","Non_spatial_results_combined.RDS"))


#Load spe object
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

#Subset for Br6660
spe <- spe[,spe$Donor == "Br6660"]

spe
# class: SpatialExperiment 
# dim: 366 2425841 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2425841): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

identical(rownames(colData(spe)),clustering_options$id_cell)
#[1] TRUE

identical(colnames(spe),clustering_options$id_cell)
#[1] TRUE

colData(spe) <- cbind(colData(spe),clustering_options)

identical(rownames(colData(spe)),colData(spe)$id_cell)
#[1] TRUE

identical(colnames(spe),colData(spe)$id_cell)
#[1] TRUE
