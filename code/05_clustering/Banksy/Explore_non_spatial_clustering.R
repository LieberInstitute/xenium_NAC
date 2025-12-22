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


#Plot each cluster on the tissue section. 
#It's best to use Br6660_Nac10_4080 becaues the islands are very obvious where the islands are
for(i in c("Lambda0_res0.25_louvain","Lambda0_res0.5_louvain","Lambda0_res0.75_louvain",
           "Lambda0_res1_louvain","Lambda0_res1.25_louvain","Lambda0_res1.50_louvain")){
  print(i)
  
  #Make a new directory
  dir.create(path = here("plots","05_clustering","Banksy",i))

  #Susbet for the one depth
  sub_spe <- spe[,spe$Sample == "Br6660_Nac10_4080"]

  for(cluster in unique(colData(sub_spe)[,i])){
    print(cluster)
    sub_spe$cluster <- ifelse(colData(sub_spe)[,i] == cluster,
                              cluster,
                              "Other")

    cluster_cols <- c("red","grey")
    names(cluster_cols) <- c(cluster,"Other")

    p <- make_escheR(sub_spe) |>
      add_fill("cluster") +
      scale_fill_manual(values = cluster_cols)
    
    ggsave(plot = p, filename = here("plots","05_clustering","Banksy",i,paste0(i,"_",cluster,".png")),
           height = 24, width = 24)
  }
  
}

#Load in annotation_df
anno_df <- read.csv(here("processed-data","05_Clustering","nonspatial_lambda0_res0.5_annotation.csv"))

spe$CellType <- anno_df$Annotation[match(spe$Lambda0_res0.5_louvain,anno_df$Cluster)]

#We want to know what percentage of the annotated celltypes 
for(i in c("Lambda0_res0.25_louvain","Lambda0_res0.5_louvain","Lambda0_res0.75_louvain",
           "Lambda0_res1_louvain","Lambda0_res1.25_louvain","Lambda0_res1.50_louvain")){
  print(i)
  x <- table(colData(spe)[,i],spe$CellType) %>% 
    as.data.frame.matrix() %>%
    sweep(MARGIN = 2,STATS = colSums(.),FUN = "/") *100
  
  
  pdf(file = here("plots","05_clustering","Banksy",paste0(i,"_byCellType_Heatmap.pdf")),
      height = 12, width = 12)
  pheatmap::pheatmap(x)
  dev.off()
}
sessionInfo()
# R version 4.5.0 Patched (2025-05-21 r88220)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] here_1.0.1                  dplyr_1.1.4                
# [3] escheR_1.8.0                ggplot2_3.5.2              
# [5] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [7] Biobase_2.68.0              GenomicRanges_1.60.0       
# [9] GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [11] S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [13] generics_0.1.4              MatrixGenerics_1.20.0      
# [15] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] SparseArray_1.8.0        lattice_0.22-7           magrittr_2.0.3          
# [4] grid_4.5.0               RColorBrewer_1.1-3       rprojroot_2.0.4         
# [7] jsonlite_2.0.0           Matrix_1.7-3             httr_1.4.7              
# [10] viridisLite_0.4.2        UCSC.utils_1.4.0         scales_1.4.0            
# [13] textshaping_1.0.1        abind_1.4-8              cli_3.6.5               
# [16] rlang_1.1.6              crayon_1.5.3             XVector_0.48.0          
# [19] withr_3.0.2              DelayedArray_0.34.1      S4Arrays_1.8.0          
# [22] tools_4.5.0              SpatialExperiment_1.18.1 GenomeInfoDbData_1.2.14 
# [25] vctrs_0.6.5              R6_2.6.1                 lifecycle_1.0.4         
# [28] magick_2.8.6             ragg_1.4.0               pkgconfig_2.0.3         
# [31] pillar_1.10.2            gtable_0.3.6             glue_1.8.0              
# [34] Rcpp_1.0.14              systemfonts_1.2.3        tibble_3.2.1            
# [37] tidyselect_1.2.1         dichromat_2.0-0.1        farver_2.1.2            
# [40] rjson_0.2.23             labeling_0.4.3           pheatmap_1.0.12         
# [43] compiler_4.5.0          

