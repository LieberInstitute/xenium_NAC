# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(escheR)
library(Seurat)
library(here)

#Load in the count normalized spe object
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe
# class: SpatialExperiment 
# dim: 366 4849373 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4849373): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#load in the snRNA-seq predictions
Br6660_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_anno_predictions.Rds"))
Br6436_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_6436_snRNA_predictions.Rds"))
predictions <- rbind(Br6660_predictions,Br6436_predictions)


#Sanity checks before adding 
identical(nrow(predictions),ncol(spe))
#[1] TRUE

identical(rownames(predictions),rownames(colData(spe)))
#[1] TRUE

identical(rownames(predictions),colnames(spe))
#[1] TRUE

#Prediction dataframes have the predicted id and the weighted probability score for each designation. 
#Add that information to the colData of the spe object. prediction.score.max
colData(spe)$snRNA_predicted_CellType <- predictions[,"predicted.id"]
colData(spe)$predicted_max_score <- predictions[,"prediction.score.max"]

table(colData(spe)$snRNA_predicted_CellType)
# Astrocyte_A Astrocyte_B  DRD1_MSN_A  DRD1_MSN_B  DRD1_MSN_C  DRD1_MSN_D 
# 556353      150080     1109211       78429       55237       26038 
# DRD2_MSN_A  DRD2_MSN_B Endothelial   Ependymal  Excitatory       Inh_A 
# 306335       19920      426859       25396      144365       27353 
# Inh_B       Inh_C       Inh_D       Inh_E       Inh_F   Microglia 
# 13465       27044       21104       26348      146711      297154 
# Oligo         OPC 
# 1198472      193499 

#load cluster colors from snRNA-seq paper. 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE)
# Loading objects:
#   cluster_cols

#Remove neuron ambig
cluster_cols <- cluster_cols[-14]
cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C  DRD1_MSN_D       Inh_B       Inh_C       Inh_D       Inh_E 
# "#2A4BF9"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281"   "#0DFBFA" 
# Excitatory       Inh_F 
# "#5C6300"     "black" 


spe$snRNA_predicted_CellType <- factor(x = spe$snRNA_predicted_CellType,
                                       levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                  "DRD2_MSN_A","DRD2_MSN_B",
                                                  "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                                  "Excitatory",
                                                  "Oligo","OPC",
                                                  "Astrocyte_A","Astrocyte_B","Ependymal",
                                                  "Microglia",
                                                  "Endothelial"))


library(ggplot2)

max_score_box <- ggplot(colData(spe),aes(x = snRNA_predicted_CellType, 
                                         y = predicted_max_score,
                                         fill = snRNA_predicted_CellType)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.50) +
  theme_bw() +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Predicted Cell Type",
       y = "Max Score")
ggsave(max_score_box,filename = here("plots","06_label_transfer","snRNA_LabelTransfer",
                                     "Max_score_boxplot.pdf"))
#Plot the predicted IDs on the tissue sections .
for(i in unique(spe$Sample)){
  print(i)
  spe_sub <- spe[,spe$Sample == i]
  p <- make_escheR(spe_sub) |>
    add_fill("snRNA_predicted_CellType") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here("plots","06_label_transfer","snRNA_LabelTransfer",paste0(i,".png")),
         height = 22, width = 22)
}


###Reproduciblity
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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] here_1.0.1                  Seurat_5.3.0                SeuratObject_5.1.0          sp_2.2-0                   
# [5] escheR_1.8.0                ggplot2_3.5.2               SpatialExperiment_1.18.1    SingleCellExperiment_1.30.1
# [9] SummarizedExperiment_1.38.1 Biobase_2.68.0              GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [13] IRanges_2.42.0              S4Vectors_0.46.0            BiocGenerics_0.54.0         generics_0.1.4             
# [17] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3      jsonlite_2.0.0          magrittr_2.0.3          spatstat.utils_3.1-4   
# [5] magick_2.8.6            farver_2.1.2            ragg_1.4.0              vctrs_0.6.5            
# [9] ROCR_1.0-11             spatstat.explore_3.4-3  htmltools_0.5.8.1       S4Arrays_1.8.0         
# [13] SparseArray_1.8.0       sctransform_0.4.2       parallelly_1.44.0       KernSmooth_2.23-26     
# [17] htmlwidgets_1.6.4       ica_1.0-3               plyr_1.8.9              plotly_4.10.4          
# [21] zoo_1.8-14              igraph_2.1.4            mime_0.13               lifecycle_1.0.4        
# [25] pkgconfig_2.0.3         Matrix_1.7-3            R6_2.6.1                fastmap_1.2.0          
# [29] GenomeInfoDbData_1.2.14 fitdistrplus_1.2-2      future_1.49.0           shiny_1.10.0           
# [33] digest_0.6.37           patchwork_1.3.0         rprojroot_2.0.4         tensor_1.5             
# [37] RSpectra_0.16-2         irlba_2.3.5.1           textshaping_1.0.1       labeling_0.4.3         
# [41] progressr_0.15.1        spatstat.sparse_3.1-0   httr_1.4.7              polyclip_1.10-7        
# [45] abind_1.4-8             compiler_4.5.0          withr_3.0.2             fastDummies_1.7.5      
# [49] MASS_7.3-65             sessioninfo_1.2.3       DelayedArray_0.34.1     rjson_0.2.23           
# [53] tools_4.5.0             lmtest_0.9-40           httpuv_1.6.16           future.apply_1.11.3    
# [57] goftest_1.2-3           glue_1.8.0              nlme_3.1-168            promises_1.3.2         
# [61] grid_4.5.0              Rtsne_0.17              cluster_2.1.8.1         reshape2_1.4.4         
# [65] gtable_0.3.6            spatstat.data_3.1-6     tidyr_1.3.1             data.table_1.17.2      
# [69] XVector_0.48.0          spatstat.geom_3.4-1     RcppAnnoy_0.0.22        ggrepel_0.9.6          
# [73] RANN_2.6.2              pillar_1.10.2           stringr_1.5.1           spam_2.11-1            
# [77] RcppHNSW_0.6.0          later_1.4.2             splines_4.5.0           dplyr_1.1.4            
# [81] lattice_0.22-7          survival_3.8-3          deldir_2.0-4            tidyselect_1.2.1       
# [85] miniUI_0.1.2            pbapply_1.7-2           gridExtra_2.3           scattermore_1.2        
# [89] stringi_1.8.7           UCSC.utils_1.4.0        lazyeval_0.2.2          codetools_0.2-20       
# [93] tibble_3.2.1            cli_3.6.5               uwot_0.2.3              systemfonts_1.2.3      
# [97] xtable_1.8-4            reticulate_1.42.0       dichromat_2.0-0.1       Rcpp_1.0.14            
# [101] globals_0.18.0          spatstat.random_3.4-1   png_0.1-8               spatstat.univar_3.1-3  
# [105] parallel_4.5.0          dotCall64_1.2           listenv_0.9.1           viridisLite_0.4.2      
# [109] scales_1.4.0            ggridges_0.5.6          purrr_1.0.4             crayon_1.5.3           
# [113] rlang_1.1.6             cowplot_1.1.3          
