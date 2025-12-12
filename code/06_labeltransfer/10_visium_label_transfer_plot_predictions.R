#cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NACy
#module load conda_R/4.5
library(SpatialExperiment)
library(ComplexHeatmap)
library(escheR)
library(here)

#Load in the annotated object --> this contains annotations from Banksy clustering. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

dim(spe)
#[1]     366 4849373

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

#Read in the predicted spatial domains 
Br6660_domains <- readRDS(here("processed-data",
                               "06_label_transfer",
                               "seurat_anno_spatial_domains_predictions.Rds"))

Br6436_domains <- readRDS(here("processed-data",
                               "06_label_transfer",
                               "seurat_6436_spatial_domains_predictions.Rds"))

domain_df <- data.frame(cell_id = c(rownames(Br6660_domains),rownames(Br6436_domains)),
                        spatial_domain = c(Br6660_domains$predicted.id,Br6436_domains$predicted.id))

identical(domain_df$cell_id,colnames(spe))
#[1] TRUE

#Add domain info into the object
spe$Visium_domain <- domain_df$spatial_domain

#List of colors for the domains. 
domain_colors <- c("D1 islands"             = "#E83E8C",  
                   "Endothelial/Ependymal"  = "#B8862B",  
                   "Excitatory"             = "#E4572E",  
                   "Inhibitory"             = "#F2B134",  
                   "MSN 1"                  = "#6DBE45",  
                   "MSN 2"                  = "#1B9E77",  
                   "MSN 3"                  = "#7570B3",  
                   "WM"                     = "#6E6E6E")


for(i in unique(spe$Sample)){
  print(i)
  spe_sub <- spe[,spe$Sample == i]
  p <- make_escheR(spe_sub) |>
    add_fill("Visium_domain") +
    scale_fill_manual(values = domain_colors)
  ggsave(filename = here("plots","06_label_transfer",
                         "Visium_label_transfer",paste0(i,".png")),
         height = 22, width = 22)
}


#Add the Banksy cell types to the domain. 
#Load the annotated banksy cell types
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

identical(Banksy_celltypes$cell_id,colnames(spe))
#[1] TRUE

spe$Banksy_celltypes <- Banksy_celltypes$CellType


library(dplyr)

pct_mat <- table(spe$Banksy_celltypes,spe$Visium_domain) %>% 
  as.data.frame.matrix() %>% t() %>% 
  sweep(MARGIN = 2,STATS = colSums(.),FUN = "/") * 100

pdf(file = here("plots","06_label_transfer",
                "Visium_label_transfer","Pct_Banksy_CellTypes_domains.pdf"),
    height = 8, width = 12)
pheatmap::pheatmap(mat = pct_mat)
dev.off()


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
#   [1] grid      stats4    stats     graphics  grDevices datasets  utils    
# [8] methods   base     
# 
# other attached packages:
#   [1] dplyr_1.1.4                 here_1.0.1                 
# [3] escheR_1.8.0                ggplot2_3.5.2              
# [5] ComplexHeatmap_2.24.0       SpatialExperiment_1.18.1   
# [7] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [9] Biobase_2.68.0              GenomicRanges_1.60.0       
# [11] GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [13] S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [15] generics_0.1.4              MatrixGenerics_1.20.0      
# [17] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6            circlize_0.4.16         shape_1.4.6.1          
# [4] rjson_0.2.23            GlobalOptions_0.1.2     lattice_0.22-7         
# [7] vctrs_0.6.5             tools_4.5.0             parallel_4.5.0         
# [10] tibble_3.2.1            cluster_2.1.8.1         pkgconfig_2.0.3        
# [13] pheatmap_1.0.12         Matrix_1.7-3            RColorBrewer_1.1-3     
# [16] lifecycle_1.0.4         GenomeInfoDbData_1.2.14 compiler_4.5.0         
# [19] farver_2.1.2            textshaping_1.0.1       codetools_0.2-20       
# [22] clue_0.3-66             pillar_1.10.2           crayon_1.5.3           
# [25] DelayedArray_0.34.1     magick_2.8.6            iterators_1.0.14       
# [28] abind_1.4-8             foreach_1.5.2           tidyselect_1.2.1       
# [31] digest_0.6.37           labeling_0.4.3          rprojroot_2.0.4        
# [34] colorspace_2.1-1        cli_3.6.5               SparseArray_1.8.0      
# [37] magrittr_2.0.3          S4Arrays_1.8.0          dichromat_2.0-0.1      
# [40] withr_3.0.2             scales_1.4.0            UCSC.utils_1.4.0       
# [43] XVector_0.48.0          httr_1.4.7              ragg_1.4.0             
# [46] png_0.1-8               GetoptLong_1.0.5        doParallel_1.0.17      
# [49] rlang_1.1.6             Rcpp_1.0.14             glue_1.8.0             
# [52] jsonlite_2.0.0          R6_2.6.1                systemfonts_1.2.3    
