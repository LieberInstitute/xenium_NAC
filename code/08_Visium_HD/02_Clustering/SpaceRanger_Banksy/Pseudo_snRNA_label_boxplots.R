# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(scuttle)
library(ggplot2)
library(here)

#Read in the filtered sfe object
#Now save.
sfe_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_snRNA_label_pseudo"
)

sfe <- loadHDF5SummarizedExperiment(sfe_dir)

sfe


#Set up colors
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE) 

ct_cols  <- cluster_cols[-14]

#Calculate log-normcalized counts
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

sfe$DRD1 <- logcounts(sfe)["DRD1",]
sfe$PDYN <- logcounts(sfe)["PDYN",]
sfe$DRD2 <- logcounts(sfe)["DRD2",]
sfe$ADORA2A <- logcounts(sfe)["ADORA2A",]
sfe$RXFP1 <- logcounts(sfe)["RXFP1",]
sfe$SEMA5B <- logcounts(sfe)["SEMA5B",]

genes <- c("DRD1","PDYN","DRD2","ADORA2A","RXFP1","SEMA5B")
# Loop through the vector
for (g in seq_along(genes)){
 print(g)
 gene_name <- genes[[g]]
  
  p <- ggplot(data = as.data.frame(colData(sfe)),aes(x = snRNA_label,y = .data[[gene_name]]),color = snRNA_label) +
    geom_boxplot(outlier.shape = NA,aes(colour = snRNA_label)) +
    geom_jitter(width = 0.2,aes(colour = snRNA_label)) +
    scale_colour_manual(values = ct_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1),
          legend.position="none")
  ggsave(plot = p,filename = here("plots","HD_Full_Analysis",
                         paste0("pseudobulk_expression_snRNA",gene_name,".pdf")),
         height = 8, width = 8)
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
#   [1] here_1.0.1                      ggplot2_3.5.2                  
# [3] scuttle_1.18.0                  HDF5Array_1.36.0               
# [5] h5mread_1.0.0                   rhdf5_2.52.0                   
# [7] DelayedArray_0.34.1             SparseArray_1.8.0              
# [9] S4Arrays_1.8.0                  abind_1.4-8                    
# [11] Matrix_1.7-3                    SpatialExperiment_1.18.1       
# [13] SingleCellExperiment_1.30.1     SummarizedExperiment_1.38.1    
# [15] Biobase_2.68.0                  GenomicRanges_1.60.0           
# [17] GenomeInfoDb_1.44.0             IRanges_2.42.0                 
# [19] S4Vectors_0.48.1                BiocGenerics_0.54.0            
# [21] generics_0.1.4                  MatrixGenerics_1.20.0          
# [23] matrixStats_1.5.0               SpatialFeatureExperiment_1.10.1
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3                 bitops_1.0-9             
# [3] deldir_2.0-4              s2_1.1.8                 
# [5] sandwich_3.1-1            rlang_1.1.6              
# [7] magrittr_2.0.3            multcomp_1.4-28          
# [9] e1071_1.7-16              compiler_4.5.0           
# [11] DelayedMatrixStats_1.30.0 systemfonts_1.2.3        
# [13] png_0.1-8                 sfheaders_0.4.4          
# [15] fftwtools_0.9-11          vctrs_0.6.5              
# [17] pkgconfig_2.0.3           wk_0.9.4                 
# [19] crayon_1.5.3              fastmap_1.2.0            
# [21] magick_2.8.6              XVector_0.48.0           
# [23] labeling_0.4.3            UCSC.utils_1.4.0         
# [25] ragg_1.4.0                beachmat_2.24.0          
# [27] jsonlite_2.0.0            rhdf5filters_1.20.0      
# [29] Rhdf5lib_1.30.0           BiocParallel_1.42.0      
# [31] jpeg_0.1-11               tiff_0.1-12              
# [33] terra_1.8-50              parallel_4.5.0           
# [35] LearnBayes_2.15.1         R6_2.6.1                 
# [37] RColorBrewer_1.1-3        limma_3.64.3             
# [39] boot_1.3-31               Rcpp_1.1.1-1.1           
# [41] zoo_1.8-14                R.utils_2.13.0           
# [43] tidyselect_1.2.1          splines_4.5.0            
# [45] dichromat_2.0-0.1         EBImage_4.50.0           
# [47] codetools_0.2-20          lattice_0.22-7           
# [49] tibble_3.2.1              withr_3.0.2              
# [51] coda_0.19-4.1             survival_3.8-3           
# [53] sf_1.0-21                 units_0.8-7              
# [55] spData_2.3.4              proxy_0.4-27             
# [57] pillar_1.10.2             KernSmooth_2.23-26       
# [59] rprojroot_2.0.4           sp_2.2-0                 
# [61] RCurl_1.98-1.17           sparseMatrixStats_1.20.0 
# [63] scales_1.4.0              class_7.3-23             
# [65] glue_1.8.0                spatialreg_1.3-6         
# [67] tools_4.5.0               BiocNeighbors_2.2.0      
# [69] data.table_1.17.2         locfit_1.5-9.12          
# [71] mvtnorm_1.3-3             grid_4.5.0               
# [73] spdep_1.3-11              DropletUtils_1.28.0      
# [75] edgeR_4.6.2               nlme_3.1-168             
# [77] GenomeInfoDbData_1.2.14   cli_3.6.5                
# [79] textshaping_1.0.1         dplyr_1.1.4              
# [81] gtable_0.3.6              R.methodsS3_1.8.2        
# [83] zeallot_0.1.0             digest_0.6.37            
# [85] classInt_0.4-11           dqrng_0.4.1              
# [87] TH.data_1.1-3             rjson_0.2.23             
# [89] htmlwidgets_1.6.4         farver_2.1.2             
# [91] htmltools_0.5.8.1         R.oo_1.27.1              
# [93] lifecycle_1.0.4           httr_1.4.7               
# [95] statmod_1.5.0             MASS_7.3-65
# 
# 
# 
# 
# 
# 
# 
# 
