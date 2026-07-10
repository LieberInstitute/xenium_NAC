# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(HDF5Array)
library(tidyverse)
library(ggplot2)
library(here)


sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_NMF")
sfe <- loadHDF5SummarizedExperiment(sfe_dir)

#Add spatial clusters and non-spatial clusters
spatial_clust <- read.csv(here("processed-data", "HD_Full_Analysis",
                               "sr_spatial_banksy_clusters_res_0.4.csv"))
rownames(spatial_clust) <- spatial_clust$V2
spatial_clust <- spatial_clust[colnames(sfe),]

#Add spatial data to sfe 
stopifnot(identical(spatial_clust$V2,colnames(sfe)))

sfe$spatial_0.4 <- spatial_clust$V1

#Annotate
anno_df <- read.csv(here("processed-data","spatial_0.4_Annotations.csv"))
anno_df[11,"Annotation"] <- "Hypo"
sfe$Spatial_Domain <- anno_df[match(sfe$spatial_0.4,anno_df$spatial_0.4),"Annotation"]

# --- Config: adjust these ---
nmf_factors <- c("nmf34", "nmf35", "nmf44") 
group_var   <- "Spatial_Domain"                     
threshold   <- 0                              

# --- Pull relevant colData into a data.frame ---
cd <- as.data.frame(colData(sfe))[, c(group_var, nmf_factors)]

# --- Center and scale each NMF factor (z-score across all spots/cells) ---
cd_scaled <- cd
cd_scaled[nmf_factors] <- scale(cd[nmf_factors],center = TRUE,scale = TRUE) #z score

# --- Reshape to long format ---
cd_long <- cd_scaled %>%
  pivot_longer(cols = all_of(nmf_factors), names_to = "Factor", values_to = "ScaledValue")

cd_raw_long <- cd %>%
  pivot_longer(cols = all_of(nmf_factors), names_to = "Factor", values_to = "RawValue")

cd_long$RawValue <- cd_raw_long$RawValue


# --- Summarize per group x factor ---
summary_df <- cd_long %>%
  group_by(.data[[group_var]], Factor) %>%
  summarise(
    MeanScaled = mean(ScaledValue, na.rm = TRUE),
    PctAbove   = mean(RawValue > threshold, na.rm = TRUE) * 100,
    .groups = "drop"
  )

# --- Dot plot ---
dotplot <- ggplot(summary_df, aes(x = Factor, y = .data[[group_var]])) +
  geom_point(aes(size = PctAbove, color = MeanScaled)) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
  scale_size(range = c(0, 8), name = "% spots > threshold") +
  labs(x = "NMF Factor", y = group_var, color = "Mean\n(scaled)") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_line(color = "grey90")
  )


ggsave(plot = dotplot,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                      "NMF_DotPlot_SpatialDomains.pdf"),
       height = 8, width = 6)

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
#   [1] here_1.0.1                      lubridate_1.9.4                
# [3] forcats_1.0.0                   stringr_1.5.1                  
# [5] dplyr_1.1.4                     purrr_1.0.4                    
# [7] readr_2.1.5                     tidyr_1.3.1                    
# [9] tibble_3.2.1                    ggplot2_3.5.2                  
# [11] tidyverse_2.0.0                 HDF5Array_1.36.0               
# [13] h5mread_1.0.0                   rhdf5_2.52.0                   
# [15] DelayedArray_0.34.1             SparseArray_1.8.0              
# [17] S4Arrays_1.8.0                  IRanges_2.42.0                 
# [19] abind_1.4-8                     S4Vectors_0.48.1               
# [21] MatrixGenerics_1.20.0           matrixStats_1.5.0              
# [23] BiocGenerics_0.54.0             generics_0.1.4                 
# [25] Matrix_1.7-3                    SpatialFeatureExperiment_1.10.1
# 
# loaded via a namespace (and not attached):
#   [1] DBI_1.2.3                   bitops_1.0-9               
# [3] deldir_2.0-4                s2_1.1.8                   
# [5] sandwich_3.1-1              rlang_1.1.6                
# [7] magrittr_2.0.3              multcomp_1.4-28            
# [9] e1071_1.7-16                compiler_4.5.0             
# [11] DelayedMatrixStats_1.30.0   systemfonts_1.2.3          
# [13] png_0.1-8                   sfheaders_0.4.4            
# [15] fftwtools_0.9-11            vctrs_0.6.5                
# [17] pkgconfig_2.0.3             SpatialExperiment_1.18.1   
# [19] wk_0.9.4                    crayon_1.5.3               
# [21] fastmap_1.2.0               magick_2.8.6               
# [23] XVector_0.48.0              labeling_0.4.3             
# [25] scuttle_1.18.0              tzdb_0.5.0                 
# [27] UCSC.utils_1.4.0            ragg_1.4.0                 
# [29] beachmat_2.24.0             GenomeInfoDb_1.44.0        
# [31] jsonlite_2.0.0              rhdf5filters_1.20.0        
# [33] Rhdf5lib_1.30.0             BiocParallel_1.42.0        
# [35] jpeg_0.1-11                 tiff_0.1-12                
# [37] terra_1.8-50                parallel_4.5.0             
# [39] LearnBayes_2.15.1           R6_2.6.1                   
# [41] stringi_1.8.7               RColorBrewer_1.1-3         
# [43] limma_3.64.3                boot_1.3-31                
# [45] GenomicRanges_1.60.0        Rcpp_1.1.1-1.1             
# [47] SummarizedExperiment_1.38.1 zoo_1.8-14                 
# [49] R.utils_2.13.0              timechange_0.3.0           
# [51] tidyselect_1.2.1            splines_4.5.0              
# [53] dichromat_2.0-0.1           EBImage_4.50.0             
# [55] codetools_0.2-20            lattice_0.22-7             
# [57] withr_3.0.2                 Biobase_2.68.0             
# [59] coda_0.19-4.1               survival_3.8-3             
# [61] sf_1.0-21                   units_0.8-7                
# [63] spData_2.3.4                proxy_0.4-27               
# [65] pillar_1.10.2               KernSmooth_2.23-26         
# [67] rprojroot_2.0.4             sp_2.2-0                   
# [69] RCurl_1.98-1.17             hms_1.1.3                  
# [71] sparseMatrixStats_1.20.0    scales_1.4.0               
# [73] class_7.3-23                glue_1.8.0                 
# [75] spatialreg_1.3-6            tools_4.5.0                
# [77] BiocNeighbors_2.2.0         data.table_1.17.2          
# [79] locfit_1.5-9.12             mvtnorm_1.3-3              
# [81] grid_4.5.0                  spdep_1.3-11               
# [83] DropletUtils_1.28.0         edgeR_4.6.2                
# [85] SingleCellExperiment_1.30.1 nlme_3.1-168               
# [87] GenomeInfoDbData_1.2.14     cli_3.6.5                  
# [89] textshaping_1.0.1           gtable_0.3.6               
# [91] R.methodsS3_1.8.2           zeallot_0.1.0              
# [93] digest_0.6.37               classInt_0.4-11            
# [95] dqrng_0.4.1                 TH.data_1.1-3              
# [97] rjson_0.2.23                htmlwidgets_1.6.4          
# [99] farver_2.1.2                htmltools_0.5.8.1          
# [101] R.oo_1.27.1                 lifecycle_1.0.4            
# [103] httr_1.4.7                  statmod_1.5.0              
# [105] MASS_7.3-65                
