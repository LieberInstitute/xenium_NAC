# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(tidyr)
library(dplyr)
library(here)

# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")

sfe <- loadHDF5SummarizedExperiment(sfe_dir)

#Make a bargraoh of cell type by spot_class
x <- as.data.frame(table(sfe$snRNA_label,sfe$rctd_spot_class))
colnames(x) <- c("Cell_Type","spot_class","Freq")
#Make a bargraph of number of cells by 
p1 <- ggplot(data = x,aes(x = Cell_Type, y = Freq,fill = spot_class))+
  geom_bar(stat = "identity",position = "dodge") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = p1,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "CellType_SpotClass_Bargraph.pdf"),
       height = 8, width = 12)


# pull the relevant columns
meta <- as.data.frame(
  colData(sfe)[, c("rctd_first_type", "rctd_second_type", "rctd_spot_class")]
)

# keep only doublets — first/second pairing is only meaningful here
meta_db <- meta %>%
  filter(rctd_spot_class %in% c("doublet_certain", "doublet_uncertain"))

# When there are doublets, what is the predicted second cell type? 
# I anticipate that doublets will be cells in which neurons are surrounded by many astros or oligos
# This will caause a shared transcriptional signature 
# count pairs, fill empty combos with 0, then convert to % within each first_type
pct <- meta_db %>%
  count(rctd_first_type, rctd_second_type, name = "n") %>%
  complete(rctd_first_type, rctd_second_type, fill = list(n = 0)) %>%
  group_by(rctd_first_type) %>%
  mutate(total = sum(n),
         pct   = ifelse(total > 0, 100 * n / total, 0)) %>%
  ungroup()

# plot
p2 <- ggplot(pct, aes(x = rctd_first_type, y = rctd_second_type, fill = pct)) +
  geom_tile(color = "grey90") +
  geom_text(aes(label = ifelse(pct > 0, round(pct), "")), size = 2.5) +
  scale_fill_gradientn(colours = c("white","lightgrey","orange","red") )+
  labs(x = "First type", y = "Second type", fill = "% of first type") +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p2,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "RCTD_First_Second_Type_Comparison_Dou.pdf"),
       height = 8, width = 12)



# Now for the doublets, does the SingleR call agree with the first_type call from RCTD? 
# This will help determine if we just go ahead and move forward with the RCTD calls. 
# pull the relevant columns
meta <- as.data.frame(
  colData(sfe)[, c("rctd_first_type", "rctd_second_type","pruned.labels","rctd_spot_class")]
)

# keep only doublets — first/second pairing is only meaningful here
meta_db <- meta %>%
  filter(rctd_spot_class %in% c("doublet_certain", "doublet_uncertain"))

# put both calls on the SAME factor levels so the diagonal = agreement
lvls <- sort(union(as.character(meta_db$rctd_first_type),
                   as.character(meta_db$pruned.labels)))
meta_db$rctd_first_type   <- factor(meta_db$rctd_first_type,   levels = lvls)
meta_db$pruned.labels <- factor(meta_db$pruned.labels, levels = lvls)

# confusion matrix -> proportion within each first_type (columns sum to 1)
prop <- meta_db %>%
  count(rctd_first_type, pruned.labels, name = "n") %>%
  complete(rctd_first_type, pruned.labels, fill = list(n = 0)) %>%
  group_by(rctd_first_type) %>%
  mutate(total = sum(n),
         prop  = ifelse(total > 0, n / total, 0)) %>%
  ungroup()

p3 <- ggplot(prop, aes(x = rctd_first_type, y = pruned.labels, fill = prop)) +
  geom_tile(color = "black") +
  # uncomment to label tiles:
  geom_text(aes(label = ifelse(prop > 0.01, round(prop, 2), "")), size = 2.5) +
  scale_fill_gradientn(colours = c("white","lightgrey","orange","red")) +
  labs(x = "First type (doublet mode)",
       y = "Singlet call",
       fill = "Proportion\nof first type") +
  coord_equal() +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(plot = p3,filename = here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                                 "RCTD_SingleR_Comparison_Doublets_Only.pdf"),
       height = 8, width = 12)


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
#   [1] here_1.0.1                      dplyr_1.1.4                    
# [3] tidyr_1.3.1                     escheR_1.8.0                   
# [5] ggplot2_3.5.2                   HDF5Array_1.36.0               
# [7] h5mread_1.0.0                   rhdf5_2.52.0                   
# [9] DelayedArray_0.34.1             SparseArray_1.8.0              
# [11] S4Arrays_1.8.0                  abind_1.4-8                    
# [13] Matrix_1.7-3                    sessioninfo_1.2.3              
# [15] SpatialExperiment_1.18.1        SingleCellExperiment_1.30.1    
# [17] SummarizedExperiment_1.38.1     Biobase_2.68.0                 
# [19] GenomicRanges_1.60.0            GenomeInfoDb_1.44.0            
# [21] IRanges_2.42.0                  S4Vectors_0.46.0               
# [23] BiocGenerics_0.54.0             generics_0.1.4                 
# [25] MatrixGenerics_1.20.0           matrixStats_1.5.0              
# [27] SpatialFeatureExperiment_1.10.1
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
# [23] labeling_0.4.3            scuttle_1.18.0           
# [25] UCSC.utils_1.4.0          ragg_1.4.0               
# [27] purrr_1.0.4               beachmat_2.24.0          
# [29] jsonlite_2.0.0            rhdf5filters_1.20.0      
# [31] Rhdf5lib_1.30.0           BiocParallel_1.42.0      
# [33] jpeg_0.1-11               tiff_0.1-12              
# [35] terra_1.8-50              parallel_4.5.0           
# [37] LearnBayes_2.15.1         R6_2.6.1                 
# [39] RColorBrewer_1.1-3        limma_3.64.3             
# [41] boot_1.3-31               Rcpp_1.0.14              
# [43] zoo_1.8-14                R.utils_2.13.0           
# [45] tidyselect_1.2.1          splines_4.5.0            
# [47] dichromat_2.0-0.1         EBImage_4.50.0           
# [49] codetools_0.2-20          lattice_0.22-7           
# [51] tibble_3.2.1              withr_3.0.2              
# [53] coda_0.19-4.1             survival_3.8-3           
# [55] sf_1.0-21                 units_0.8-7              
# [57] spData_2.3.4              proxy_0.4-27             
# [59] pillar_1.10.2             KernSmooth_2.23-26       
# [61] rprojroot_2.0.4           sp_2.2-0                 
# [63] RCurl_1.98-1.17           sparseMatrixStats_1.20.0 
# [65] scales_1.4.0              class_7.3-23             
# [67] glue_1.8.0                spatialreg_1.3-6         
# [69] tools_4.5.0               BiocNeighbors_2.2.0      
# [71] data.table_1.17.2         locfit_1.5-9.12          
# [73] mvtnorm_1.3-3             grid_4.5.0               
# [75] spdep_1.3-11              DropletUtils_1.28.0      
# [77] edgeR_4.6.2               nlme_3.1-168             
# [79] GenomeInfoDbData_1.2.14   cli_3.6.5                
# [81] textshaping_1.0.1         gtable_0.3.6             
# [83] R.methodsS3_1.8.2         zeallot_0.1.0            
# [85] digest_0.6.37             classInt_0.4-11          
# [87] dqrng_0.4.1               TH.data_1.1-3            
# [89] rjson_0.2.23              htmlwidgets_1.6.4        
# [91] farver_2.1.2              htmltools_0.5.8.1        
# [93] R.oo_1.27.1               lifecycle_1.0.4          
# [95] httr_1.4.7                statmod_1.5.0            
# [97] MASS_7.3-65        
# 
# 
# 
#  
# 
# 
