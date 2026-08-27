# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(ggplot2)
library(escheR)
library(here)

spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))


#Subset object to have the specific area. 
spe <- spe[,spe$Sample == "Br6660_Nac10_4080"]
coords <- spatialCoords(spe)
condition <-  coords[, "x_final"] <= 13000 & coords[, "x_final"] >= 8000 & coords[, "y_final"] <= 10000
spe <- spe[, condition]


#Read in the colors
CellType_cols <- readRDS(here("processed-data","05_Clustering","CellType_cols.Rds"))
names(CellType_cols)[18] <- "Astrocyte_Oligo"
names(CellType_cols)[22] <- "Microglia_Oligo"
CellType_cols <- CellType_cols[-c(3:5)]
CellType_cols <- c(CellType_cols,"orange")
names(CellType_cols)[20] <- "WM"

p <- make_escheR(spe) |>
  add_fill("CellTypes") +
  scale_fill_manual(values = CellType_cols)
ggsave(plot = p,filename = here("plots","05_clustering","coord_subsetted_xenium_celltypes.png"),
       height = 12, width = 12, dpi = 300)

spe$CellTypes <- ifelse(spe$CellTypes %in% c("D1_Island_A","D1_Island_B"),
                        spe$CellTypes,
                        "Other")

p <- make_escheR(spe) |>
  add_fill("CellTypes") +
  scale_fill_manual(values = CellType_cols)
ggsave(plot = p,filename = here("plots","05_clustering","coord_subsetted_xenium_IslandsOnly.png"),
       height = 12, width = 12, dpi = 300)



spe$VIP <- assay(spe,"nucleus_normcounts")["VIP",]
spe$PROK2 <- assay(spe,"nucleus_normcounts")["PROK2",]
spe$SEMA5B <- assay(spe,"nucleus_normcounts")["SEMA5B",]
spe$RXFP1 <- assay(spe,"nucleus_normcounts")["RXFP1",]


VIP <- make_escheR(spe) |>
  add_fill("VIP") +
  scale_fill_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = VIP,filename = here("plots","05_clustering","coord_subsetted_xenium_VIP_expression.png"),
       height = 12, width = 12, dpi = 300)


PROK2 <- make_escheR(spe) |>
  add_fill("PROK2") +
  scale_fill_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = PROK2,filename = here("plots","05_clustering","coord_subsetted_xenium_PROK2_expression.png"),
       height = 12, width = 12, dpi = 300)

SEMA5B <- make_escheR(spe) |>
  add_fill("SEMA5B") +
  scale_fill_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = SEMA5B,filename = here("plots","05_clustering","coord_subsetted_xenium_SEMA5B_expression.png"),
       height = 12, width = 12, dpi = 300)

RXFP1 <- make_escheR(spe) |>
  add_fill("RXFP1") +
  scale_fill_gradientn(colours = c("lightgrey","orange","red"))
ggsave(plot = RXFP1,filename = here("plots","05_clustering","coord_subsetted_xenium_RXFP1_expression.png"),
       height = 12, width = 12, dpi = 300)



legend_only <- cowplot::get_legend(RXFP1)


ggsave(here("plots","05_clustering","RXFP1_legend.pdf"),
       plot = legend_only, width = 2, height = 3)


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
#   [1] here_1.0.1                  escheR_1.8.0               
# [3] ggplot2_3.5.2               SpatialExperiment_1.18.1   
# [5] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [7] Biobase_2.68.0              GenomicRanges_1.60.0       
# [9] GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [11] S4Vectors_0.48.1            BiocGenerics_0.54.0        
# [13] generics_0.1.4              MatrixGenerics_1.20.0      
# [15] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] SparseArray_1.8.0       lattice_0.22-7          magrittr_2.0.3         
# [4] grid_4.5.0              RColorBrewer_1.1-3      rprojroot_2.0.4        
# [7] jsonlite_2.0.0          Matrix_1.7-3            httr_1.4.7             
# [10] UCSC.utils_1.4.0        viridisLite_0.4.2       scales_1.4.0           
# [13] textshaping_1.0.1       abind_1.4-8             cli_3.6.5              
# [16] rlang_1.1.6             crayon_1.5.3            XVector_0.48.0         
# [19] cowplot_1.1.3           withr_3.0.2             DelayedArray_0.34.1    
# [22] S4Arrays_1.8.0          tools_4.5.0             dplyr_1.1.4            
# [25] GenomeInfoDbData_1.2.14 vctrs_0.6.5             R6_2.6.1               
# [28] lifecycle_1.0.4         magick_2.8.6            ragg_1.4.0             
# [31] pkgconfig_2.0.3         pillar_1.10.2           gtable_0.3.6           
# [34] glue_1.8.0              Rcpp_1.1.1-1.1          systemfonts_1.2.3      
# [37] tibble_3.2.1            tidyselect_1.2.1        dichromat_2.0-0.1      
# [40] farver_2.1.2            rjson_0.2.23            labeling_0.4.3         
# [43] compiler_4.5.0   
