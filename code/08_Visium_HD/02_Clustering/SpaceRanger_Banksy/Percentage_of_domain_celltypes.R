#cat  cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(RColorBrewer)
library(HDF5Array)
library(circlize)
library(here)

#############################################
###########     Human   #####################
#############################################

###### Visium-HD
# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")

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

#What is the cell type composition of each domain? 
dc <- prop.table(table(sfe$snRNA_label,sfe$Spatial_Domain),margin = 2) * 100
cts  <- rownames(dc)
doms <- colnames(dc)

# distinct categorical palettes (hcl.colors is base R, no extra deps)

#Load cluster colors
#Plot the max type
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE) 
# Loading objects:
#   cluster_cols

ct_cols  <- cluster_cols[-14]
dom_cols <- readRDS(here("processed-data", "HD_Full_Analysis", "Cluster_colors",
                         "clust_M1_lam0.8_k50_res0.4_cell_level_spatial_colors_sr.Rds"))
anno_df$spatial_0.4 <- as.character(anno_df$spatial_0.4)
names(dom_cols) <- anno_df$Annotation


# continuous fill for the percentages
col_fun <- colorRamp2(
  c(0, max(dc) / 2, max(dc)),
  c("white", "orange", "red")
)

# annotation bars
row_ha <- rowAnnotation(
  `Cell type` = cts,
  col = list(`Cell type` = ct_cols),
  show_annotation_name = FALSE
)

col_ha <- HeatmapAnnotation(
  Domain = doms,
  col = list(Domain = dom_cols),
  show_annotation_name = FALSE
)

hm <- Heatmap(
  dc,
  rect_gp = gpar(col = "black", lwd = 0.5),
  name             = "% of domain",
  col              = col_fun,
  left_annotation  = row_ha,
  top_annotation   = col_ha,
  cluster_rows     = TRUE,
  cluster_columns  = TRUE,
  row_names_side   = "right",
  column_names_rot = 45,
  heatmap_legend_param = list(title = "% of domain"),
  column_title = "Cell type composition per spatial domain"
)

pdf(file = here("plots","HD_Full_Analysis",
                "Banksy_sr","Annotated_snRNA-seq_spatial_domain_heatmap.pdf"),height = 12, width = 16)
draw(hm)
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
#   [1] here_1.0.1                      circlize_0.4.16                
# [3] HDF5Array_1.36.0                h5mread_1.0.0                  
# [5] rhdf5_2.52.0                    DelayedArray_0.34.1            
# [7] SparseArray_1.8.0               S4Arrays_1.8.0                 
# [9] abind_1.4-8                     Matrix_1.7-3                   
# [11] RColorBrewer_1.1-3              ComplexHeatmap_2.24.0          
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
# [7] magrittr_2.0.3            clue_0.3-66              
# [9] GetoptLong_1.0.5          multcomp_1.4-28          
# [11] e1071_1.7-16              compiler_4.5.0           
# [13] DelayedMatrixStats_1.30.0 png_0.1-8                
# [15] sfheaders_0.4.4           fftwtools_0.9-11         
# [17] SpatialExperiment_1.18.1  wk_0.9.4                 
# [19] shape_1.4.6.1             crayon_1.5.3             
# [21] fastmap_1.2.0             magick_2.8.6             
# [23] XVector_0.48.0            scuttle_1.18.0           
# [25] UCSC.utils_1.4.0          beachmat_2.24.0          
# [27] jsonlite_2.0.0            rhdf5filters_1.20.0      
# [29] Rhdf5lib_1.30.0           BiocParallel_1.42.0      
# [31] jpeg_0.1-11               tiff_0.1-12              
# [33] terra_1.8-50              cluster_2.1.8.1          
# [35] parallel_4.5.0            LearnBayes_2.15.1        
# [37] R6_2.6.1                  limma_3.64.3             
# [39] boot_1.3-31               Rcpp_1.1.1-1.1           
# [41] iterators_1.0.14          zoo_1.8-14               
# [43] R.utils_2.13.0            splines_4.5.0            
# [45] EBImage_4.50.0            doParallel_1.0.17        
# [47] codetools_0.2-20          lattice_0.22-7           
# [49] coda_0.19-4.1             survival_3.8-3           
# [51] sf_1.0-21                 units_0.8-7              
# [53] spData_2.3.4              proxy_0.4-27             
# [55] KernSmooth_2.23-26        foreach_1.5.2            
# [57] rprojroot_2.0.4           sp_2.2-0                 
# [59] RCurl_1.98-1.17           sparseMatrixStats_1.20.0 
# [61] class_7.3-23              spatialreg_1.3-6         
# [63] tools_4.5.0               BiocNeighbors_2.2.0      
# [65] data.table_1.17.2         locfit_1.5-9.12          
# [67] mvtnorm_1.3-3             Cairo_1.6-2              
# [69] spdep_1.3-11              DropletUtils_1.28.0      
# [71] edgeR_4.6.2               colorspace_2.1-1         
# [73] nlme_3.1-168              GenomeInfoDbData_1.2.14  
# [75] cli_3.6.5                 R.methodsS3_1.8.2        
# [77] zeallot_0.1.0             digest_0.6.37            
# [79] classInt_0.4-11           dqrng_0.4.1              
# [81] TH.data_1.1-3             rjson_0.2.23             
# [83] htmlwidgets_1.6.4         htmltools_0.5.8.1        
# [85] R.oo_1.27.1               lifecycle_1.0.4          
# [87] httr_1.4.7                GlobalOptions_0.1.2      
# [89] statmod_1.5.0             MASS_7.3-65              
