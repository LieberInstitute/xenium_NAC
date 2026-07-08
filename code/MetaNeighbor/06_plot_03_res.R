# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SingleCellExperiment)
library(ComplexHeatmap)
library(MetaNeighbor)
library(RColorBrewer)
library(circlize)
library(here)

aurocs <- readRDS(here("processed-data","MetaNeighbor","All_Species_MSNs_Only_aurocs.Rds"))


#Save the heatmap to pull the row and column dendrogram 
hm_dend <- plotHeatmap(aurocs)

#Colors for plot 
aurocs_cols <- rev(colorRampPalette(brewer.pal(11, "RdYlBu"))(100))

parts <- strsplit(rownames(aurocs), "\\|")
prefix   <- sapply(parts, `[`, 1)   
celltype <- sapply(parts, `[`, 2)     

species <- sub("_.*", "", prefix)         
modality <- sub("^[^_]+_", "", prefix)     

# sanity check
table(species, modality)

# --- Color palettes for annotations ---
species_colors <- setNames(
  brewer.pal(length(unique(species)), "Dark2"),
  unique(species)
)

modality_colors <- setNames(
  brewer.pal(length(unique(modality)), "Set2"),
  unique(modality)
)

category <- c("Island","D1","D1","D2","D2",
              "Island","Island","Island","D1","Island",
              "D1","D1","D1/D2","D2","D2",
              "D2","D1","Island","D1","Island",
              "D2","D2","Island","Island","D1","D2","Island",
              "Island","D1/D2","D1/D2","D1/D2")
names(category) <- celltype

category_colors <- c(
  "Island" = "darkorchid4",
  "D1"     = "dodgerblue",
  "D2"     = "orange",
  "D1/D2"  = "seagreen"
)

ha_row <- rowAnnotation(
  Category = category,
  Species = species,
  modality = modality,
  col = list(Species = species_colors, modality = modality_colors,Category = category_colors),
  annotation_name_side = "top"
)

ha_col <- HeatmapAnnotation(
  Category = category,
  Species = species,
  modality = modality,
  col = list(Species = species_colors, modality = modality_colors,Category = category_colors),
  show_legend = FALSE
)


# Make the heatmap
hm <- Heatmap(
  aurocs,
  name = "AUROC",
  col = colorRamp2(seq(0, 1, length = 100), aurocs_cols),
  na_col = gray(0.95),
  left_annotation = ha_row,
  top_annotation = ha_col,
  cluster_rows = hm_dend$rowDendrogram,
  cluster_columns = hm_dend$colDendrogram,
  row_names_side = "right",
  column_names_rot = 90,
  column_names_gp = gpar(fontsize = 7),
  row_names_gp = gpar(fontsize = 7),
  rect_gp = gpar(col = "black", lwd = 0.5),
)


pdf(file = here("plots","MetaNeighbor","All_Species_MSNs_All_Modalities.pdf"),
    width = 16,
    height = 12)
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
#   [1] here_1.0.1                  circlize_0.4.16            
# [3] RColorBrewer_1.1-3          MetaNeighbor_1.28.0        
# [5] ComplexHeatmap_2.24.0       SingleCellExperiment_1.30.1
# [7] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [11] IRanges_2.42.0              S4Vectors_0.48.1           
# [13] BiocGenerics_0.54.0         generics_0.1.4             
# [15] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] gplots_3.2.0            SparseArray_1.8.0       bitops_1.0-9           
# [4] KernSmooth_2.23-26      shape_1.4.6.1           gtools_3.9.5           
# [7] lattice_0.22-7          digest_0.6.37           magrittr_2.0.3         
# [10] caTools_1.18.3          iterators_1.0.14        foreach_1.5.2          
# [13] doParallel_1.0.17       rprojroot_2.0.4         jsonlite_2.0.0         
# [16] Matrix_1.7-3            GlobalOptions_0.1.2     httr_1.4.7             
# [19] UCSC.utils_1.4.0        codetools_0.2-20        abind_1.4-8            
# [22] crayon_1.5.3            XVector_0.48.0          DelayedArray_0.34.1    
# [25] S4Arrays_1.8.0          tools_4.5.0             parallel_4.5.0         
# [28] colorspace_2.1-1        GenomeInfoDbData_1.2.14 GetoptLong_1.0.5       
# [31] R6_2.6.1                png_0.1-8               magick_2.8.6           
# [34] clue_0.3-66             cluster_2.1.8.1         Rcpp_1.1.1-1.1         
# [37] rjson_0.2.23            Cairo_1.6-2             compiler_4.5.0         
