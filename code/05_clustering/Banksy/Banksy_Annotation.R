##Goal: Investigate Banksy louvain clustering. 
#lambda 0.8 + k_geom = 15
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
library(scDotPlot)
library(ggplot2)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

#Subset for donor Br6660
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

#Read in the louvain clustering that was generatd with a resolution value of 1.0 
banksy_louvain_clusters <- read.csv(here("processed-data","05_Clustering","Banksy_Results",
                                         "NonSpatial","Lambda0_res0.5_louvain.csv"))

head(banksy_louvain_clusters)
# X V1                V2
# 1 1  1 Br6660_NAc1_580_1
# 2 2  1 Br6660_NAc1_580_2
# 3 3  1 Br6660_NAc1_580_3
# 4 4  2 Br6660_NAc1_580_4
# 5 5  1 Br6660_NAc1_580_5
# 6 6  3 Br6660_NAc1_580_6

banksy_louvain_clusters <- banksy_louvain_clusters[,-1]

colnames(banksy_louvain_clusters) <- c("Cluster","Key")

#How many clusters? How many cells in each cluster
table(banksy_louvain_clusters$Cluster)
# 1      2      3      4      5      6      7      8      9     10     11 
# 70156 169703 244702 235062 165022 131815 170163 165626  99478 130626 128682 
# 12     13     14     15     16     17     18     19     20     21     22 
# 202949  99490 117740  60232  31380  18021  77964  33911  18962  12110  42047 


identical(banksy_louvain_clusters$Key,rownames(colData(spe)))
#[1] TRUE

identical(colnames(spe),banksy_louvain_clusters$Key)
#[1] TRUE

#Add the clusters to the object
spe$Banksy_louvain_res0.5 <- as.factor(banksy_louvain_clusters$Cluster)

spe$Sample <- factor(spe$Sample,levels = unique(spe$Sample)[1:11])

#Load in annotation_df
anno_df <- read.csv(here("processed-data","05_Clustering","nonspatial_lambda0_res0.5_annotation.csv"))

spe$CellType <- anno_df$Annotation[match(spe$Banksy_louvain_res0.5,anno_df$Cluster)]

#Make colors for CellType
CellType_cols <- Polychrome::createPalette(length(unique(spe$CellType)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(CellType_cols) <- unique(spe$CellType)
CellType_cols
# Excitatory  Microglia_A         WM_A         WM_B         WM_C          OPC
# "#DB1C5F"    "#0D87E4"    "#FDCB22"    "#F3700D"    "#385A51"    "#00FF0D"
# Astro_A Fibroblast_B Fibroblast_C     DRD2_MSN     DRD1_MSN      Astro_B
# "#FD00FD"    "#FFA7E2"    "#2AFECA"    "#7AAA16"    "#9400FF"    "#823526"
# MSN_Oligo    Inh_PVALB  D1_Island_A  Microglia_B      Inh_SST         WM_D
# "#F5DEC0"    "#93E5FF"    "#7A1699"    "#FF0DBC"    "#C4B3FB"    "#F80D2A"
# WM_E         CHAT    Ependymal  D1_Island_B
# "#224B82"    "#FBA475"    "#73690D"    "#B62A7A"

names(CellType_cols)[c(6,15,19,22)] <- c("D1_Island_B","OPC","D1_Island_A","WM_E")

#Fix the Fibroblast names as they are B,C and not A,B
spe[,spe$CellType == "Fibroblast_B"]$CellType <- "Fibroblast_A"
spe[,spe$CellType == "Fibroblast_C"]$CellType <- "Fibroblast_B"


names(CellType_cols)[c(8,9)] <- c("Fibroblast_A","Fibroblast_B")


#Save the colors
saveRDS(object = CellType_cols,file = here("processed-data","05_Clustering","CellType_cols.Rds"))


#Bargraph of celltypes by sample
celltype_sample <- as.data.frame.matrix(table(spe$CellType,spe$Sample))
celltype_sample_prop <- sweep(x = celltype_sample,MARGIN = 2,STATS = colSums(celltype_sample),FUN = "/") * 100
celltype_sample_prop$CellType <- rownames(celltype_sample_prop)
#Using CellType as id variables
celltype_sample_prop <- reshape2::melt(celltype_sample_prop)
celltype_bar <- ggplot(data = celltype_sample_prop,aes(x = variable,y=value,fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = CellType_cols) +
  labs(x = "Sample",
       y = "% of Sample") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(plot = celltype_bar,filename = here("plots","05_clustering","Banksy","celltype_bar.pdf"))

##########
#plot the clusters on the tissue
for(i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[,spe$Sample == i]
  p <- make_escheR(sub_spe) |>
    add_fill("CellType") +
    scale_fill_manual(values = CellType_cols)
  ggsave(filename = here("plots","05_clustering","Banksy","Post_Annotation",paste0(i,".png")),
         height = 20, width = 20)
}


#Save the object
saveRDS(spe, here("processed-data", "05_Clustering", "SPEs", "spe_annotated.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2025-10-30 17:00:40 EDT"
# user   system  elapsed 
# 636.358   43.219 3024.683 
# ─ Session info ──────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.0 Patched (2025-05-21 r88220)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2025-10-30
# pandoc   3.7.0.1 @ /jhpce/shared/community/core/conda_R/4.5/bin/pandoc
# quarto   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-8   2024-09-12 [2] CRAN (R 4.5.0)
# ape                    5.8-1   2024-12-16 [1] CRAN (R 4.5.0)
# aplot                  0.2.7   2025-06-21 [1] CRAN (R 4.5.0)
# beachmat               2.24.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# beeswarm               0.4.0   2021-06-01 [2] CRAN (R 4.5.0)
# Biobase              * 2.68.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocGenerics         * 0.54.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocNeighbors          2.2.0   2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocParallel           1.42.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocSingular           1.24.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# cli                    3.6.5   2025-04-23 [2] CRAN (R 4.5.0)
# cluster                2.1.8.1 2025-03-12 [3] CRAN (R 4.5.0)
# codetools              0.2-20  2024-03-31 [3] CRAN (R 4.5.0)
# colorspace             2.1-1   2024-07-26 [2] CRAN (R 4.5.0)
# cowplot                1.1.3   2024-01-22 [2] CRAN (R 4.5.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.5.0)
# data.table             1.17.2  2025-05-12 [2] CRAN (R 4.5.0)
# DelayedArray           0.34.1  2025-04-17 [2] Bioconductor 3.21 (R 4.5.0)
# deldir                 2.0-4   2024-02-28 [2] CRAN (R 4.5.0)
# dichromat              2.0-0.1 2022-05-02 [2] CRAN (R 4.5.0)
# digest                 0.6.37  2024-08-19 [2] CRAN (R 4.5.0)
# dotCall64              1.2     2024-10-04 [2] CRAN (R 4.5.0)
# dplyr                  1.1.4   2023-11-17 [2] CRAN (R 4.5.0)
# escheR               * 1.8.0   2025-04-15 [1] Bioconductor 3.21 (R 4.5.0)
# farver                 2.1.2   2024-05-13 [2] CRAN (R 4.5.0)
# fastDummies            1.7.5   2025-01-20 [2] CRAN (R 4.5.0)
# fastmap                1.2.0   2024-05-15 [2] CRAN (R 4.5.0)
# fitdistrplus           1.2-2   2025-01-07 [2] CRAN (R 4.5.0)
# fs                     1.6.6   2025-04-12 [2] CRAN (R 4.5.0)
# future                 1.49.0  2025-05-09 [2] CRAN (R 4.5.0)
# future.apply           1.11.3  2024-10-27 [2] CRAN (R 4.5.0)
# generics             * 0.1.4   2025-05-09 [2] CRAN (R 4.5.0)
# GenomeInfoDb         * 1.44.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# GenomeInfoDbData       1.2.14  2025-05-21 [2] Bioconductor
# GenomicRanges        * 1.60.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# ggbeeswarm             0.7.2   2023-04-29 [2] CRAN (R 4.5.0)
# ggfun                  0.1.8   2024-12-03 [2] CRAN (R 4.5.0)
# ggplot2              * 3.5.2   2025-04-09 [2] CRAN (R 4.5.0)
# ggplotify              0.1.2   2023-08-09 [1] CRAN (R 4.5.0)
# ggrepel                0.9.6   2024-09-07 [2] CRAN (R 4.5.0)
# ggridges               0.5.6   2024-01-23 [2] CRAN (R 4.5.0)
# ggsci                  3.2.0   2024-06-18 [2] CRAN (R 4.5.0)
# ggtree                 3.16.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.0)
# globals                0.18.0  2025-05-08 [2] CRAN (R 4.5.0)
# glue                   1.8.0   2024-09-30 [2] CRAN (R 4.5.0)
# goftest                1.2-3   2021-10-07 [2] CRAN (R 4.5.0)
# gridExtra              2.3     2017-09-09 [2] CRAN (R 4.5.0)
# gridGraphics           0.5-1   2020-12-13 [1] CRAN (R 4.5.0)
# gtable                 0.3.6   2024-10-25 [2] CRAN (R 4.5.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.5.0)
# htmltools              0.5.8.1 2024-04-04 [2] CRAN (R 4.5.0)
# htmlwidgets            1.6.4   2023-12-06 [2] CRAN (R 4.5.0)
# httpuv                 1.6.16  2025-04-16 [2] CRAN (R 4.5.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.5.0)
# ica                    1.0-3   2022-07-08 [2] CRAN (R 4.5.0)
# igraph                 2.1.4   2025-01-23 [2] CRAN (R 4.5.0)
# IRanges              * 2.42.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# irlba                  2.3.5.1 2022-10-03 [2] CRAN (R 4.5.0)
# jsonlite               2.0.0   2025-03-27 [2] CRAN (R 4.5.0)
# KernSmooth             2.23-26 2025-01-01 [3] CRAN (R 4.5.0)
# labeling               0.4.3   2023-08-29 [2] CRAN (R 4.5.0)
# later                  1.4.2   2025-04-08 [2] CRAN (R 4.5.0)
# lattice                0.22-7  2025-04-02 [3] CRAN (R 4.5.0)
# lazyeval               0.2.2   2019-03-15 [2] CRAN (R 4.5.0)
# lifecycle              1.0.4   2023-11-07 [2] CRAN (R 4.5.0)
# listenv                0.9.1   2024-01-29 [2] CRAN (R 4.5.0)
# lmtest                 0.9-40  2022-03-21 [2] CRAN (R 4.5.0)
# magick                 2.8.6   2025-03-23 [2] CRAN (R 4.5.0)
# magrittr               2.0.3   2022-03-30 [2] CRAN (R 4.5.0)
# MASS                   7.3-65  2025-02-28 [3] CRAN (R 4.5.0)
# Matrix                 1.7-3   2025-03-11 [3] CRAN (R 4.5.0)
# MatrixGenerics       * 1.20.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# matrixStats          * 1.5.0   2025-01-07 [2] CRAN (R 4.5.0)
# mime                   0.13    2025-03-17 [2] CRAN (R 4.5.0)
# miniUI                 0.1.2   2025-04-17 [2] CRAN (R 4.5.0)
# nlme                   3.1-168 2025-03-31 [3] CRAN (R 4.5.0)
# parallelly             1.44.0  2025-05-07 [2] CRAN (R 4.5.0)
# patchwork              1.3.0   2024-09-16 [2] CRAN (R 4.5.0)
# pbapply                1.7-2   2023-06-27 [2] CRAN (R 4.5.0)
# pillar                 1.10.2  2025-04-05 [2] CRAN (R 4.5.0)
# pkgconfig              2.0.3   2019-09-22 [2] CRAN (R 4.5.0)
# plotly                 4.10.4  2024-01-13 [2] CRAN (R 4.5.0)
# plyr                   1.8.9   2023-10-02 [2] CRAN (R 4.5.0)
# png                    0.1-8   2022-11-29 [2] CRAN (R 4.5.0)
# Polychrome             1.5.4   2025-04-06 [1] CRAN (R 4.5.0)
# polyclip               1.10-7  2024-07-23 [2] CRAN (R 4.5.0)
# progressr              0.15.1  2024-11-22 [2] CRAN (R 4.5.0)
# promises               1.3.2   2024-11-28 [2] CRAN (R 4.5.0)
# purrr                  1.0.4   2025-02-05 [2] CRAN (R 4.5.0)
# R6                     2.6.1   2025-02-15 [2] CRAN (R 4.5.0)
# ragg                   1.4.0   2025-04-10 [2] CRAN (R 4.5.0)
# RANN                   2.6.2   2024-08-25 [2] CRAN (R 4.5.0)
# RColorBrewer           1.1-3   2022-04-03 [2] CRAN (R 4.5.0)
# Rcpp                   1.0.14  2025-01-12 [2] CRAN (R 4.5.0)
# RcppAnnoy              0.0.22  2024-01-23 [2] CRAN (R 4.5.0)
# RcppHNSW               0.6.0   2024-02-04 [2] CRAN (R 4.5.0)
# reshape2               1.4.4   2020-04-09 [2] CRAN (R 4.5.0)
# reticulate             1.42.0  2025-03-25 [2] CRAN (R 4.5.0)
# rjson                  0.2.23  2024-09-16 [2] CRAN (R 4.5.0)
# rlang                  1.1.6   2025-04-11 [2] CRAN (R 4.5.0)
# ROCR                   1.0-11  2020-05-02 [2] CRAN (R 4.5.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.5.0)
# RSpectra               0.16-2  2024-07-18 [2] CRAN (R 4.5.0)
# rsvd                   1.0.5   2021-04-16 [2] CRAN (R 4.5.0)
# Rtsne                  0.17    2023-12-07 [2] CRAN (R 4.5.0)
# S4Arrays               1.8.0   2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# S4Vectors            * 0.46.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# ScaledMatrix           1.16.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# scales                 1.4.0   2025-04-24 [2] CRAN (R 4.5.0)
# scater                 1.36.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# scattermore            1.2     2023-06-12 [2] CRAN (R 4.5.0)
# scatterplot3d          0.3-44  2023-05-05 [1] CRAN (R 4.5.0)
# scDotPlot            * 1.2.1   2025-07-13 [1] Bioconductor 3.21 (R 4.5.0)
# sctransform            0.4.2   2025-04-30 [2] CRAN (R 4.5.0)
# scuttle                1.18.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# sessioninfo          * 1.2.3   2025-02-05 [2] CRAN (R 4.5.0)
# Seurat                 5.3.0   2025-04-23 [2] CRAN (R 4.5.0)
# SeuratObject           5.1.0   2025-04-22 [2] CRAN (R 4.5.0)
# shiny                  1.10.0  2024-12-14 [2] CRAN (R 4.5.0)
# SingleCellExperiment * 1.30.1  2025-05-07 [2] Bioconductor 3.21 (R 4.5.0)
# sp                     2.2-0   2025-02-01 [2] CRAN (R 4.5.0)
# spam                   2.11-1  2025-01-20 [2] CRAN (R 4.5.0)
# SparseArray            1.8.0   2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# SpatialExperiment    * 1.18.1  2025-05-11 [2] Bioconductor 3.21 (R 4.5.0)
# spatstat.data          3.1-6   2025-03-17 [2] CRAN (R 4.5.0)
# spatstat.explore       3.4-3   2025-05-21 [2] CRAN (R 4.5.0)
# spatstat.geom          3.4-1   2025-05-20 [2] CRAN (R 4.5.0)
# spatstat.random        3.4-1   2025-05-20 [2] CRAN (R 4.5.0)
# spatstat.sparse        3.1-0   2024-06-21 [2] CRAN (R 4.5.0)
# spatstat.univar        3.1-3   2025-05-08 [2] CRAN (R 4.5.0)
# spatstat.utils         3.1-4   2025-05-15 [2] CRAN (R 4.5.0)
# stringi                1.8.7   2025-03-27 [2] CRAN (R 4.5.0)
# stringr                1.5.1   2023-11-14 [2] CRAN (R 4.5.0)
# SummarizedExperiment * 1.38.1  2025-04-30 [2] Bioconductor 3.21 (R 4.5.0)
# survival               3.8-3   2024-12-17 [3] CRAN (R 4.5.0)
# systemfonts            1.2.3   2025-04-30 [2] CRAN (R 4.5.0)
# tensor                 1.5     2012-05-05 [2] CRAN (R 4.5.0)
# textshaping            1.0.1   2025-05-01 [2] CRAN (R 4.5.0)
# tibble                 3.2.1   2023-03-20 [2] CRAN (R 4.5.0)
# tidyr                  1.3.1   2024-01-24 [2] CRAN (R 4.5.0)
# tidyselect             1.2.1   2024-03-11 [2] CRAN (R 4.5.0)
# tidytree               0.4.6   2023-12-12 [1] CRAN (R 4.5.0)
# treeio                 1.32.0  2025-04-15 [1] Bioconductor 3.21 (R 4.5.0)
# UCSC.utils             1.4.0   2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# uwot                   0.2.3   2025-02-24 [2] CRAN (R 4.5.0)
# vctrs                  0.6.5   2023-12-01 [2] CRAN (R 4.5.0)
# vipor                  0.4.7   2023-12-18 [2] CRAN (R 4.5.0)
# viridis                0.6.5   2024-01-29 [2] CRAN (R 4.5.0)
# viridisLite            0.4.2   2023-05-02 [2] CRAN (R 4.5.0)
# withr                  3.0.2   2024-10-28 [2] CRAN (R 4.5.0)
# xtable                 1.8-4   2019-04-21 [2] CRAN (R 4.5.0)
# XVector                0.48.0  2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# yulab.utils            0.2.0   2025-01-29 [2] CRAN (R 4.5.0)
# zoo                    1.8-14  2025-04-10 [2] CRAN (R 4.5.0)
# 
# [1] /users/rphillip/R/4.5
# [2] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/library
# * ── Packages attached to the search path.
# 
# ─────────────────────────────────────────────────────────────────────────────────────────────────────────────
# 
# 
# 
# 
