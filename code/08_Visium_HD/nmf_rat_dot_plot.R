library(SingleCellExperiment)
library(HDF5Array) 
library(ggplot2)
library(tidyr)
library(dplyr)
library(here)

###### Visium-HD
# Load object 
sfe_dir <-  here("processed-data", "HD_Full_Analysis", "sfe_with_rat_NMF")
sfe <- loadHDF5SummarizedExperiment(sfe_dir)


cd <- as.data.frame(colData(sfe))

domain_col <- "Spatial_Domain"
pattern_cols <- c("nmf10_rat","nmf18_rat")
col_min      <- -2.5            
col_max      <- 2.5
dot_scale    <- 6  


# detect NMF pattern columns (nmf1, nmf2, ... or nmf1_norm, nmf2_norm, ...),
# sorted numerically rather than alphabetically (so nmf2 sorts before nmf10

# ---- Build dot plot summary stats (mirrors Seurat::DotPlot's logic) --------
# for each pattern x domain: avg_exp = mean value across spots in that domain,
# pct_exp = % of spots with a nonzero pattern value ("percent expressing"
# analog), avg_exp_scaled = per-pattern z-score across domains, clipped to
# [col_min, col_max] the same way Seurat::DotPlot clips by default.
dot_df <- cd[, c(domain_col, pattern_cols)] |>
  pivot_longer(cols = all_of(pattern_cols), names_to = "pattern", values_to = "value") |>
  group_by(.data[[domain_col]], pattern) |>
  summarise(
    avg_exp = mean(value, na.rm = TRUE),
    pct_exp = 100 * mean(value > 0, na.rm = TRUE),
    .groups = "drop"
  ) |>
  group_by(pattern) |>
  mutate(avg_exp_scaled = as.numeric(scale(avg_exp))) |>
  ungroup() |>
  mutate(avg_exp_scaled = pmin(pmax(avg_exp_scaled, col_min), col_max))

dot_df$pattern <- factor(dot_df$pattern, levels = pattern_cols)

# ---- Plot --------------------------------------------------------------------
p <- ggplot(dot_df, aes(y = pattern, x = .data[[domain_col]])) +
  geom_point(aes(size = pct_exp, color = avg_exp_scaled)) +
  scale_size(range = c(0, dot_scale), name = "% spots\n> 0") +
  scale_color_gradient(low = "lightgrey", high = "darkred", name = "Scaled\navg value") +
  labs(x = "NMF pattern", y = "Spatial domain") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


ggsave(here("plots","HD_Full_Analysis","rat_nmf_spatial_domain_dotplot.pdf"), plot = p,
       width = 10, height = 6)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# Reproducibility information:
# [1] "2026-07-17 11:37:15 EDT"
# user   system  elapsed 
# 230.851   14.397 1015.269 
# ─ Session info ──────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.5.0 Patched (2025-05-21 r88220)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2026-07-17
# pandoc   3.7.0.1 @ /jhpce/shared/community/core/conda_R/4.5/bin/pandoc
# quarto   NA
# 
# ─ Packages ──────────────────────────────────────────────────────────────────────────
# package                  * version   date (UTC) lib source
# abind                    * 1.4-8     2024-09-12 [2] CRAN (R 4.5.0)
# beachmat                   2.24.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# Biobase                  * 2.68.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocGenerics             * 0.54.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocNeighbors              2.2.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# BiocParallel               1.42.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# bitops                     1.0-9     2024-10-03 [2] CRAN (R 4.5.0)
# boot                       1.3-31    2024-08-28 [3] CRAN (R 4.5.0)
# class                      7.3-23    2025-01-01 [3] CRAN (R 4.5.0)
# classInt                   0.4-11    2025-01-08 [2] CRAN (R 4.5.0)
# cli                        3.6.5     2025-04-23 [2] CRAN (R 4.5.0)
# coda                       0.19-4.1  2024-01-31 [2] CRAN (R 4.5.0)
# codetools                  0.2-20    2024-03-31 [3] CRAN (R 4.5.0)
# crayon                     1.5.3     2024-06-20 [2] CRAN (R 4.5.0)
# data.table                 1.17.2    2025-05-12 [2] CRAN (R 4.5.0)
# DBI                        1.2.3     2024-06-02 [2] CRAN (R 4.5.0)
# DelayedArray             * 0.34.1    2025-04-17 [2] Bioconductor 3.21 (R 4.5.0)
# DelayedMatrixStats         1.30.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# deldir                     2.0-4     2024-02-28 [2] CRAN (R 4.5.0)
# dichromat                  2.0-0.1   2022-05-02 [2] CRAN (R 4.5.0)
# digest                     0.6.37    2024-08-19 [2] CRAN (R 4.5.0)
# dplyr                    * 1.1.4     2023-11-17 [2] CRAN (R 4.5.0)
# dqrng                      0.4.1     2024-05-28 [2] CRAN (R 4.5.0)
# DropletUtils               1.28.0    2025-04-17 [2] Bioconductor 3.21 (R 4.5.0)
# e1071                      1.7-16    2024-09-16 [2] CRAN (R 4.5.0)
# EBImage                    4.50.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# edgeR                      4.6.2     2025-05-07 [2] Bioconductor 3.21 (R 4.5.0)
# farver                     2.1.2     2024-05-13 [2] CRAN (R 4.5.0)
# fastmap                    1.2.0     2024-05-15 [2] CRAN (R 4.5.0)
# fftwtools                  0.9-11    2021-03-01 [2] CRAN (R 4.5.0)
# generics                 * 0.1.4     2025-05-09 [2] CRAN (R 4.5.0)
# GenomeInfoDb             * 1.44.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# GenomeInfoDbData           1.2.14    2025-05-21 [2] Bioconductor
# GenomicRanges            * 1.60.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# ggplot2                  * 3.5.2     2025-04-09 [2] CRAN (R 4.5.0)
# glue                       1.8.0     2024-09-30 [2] CRAN (R 4.5.0)
# gtable                     0.3.6     2024-10-25 [2] CRAN (R 4.5.0)
# h5mread                  * 1.0.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# HDF5Array                * 1.36.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# here                     * 1.0.1     2020-12-13 [2] CRAN (R 4.5.0)
# htmltools                  0.5.8.1   2024-04-04 [2] CRAN (R 4.5.0)
# htmlwidgets                1.6.4     2023-12-06 [2] CRAN (R 4.5.0)
# httr                       1.4.7     2023-08-15 [2] CRAN (R 4.5.0)
# IRanges                  * 2.42.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# jpeg                       0.1-11    2025-03-21 [2] CRAN (R 4.5.0)
# jsonlite                   2.0.0     2025-03-27 [2] CRAN (R 4.5.0)
# KernSmooth                 2.23-26   2025-01-01 [3] CRAN (R 4.5.0)
# labeling                   0.4.3     2023-08-29 [2] CRAN (R 4.5.0)
# lattice                    0.22-7    2025-04-02 [3] CRAN (R 4.5.0)
# LearnBayes                 2.15.1    2018-03-18 [2] CRAN (R 4.5.0)
# lifecycle                  1.0.4     2023-11-07 [2] CRAN (R 4.5.0)
# limma                      3.64.3    2025-08-03 [1] Bioconductor 3.21 (R 4.5.0)
# locfit                     1.5-9.12  2025-03-05 [2] CRAN (R 4.5.0)
# magick                     2.8.6     2025-03-23 [2] CRAN (R 4.5.0)
# magrittr                   2.0.3     2022-03-30 [2] CRAN (R 4.5.0)
# MASS                       7.3-65    2025-02-28 [3] CRAN (R 4.5.0)
# Matrix                   * 1.7-3     2025-03-11 [3] CRAN (R 4.5.0)
# MatrixGenerics           * 1.20.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# matrixStats              * 1.5.0     2025-01-07 [2] CRAN (R 4.5.0)
# multcomp                   1.4-28    2025-01-29 [2] CRAN (R 4.5.0)
# mvtnorm                    1.3-3     2025-01-10 [2] CRAN (R 4.5.0)
# nlme                       3.1-168   2025-03-31 [3] CRAN (R 4.5.0)
# pillar                     1.10.2    2025-04-05 [2] CRAN (R 4.5.0)
# pkgconfig                  2.0.3     2019-09-22 [2] CRAN (R 4.5.0)
# png                        0.1-8     2022-11-29 [2] CRAN (R 4.5.0)
# proxy                      0.4-27    2022-06-09 [2] CRAN (R 4.5.0)
# purrr                      1.0.4     2025-02-05 [2] CRAN (R 4.5.0)
# R.methodsS3                1.8.2     2022-06-13 [2] CRAN (R 4.5.0)
# R.oo                       1.27.1    2025-05-02 [2] CRAN (R 4.5.0)
# R.utils                    2.13.0    2025-02-24 [2] CRAN (R 4.5.0)
# R6                         2.6.1     2025-02-15 [2] CRAN (R 4.5.0)
# ragg                       1.4.0     2025-04-10 [2] CRAN (R 4.5.0)
# RColorBrewer               1.1-3     2022-04-03 [2] CRAN (R 4.5.0)
# Rcpp                       1.1.1-1.1 2026-04-24 [1] CRAN (R 4.5.0)
# RCurl                      1.98-1.17 2025-03-22 [2] CRAN (R 4.5.0)
# rhdf5                    * 2.52.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# rhdf5filters               1.20.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# Rhdf5lib                   1.30.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# rjson                      0.2.23    2024-09-16 [2] CRAN (R 4.5.0)
# rlang                      1.1.6     2025-04-11 [2] CRAN (R 4.5.0)
# rprojroot                  2.0.4     2023-11-05 [2] CRAN (R 4.5.0)
# s2                         1.1.8     2025-05-12 [2] CRAN (R 4.5.0)
# S4Arrays                 * 1.8.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# S4Vectors                * 0.48.1    2026-04-05 [1] Bioconductor 3.22 (R 4.5.0)
# sandwich                   3.1-1     2024-09-15 [2] CRAN (R 4.5.0)
# scales                     1.4.0     2025-04-24 [2] CRAN (R 4.5.0)
# scuttle                    1.18.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# sessioninfo                1.2.3     2025-02-05 [2] CRAN (R 4.5.0)
# sf                         1.0-21    2025-05-15 [2] CRAN (R 4.5.0)
# sfheaders                  0.4.4     2024-01-17 [2] CRAN (R 4.5.0)
# SingleCellExperiment     * 1.30.1    2025-05-07 [2] Bioconductor 3.21 (R 4.5.0)
# sp                         2.2-0     2025-02-01 [2] CRAN (R 4.5.0)
# SparseArray              * 1.8.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# sparseMatrixStats          1.20.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# SpatialExperiment          1.18.1    2025-05-11 [2] Bioconductor 3.21 (R 4.5.0)
# SpatialFeatureExperiment   1.10.1    2025-05-11 [2] Bioconductor 3.21 (R 4.5.0)
# spatialreg                 1.3-6     2024-12-02 [2] CRAN (R 4.5.0)
# spData                     2.3.4     2025-01-08 [2] CRAN (R 4.5.0)
# spdep                      1.3-11    2025-04-24 [2] CRAN (R 4.5.0)
# statmod                    1.5.0     2023-01-06 [2] CRAN (R 4.5.0)
# SummarizedExperiment     * 1.38.1    2025-04-30 [2] Bioconductor 3.21 (R 4.5.0)
# survival                   3.8-3     2024-12-17 [3] CRAN (R 4.5.0)
# systemfonts                1.2.3     2025-04-30 [2] CRAN (R 4.5.0)
# terra                      1.8-50    2025-05-09 [2] CRAN (R 4.5.0)
# textshaping                1.0.1     2025-05-01 [2] CRAN (R 4.5.0)
# TH.data                    1.1-3     2025-01-17 [2] CRAN (R 4.5.0)
# tibble                     3.2.1     2023-03-20 [2] CRAN (R 4.5.0)
# tidyr                    * 1.3.1     2024-01-24 [2] CRAN (R 4.5.0)
# tidyselect                 1.2.1     2024-03-11 [2] CRAN (R 4.5.0)
# tiff                       0.1-12    2023-11-28 [2] CRAN (R 4.5.0)
# UCSC.utils                 1.4.0     2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# units                      0.8-7     2025-03-11 [2] CRAN (R 4.5.0)
# utf8                       1.2.5     2025-05-01 [2] CRAN (R 4.5.0)
# vctrs                      0.6.5     2023-12-01 [2] CRAN (R 4.5.0)
# withr                      3.0.2     2024-10-28 [2] CRAN (R 4.5.0)
# wk                         0.9.4     2024-10-11 [2] CRAN (R 4.5.0)
# XVector                    0.48.0    2025-04-15 [2] Bioconductor 3.21 (R 4.5.0)
# zeallot                    0.1.0     2018-01-28 [2] CRAN (R 4.5.0)
# zoo                        1.8-14    2025-04-10 [2] CRAN (R 4.5.0)
# 
# [1] /users/rphillip/R/4.5
# [2] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/library
# * ── Packages attached to the search path.
# 
# ─────────────────────────────────────────────────────────────────────────────────────
# 
