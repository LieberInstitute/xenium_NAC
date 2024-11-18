#Lines 5-163 from: https://github.com/LieberInstitute/spatial_NAc/blob/7da8b45efcf758834124061a3ab54ca9fc4f2f80/code/14_pseudobulk_spatial/04-plot_markers_visium.R#L18
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc
#module load r_nac

library(here)
library(PRECAST)
library(HDF5Array)
library(sessioninfo)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(spatialNAcUtils)
library(purrr)
library(ggpubr)
library(ggsci)
library(dittoSeq)
library(getopt)
library(pheatmap)
library(cowplot)

spot_plot2 <- function(spe, sample_id, image_id = "lowres",
                       title = sprintf("%s_%s", sample_id, var_name), var_name,
                       multi_gene_method = c("z_score", "pca", "sparsity"),
                       include_legend = TRUE, is_discrete, colors = NULL,
                       assayname = "logcounts", minCount = 0.5, spatial = FALSE) {
  #   This value was determined empirically, and results in good spot sizes.
  #   Note that it's sample-independent, and the final spot size to pass to
  #   'vis_gene' or 'vis_clus' uses this value along with the image
  #   dimensions, scale factors, and spatial coordinates for this particular
  #   sample
  IDEAL_POINT_SIZE <- 200
  
  ############################################################################
  #   Check validity of arguments
  ############################################################################
  
  # (Note that 'sample_id', 'var_name', 'assayname', 'minCount', and 'colors'
  # are not checked for validity here, since spatialLIBD functions handle
  # their validity)
  
  multi_gene_method <- rlang::arg_match(multi_gene_method)
  
  #   Check validity of spatial coordinates
  if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") == sort(colnames(spatialCoords(spe))))) {
    stop("Abnormal spatial coordinates: should have 'pxl_row_in_fullres' and 'pxl_col_in_fullres' columns.")
  }
  
  #   State assumptions about columns expected to be in the colData
  expected_cols <- c("array_row", "array_l", "sample_id", "exclude_overlapping")
  if (!all(expected_cols %in% colnames(colData(spe)))) {
    stop(
      sprintf(
        'Missing at least one of the following colData columns: "%s"',
        paste(expected_cols, collapse = '", "')
      )
    )
  }
  
  #   Subset to specific sample ID, and ensure overlapping spots are dropped
  subset_cols <- (spe$sample_id == sample_id) &
    (is.na(spe$exclude_overlapping) | !spe$exclude_overlapping)
  if (length(which(subset_cols)) == 0) {
    stop("No non-excluded spots belong to this sample. Perhaps check spe$exclude_overlapping for issues.")
  }
  spe_small <- spe[, subset_cols]
  
  ############################################################################
  #   Compute an appropriate spot size for this sample
  ############################################################################
  
  #   Determine some pixel values for the horizontal bounds of the spots
  MIN_COL <- min(spatialCoords(spe_small)[, "pxl_row_in_fullres"])
  MAX_COL <- max(spatialCoords(spe_small)[, "pxl_row_in_fullres"])
  
  #   The distance between spots (in pixels) is double the average distance
  #   between array columns
  INTER_SPOT_DIST_PX <- 2 * (MAX_COL - MIN_COL) /
    (max(spe_small$array_col) - min(spe_small$array_col))
  
  #   Find the appropriate spot size for this donor. This can vary because
  #   ggplot downscales a plot the fit desired output dimensions (in this
  #   case presumably a square region on a PDF), and stitched images can vary
  #   in aspect ratio. Also, lowres images always have a larger image
  #   dimension of 1200, no matter how many spots fit in either dimension.
  small_image_data <- imgData(spe_small)[
    imgData(spe_small)$image_id == image_id,
  ]
  spot_size <- IDEAL_POINT_SIZE * INTER_SPOT_DIST_PX *
    small_image_data$scaleFactor / max(dim(small_image_data$data[[1]]))
  
  
  ############################################################################
  #   Produce the plot
  ############################################################################
  
  #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
  #   'vis_gene'.
  if (is_discrete) {
    #   Supply a color scale if 'color' is not NULL. Otherwise, fall back
    #   upon 'vis_clus' defaults
    if (is.null(colors)) {
      p <- vis_clus(
        spe_small,
        sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
        return_plots = TRUE, spatial = spatial, point_size = spot_size
      )
    } else {
      p <- vis_clus(
        spe_small,
        sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
        return_plots TRUE, spatial = spatial, colors = colors,
        point_size = spot_size
      )
    }
  } else {
    if (is.null(colors)) {
      p <- vis_gene(
        spe_small,
        sampleid = sample_id, image_id = image_id, geneid = var_name,
        multi_gene_method = multi_gene_method, return_plots = TRUE,
        spatial = spatial, point_size = spot_size, assayname = assayname,
        alpha = NA, auto_crop = FALSE,
        minCount = minCount, na_color = NA
      )
    } else {
      p <- vis_gene(
        spe_small,
        sampleid = sample_id, image_id = image_id, geneid = var_name,
        multi_gene_method = multi_gene_method, return_plots = TRUE,
        spatial = spatial, point_size = spot_size, assayname = assayname,
        cont_colors = colors, alpha = NA, auto_crop = FALSE,
        minCount = minCount, na_color = NA
      )
    }
  }
  
  #   Remove the legend if requested
  if (!include_legend) {
    p <- p + theme(legend.position = "none")
  }
  
  #   Overwrite the title
  p <- p + labs(title = title)
  
  return(p)
}


select_samples <- c("Br6432", "Br6522", "Br3942","Br2720")

spe_dir <- here(
  "processed-data", "05_harmony_BayesSpace", "03-filter_normalize_spe", "spe_filtered_hdf5"
)

spe <- loadHDF5SummarizedExperiment(spe_dir)

# Subset to select samples
spe <- spe[ ,spe$donor %in% select_samples]
spe$exclude_overlapping[spe$sample_id_original == "V11D01-384_A1" & spe$overlap_slide == "V11D01-384_D1"] <- FALSE
spe$exclude_overlapping[spe$sample_id_original == "V11D01-384_D1" & spe$overlap_slide == "V11D01-384_A1"] <- TRUE

#cd
rownames(spe) <- rowData(spe)$gene_name

#load the sce object
library(scater)
library(scran)

#load the single cell objects
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

####Add some additional information to the panel document
#load in nnSVG results
nnSVG_precast_res <- read.csv(file = "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/07-run_nnSVG/nnSVG_precast_out/summary_across_samples.csv")

#add gene_name informaton to top100 df
nnSVG_precast_res <- dplyr::left_join(x = nnSVG_precast_res,
                                      y = as.data.frame(rowData(spe)[,c("gene_id","gene_name")]),
                                      by = "gene_id")


colnames(nnSVG_precast_res)[6] <- "Gene"

#Read in the version 3 panel
V3_Panel <- as.data.frame(readxl::read_excel(path = "~/NAc_Xenium_Panel/Tables/NAc_Xenium_Panel_v3.xlsx",sheet = "Panel_v3"))


#Add nnSVG stats information about each gene
#From: https://github.com/LieberInstitute/spatial_NAc/blob/main/code/05_harmony_BayesSpace/08-gather_nnSVG.R
# nnsvg_Prop_sig_adj = proportion of samples where the gene was significant
# nnsvg_prop_top_100 = proportion in the top 100 ranks
# nnsvg_avg_rank = average rank
# nnsvg_avg_rank_rank = rank of average ranks
V3_Panel <- dplyr::left_join(x = V3_Panel,
                             y = nnSVG_precast_res[,c(2:6)],
                             by = "Gene")

#Any gene not included in the nnSVG results will have NA for the nnSVG colummn 


#Load in the differential expression results from single cell analysis 
load("~/NAc_Xenium_Panel/Tables/markers_1vAll_CellType_Final.rda",
     verbose = TRUE)
# Loading objects:
#   marrs_1vALL_df


V3_Panel$Cluster_Marked <- NA
rownames(V3_Panel) <- V3_Panel$Gene
for(i in V3_Panel$Gene){
  print(i)
  x <- subset(markers_1vALL_df,subset=(gene_name == i & log.FDR <= log(0.05)) & std.logFC >= 0.5) 
  if(nrow(x)>0){
    V3_Panel[i,"Cluster_DE"] <- "Yes"
    V3_Panel[i,"Cluster_Marked"] <- ifelse(nrow(x) == 1,
                                           x$cellType.target,
                                           paste(x$cellType.target,collapse = ","))
  }else{
    V3_Panel[i,"Cluster_DE"] <- "No"
    V3_Panel[i,"Cluster_Marked"] <- NA
  }
}


#Also load in the base panel 
Base_Panel <- as.data.frame(readxl::read_excel(path = "~/NAc_Xenium_Panel/Tables/Xenium_Human_Brain_Panel_Gene_List.xlsx"))

#Are any of the genes in the custom panel included in the 
V3_Panel$Base_Panel <- ifelse(V3_Panel$Gene %in% Base_Panel$Gene,
                              "Yes",
                              "No")

#Ensure that none of the genes included are already in the base panel
table(V3_Panel$Base_Panel)
# No 
# 100 

#Now make spot plots and violins for all of the genes. 

######### VIOLIN PLOTS #########
#Load the cluster colors
load(here("processed-data","12_snRNA","070924_21colors_celltypeFinal.rda"),verbose = TRUE)
# Loading objects:
#   cluster_cols

for(i in V3_Panel$Gene){
  print(i)
  x <- plotExpression(object = sce,
                      features = i,
                      x = "CellType.Final",
                      colour_by = "CellType.Final",
                      swap_rownames = "gene_name") +
    scale_color_manual(values = cluster_cols) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3)
  ggplot2::ggsave(filename = paste0("~/NAc_Xenium_Panel/Plots/Panel_v3/Violin_Plots/",i,"_panelv3.png"),
                  plot = x,dpi = 100,width = 12,height = 8)
}


######### SPOT PLOTS #########
#-50 because CHAT not found in the object? 
for(i in V3_Panel$Gene[-50]){
  print(which(i == V3_Panel$Gene))
  print(i)
  #Make the plots for 3 donors that are anterior middle posterior 
  Br6432 <- spot_plot2(spe, sample_id = "Br6432", var_name = i, 
                       is_discrete = FALSE, spatial = TRUE, minCount = 0.5) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  Br6522 <- spot_plot2(spe, sample_id = "Br6522", var_name = i, 
                       is_discrete = FALSE, spatial = TRUE, minCount = 0.5) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
          panel.bord = element_rect(colour = "black", fill=NA, linewidth=0.5))
  Br3942 <- spot_plot2(spe, sample_id = "Br3942", var_name = i, 
                       is_discrete = FALSE, spatial = TRUE, minCount = 0.5) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  Br2720 <- spot_plot2(spe, sample_id = "Br2720", var_name = i, 
                       is_discrete = FALSE, spatial = TRUE, minCount = 0.5) + 
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm"), title =element_text(size=12, face='bold'), 
          panel.border = element_rect(colour = "black", fill=NA, linewidth=0.5))
  #Save the plot 
  x <- cowplot::plot_grid(plotlist = list(Br6432,Br6522,Br3942,Br2720),ncol = 2)
  ggplot2::ggsave(filename = paste0("~/NAc_Xenium_Panel/Plots/Panel_v3/Spot_Plots/",i,"_panelv3.png"),plot = x,dpi = 100,width = 12,height = 8)
}

######Save the panel
#reorder the datafarme
V3_Panel <- V3_Panel[,c("Gene","Class",
                        "nnsvg_prop_sig_adj","nnsvg_prop_top_100","nnsvg_avg_rank","nnsvg_avg_rank_rank",
                        "Cluster_DE","Cluster_Marked","Base_Panel","Notes")]

#Write out the panel
write.table(x = V3_Panel,file = "~/NAc_Xenium_Panel/Tables/NAc_Xenium_Panel_v3.txt",
          row.names = FALSE,quote = FALSE,sep = "\t")

sessioninfo::session_info()
# ─ Session info ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-11-08
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# package                * version    date (UTC) lib source
# abind                  * 1.4-5      2016-07-21 [1] CRAN (R 4.3.2)
# AnnotationDbi            1.64.1     2023-11-03 [1] Bioconductor
# AnnotationHub            3.10.0     2023-10-24 [1] Bioconductor
# attempt                  0.3.1      2020-05-03 [1] CRAN (R 4.3.2)
# backports                1.4.1      2021-12-13 [1] CRAN (R 4.3.0)
# beachmat                 2.18.0     2023-10-24 [1] Bioconductor
# beeswarm                 0.4.0      2021-06-01 [1] CRAN (R 4.3.2)
# benchmarkme              1.0.8      2022-06-12 [1] CRAN (R 4.3.2)
# benchmarkmeData          1.0.4      2020-04-23 [1] CRAN (R 4.3.2)
# Biobase                * 2.62.0     2023-10-24 [1] Bioconductor
# BiocFileCache            2.10.1     2023-10-26 [1] Bioconductor
# BiocGenerics           * 0.48.1     2023-11-01 [1] Bioconductor
# BiocIO                   1.12.0     2023-10-24 [1] Bioconductor
# BiocManager              1.30.22    2023-08-08 [1] CRAN (R 4.3.2)
# BiocNeighbors            1.20.2     2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel             1.36.0     2023-10-24 [1] Bioconductor
# BiocSingular             1.18.0     2023-10-24 [1] Bioconductor
# BiocVersion              3.18.1     2023-11-15 [1] Bioconductor
# Biostrings               2.70.1     2023-10-25 [1] Bioconductor
# bit                      4.0.5      2022-11-15 [1] CRAN (R 4.3.0)
# bit64                    4.0.5      2020-08-30 [1] CRAN (R 4.3.0)
# bitops                   1.0-7      2021-04-24 [1] CRAN (R 4.3.2)
# blob                     1.2.4      2023-03-17 [1] CRAN (R 4.3.0)
# bluster                  1.11.4     2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
# broom                    1.0.5      2023-06-09 [1] CRAN (R 4.3.0)
# bslib                    0.6.1      2023-11-28 [1] CRAN (R 4.3.2)
# cachem                   1.0.8      2023-05-01 [1] CRAN (R 4.3.0)
# car                      3.1-2      2023-03-30 [1] CRAN (R 4.3.2)
# carData                  3.0-5      2022-01-06 [1] CRAN (R 4.3.2)
# cellranger               1.1.0      2016-07-27 [1] CRAN (R 4.3.0)
# cli                      3.6.2      2023-12-11 [1] CRAN (R 4.3.2)
# cluster                  2.1.6      2023-12-01 [1] CRAN (R 4.3.2)
# codetools                0.2-19     2023-02-01 [1] CRAN (R 4.3.0)
# colorspace               2.1-0      2023-01-23 [1] CRAN (R 4.3.0)
# CompQuadForm             1.4.3      2017-04-12 [1] CRAN (R 4.3.2)
# config                   0.3.2      2023-08-30 [1] CRAN (R 4.3.2)
# cowplot                * 1.1.1      2020-12-30 [1] CRAN (R 4.3.2)
# crayon                   1.5.2      2022-09-29 [1] CRAN (R 4.3.0)
# curl                     5.2.0      2023-12-08 [1] CRAN (R 4.3.2)
# data.table               1.14.10    2023-12-08 [1] CRAN (R 4.3.2)
# DBI                    1.1.3      2022-06-18 [1] CRAN (R 4.3.0)
# dbplyr                   2.4.0      2023-10-26 [1] CRAN (R 4.3.1)
# DelayedArray           * 0.28.0     2023-10-24 [1] Bioconductor
# DelayedMatrixStats       1.24.0     2023-10-24 [1] Bioconductor
# deldir                   2.0-2      2023-11-23 [1] CRAN (R 4.3.2)
# digest                   0.6.33     2023-07-07 [1] CRAN (R 4.3.0)
# dittoSeq               * 1.14.3     2024-03-20 [1] Bioconductor 3.18 (R 4.3.2)
# doParallel               1.0.17     2022-02-07 [1] CRAN (R 4.3.2)
# dotCall64                1.1-1      2023-11-28 [1] CRAN (R 4.3.2)
# dplyr                  * 1.1.4      2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                    0.3.2      2023-11-29 [1] CRAN (R 4.3.2)
# DR.SC                    3.3        2023-08-09 [1] CRAN (R 4.3.2)
# DT                       0.31       2023-12-09 [1] CRAN (R 4.3.2)
# edgeR                    4.0.3      2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# ellipsis                 0.3.2      2021-04-29 [1] CRAN (R 4.3.0)
# ExperimentHub            2.10.0     2023-10-24 [1] Bioconductor
# fansi                    1.0.6      2023-12-08 [1] CRAN (R 4.3.2)
# farver                   2.1.1      2022-07-06 [1] CRAN (R 4.3.0)
# fastDummies              1.7.3      2023-07-06 [1] CRAN (R 4.3.2)
# fastmap                  1.1.1      2023-02-24 [1] CRAN (R 4.3.0)
# fields                   15.2       2023-08-17 [1] CRAN (R 4.3.2)
# filelock                 1.0.3      2023-12-11 [1] CRAN (R 4.3.2)
# fitdistrplus             1.1-11     20204-25 [1] CRAN (R 4.3.2)
# forcats                * 1.0.0      2023-01-29 [1] CRAN (R 4.3.0)
# foreach                  1.5.2      2022-02-02 [1] CRAN (R 4.3.0)
# future                   1.33.0     2023-07-01 [1] CRAN (R 4.3.0)
# future.apply             1.11.0     2023-05-21 [1] CRAN (R 4.3.0)
# generics                 0.1.3      2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb           * 1.38.1     2023-11-08 [1] Bioconductor
# GenomeInfoDbData         1.2.11     2023-12-12 [1] Bioconductor
# GenomicAlignments        1.38.0     2023-10-24 [1] Bioconductor
# GenomicRanges          * 1.54.1     2023-10-29 [1] Bioconductor
# getopt                 * 1.20.4     2023-10-01 [1] CRAN (R 4.3.2)
# ggbeeswarm               0.7.2      2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2                * 3.5.1      2024-04-23 [1] CRAN (R 4.3.2)
# ggpubr                 * 0.6.0      2023-02-10 [1] CRAN (R 4.3.2)
# ggrepel                  0.9.4      2023-10-13 [1] CRAN (R 4.3.2)
# ggridges                 0.5.4      2022-09-26 [1] CRAN (R 4.3.2)
# ggsci                  * 3.0.1      2024-03-02 [1] CRAN (R 4.3.2)
# ggsignif                 0.6.4      2022-10-13 [1] CRAN (R 4.3.2)
# ggthemes                 5.1.0      2024-02-10 [1] CRAN (R 4.3.2)
# GiRaF                    1.0.1      2020-10-14 [1] CRAN (R 4.3.2)
# globals                  0.16.2     2022-11-21 [1] CRAN (R 4.3.0)
# glue                     1.7.0      2024-01-09 [1] CRAN (R 4.3.2)
# goftest                  1.2-3      2021-10-07 [1] CRAN (R 4.3.2)
# golem                    0.4.1      2023-06-05 [1] CRAN (R 4.3.2)
# gridExtra                2.3        2017-09-09 [1] CRAN (R 4.3.2)
# gtable                   0.3.4      2023-08-21 [1] CRAN (R 4.3.1)
# gtools                 * 3.9.5      2023-11-20 [1] CRAN (R 4.3.2)
# HDF5Array              * 1.30.0     2023-10-24 [1] Bioconductor
# here                   * 1.0.1      2020-12-13 [1] CRAN (R 4.3.2)
# hms                      1.1.3      2023-03-21 [1] CRAN (R 4.3.0)
# htmltools                0.5.7      2023-11-03 [1] CRAN (R 4.3.2)
# htmlwidgets              1.6.4      2023-12-06 [1] CRAN (R 4.3.2)
# httpuv                   1.6.13     2023-12-06 [1] CRAN (R 4.3.2)
# httr                     1.4.7      2023-08-15 [1] CRAN (R 4.3.1)
# ica                      1.0-3      2022-07-08 [1] CRAN (R 4.3.2)
# igraph                   2.0.3      2024-03-13 [1] CRAN (R 4.3.2)
# interactiveDisplayBase   1.40.0   2023-10-24 [1] Bioconductor
# IRanges                * 2.36.0     2023-10-24 [1] Bioconductor
# irlba                    2.3.5.1    2022-10-03 [1] CRAN (R 4.3.2)
# iterators                1.0.14     2022-02-05 [1] CRAN (R 4.3.0)
# jquerylib                0.1.4      2021-04-26 [1] CRAN (R 4.3.0)
# jsonlite                 1.8.8      2023-12-04 [1] CRAN (R 4.3.2)
# KEGGREST                 1.42.0     2023-10-24 [1] Bioconductor
# KernSmooth               2.23-22    2023-07-10 [1] CRAN (R 4.3.0)
# labeling                 0.4.3      2023-08-29 [1] CRAN (R 4.3.1)
# later                    1.3.2      2023-12-06 [1] CRAN (R 4.3.2)
# lattice                  0.22-5     2023-10-24 [1] CRAN (R 4.3.1)
# lazyeval                 0.2.2      2019-03-15 [1] CRAN (R 4.3.0)
# leiden                   0.4.3.1    2023-11-17 [1] CRAN (R 4.3.2)
# lifecycle                1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
# limma                    3.58.1     2023-10-31 [1] Bioconductor
# listenv                  0.9.0      2022-12-16 [1] CRAN (R 4.3.0)
# lmtest                   0.9-40     2022-03-21 [1] CRAN (R 4.3.2)
# locfit                   1.5-9.8    2023-06-11 [1] CRAN (R 4.3.2)
# lubridate              * 1.9.3      2023-09-27 [1] CRAN (R 4.3.1)
# magick                   2.8.1      2023-10-22 [1] CRAN (R 4.3.2)
# magrittr                 2.0.3      2022-03-30 [1] CRAN (R 4.3.0)
# maps                   3.4.1.1    2023-11-03 [1] CRAN (R 4.3.2)
# MASS                     7.3-60     2023-05-04 [1] CRAN (R 4.3.0)
# Matrix                 * 1.6-4      2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics         * 1.14.0     2023-10-24 [1] Bioconductor
# matrixStats            * 1.2.0      2023-12-11 [1] CRAN (R 4.3.2)
# mclust                   6.0.1      2023-11-15 [1] CRAN (R 4.3.2)
# memoise                  2.0.1      2021-11-26 [1] CRAN (R 4.3.0)
# metapod                  1.10.0     2023-10-24 [1] Bioconductor
# mime                     0.12       2021-09-28 [1] CRAN (R 4.3.0)
# miniUI                   0.1.1.1    2018-05-18 [1] CRAN (R 4.3.2)
# munsell                  0.5.0      2018-06-12 [1] CRAN (R 4.3.0)
# nlme                     3.1-164    2023-11-27 [1] CRAN (R 4.3.2)
# paletteer                1.6.0      2024-01-21 [1] CRAN (R 4.3.2)
# parallelly               1.36.0     2023-05-26 [1] CRAN (R 4.3.0)
# patchwork                1.3.0.9000 2024-11-06 [1] Github (thomasp85/patchwork@2695a9f)
# pbapply                  1.7-2      2023-06-27 [1] CRAN (R 4.3.2)
# pheatmap               * 1.0.12     2019-01-04 [1] CRAN (R 4.3.2)
# pillar                   1.9.0      2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig                2.0.3      2019-09-22 [1] CRAN (R 4.3.0)
# plotly                   4.10.3     2023-10-21 [1] CRAN (R 4.3.2)
# plyr                     1.8.9      2023-10-02 [1] CRAN (R 4.3.1)
# png                      0.1-8      2022-11-29 [1] CRAN (R 4.3.2)
# polyclip                 1.10-6     2023-09-27 [1] CRAN (R 4.3.2)
# PRECAST                * 1.6.4      2024-01-25 [1] CRAN (R 4.3.2)
# progressr                0.14.0     2023-08-10 [1] CRAN (R 4.3.1)
# promises                 1.2.1      2023-08-10 [1] CRAN (R 4.3.1)
# purrr                  * 1.0.2      2023-08-10 [1] CRAN (R 4.3.1)
# R6                       2.5.1      2021-08-19 [1] CRAN (R 4.3.0)
# ragg                     1.2.7      2023-12-11 [1] CRAN (R 4.3.2)
# RANN                     2.6.1      2019-01-08 [1] CRAN (R 4.3.2)
# rappdirs                 0.3.3      2021-01-31 [1] CRAN (R 4.3.0)
# RColorBrewer             1.1-3      2022-04-03 [1] CRAN (R 4.3.0)
# Rcpp                     1.0.12     2024-01-09 [1] CRAN (R 4.3.2)
# RcppAnnoy                0.0.21     2023-07-02 [1] CRAN (R 4.3.2)
# RcppHNSW                 0.5.0      2023-09-19 [1] CRAN (R 4.3.2)
# RCurl                    1.98-1.13  2023-11-02 [1] CRAN (R 4.3.2)
# readr                  * 2.1.4      2023-02-10 [1] CRAN (R 4.3.0)
# readxl                   1.4.3      2023-07-06 [1] CRAN (R 4.3.0)
# rematch2                 2.1.2      2020-05-01 [1] CRAN (R 4.3.0)
# reshape2                 1.4.4      2020-04-09 [1] CRAN (R 4.3.0)
# restfulr                 0.0.15     2022-06-16 [1] CRAN (R 4.3.2)
# reticulate               1.34.0     2023-10-12 [1] CRAN (R 4.3.2)
# rhdf5                  * 2.46.1     2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters             1.14.1     2023-11-06 [1] Bioconductor
# Rhdf5lib                 1.24.1     2023-12-11 [1] Bioconductor 3.18 (R 4.3.2)
# rjson                    0.2.21     2022-01-09 [1] CRAN (R 4.3.2)
# rlang                    1.1.3      2024-01-10 [1] CRAN (R 4.3.2)
# ROCR                     1.0-11     2020-05-02 [1] CRAN (R 4.3.2)
# rprojroot                2.0.4      2023-11-05 [1] CRAN (R 4.3.2)
# Rsamtools                2.18.0     2023-10-24 [1] Bioconductor
# RSpectra                 0.16-1     2022-04-24 [1] CRAN (R 4.3.2)
# RSQLite                  2.3.4      2023-12-08 [1] CRAN (R 4.3.2)
# rstatix                  0.7.2      2023-02-01 [1] CRAN (R 4.3.2)
# rsvd                     1.0.5      2021-04-16 [1] CRAN (R 4.3.2)
# rtracklayer              1.62.0     2023-10-24 [1] Bioconductor
# Rtsne                    0.17       2023-12-07 [1] CRAN (R 4.3.2)
# S4Arrays               * 1.2.0      2023-10-24 [1] Bioconductor
# S4Vectors              * 0.40.2     2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# sass                     0.4.8      2023-12-06 [1] CRAN (R 4.3.2)
# ScaledMatrix             1.10.0     2023-10-24 [1] Bioconductor
# scales                   1.3.0      2023-11-28 [1] CRAN (R 4.3.2)
# scater                 * 1.30.1     2023-11-16 [1] Bioconductor
# scattermore              1.2        2023-06-12 [1] CRAN (R 4.3.2)
# scran                  * 1.30.0     2023-10-24 [1] Bioconductor
# sctransform              0.4.1      2023-10-19 [1] CRAN (R 4.3.2)
# scuttle                * 1.12.0     2023-10-24 [1] Bioconductor
# sessioninfo            * 1.2.2      2021-12-06 [1] CRAN (R 4.3.2)
# Seurat                   5.0.1      2023-11-17 [1] CRAN (R 4.3.2)
# SeuratObject             5.0.1      2023-11-17 [1] CRAN (R 4.3.2)
# shiny                    1.8.0      2023-11-17 [1] CRAN (R 4.3.2)
# shinyWidgets             0.8.0      2023-08-30 [1] CRAN (R 4.3.2)
# SingleCellExperiment   * 1.24.0     2023-10-24 [1] Bioconductor
# sp                       2.1-2      2023-11-26 [1] CRAN (R 4.3.2)
# spam                     2.10-0     2023-10-23 [1] CRAN (R 4.3.2)
# SparseArray            * 1.2.2      2023-11-07 [1] Bioconductor
# sparseMatrixStats        1.14.0     2023-10-24 [1] Bioconductor
# SpatialExperiment      * 1.12.0     2023-10-24 [1] Bioconductor
# spatialLIBD            * 1.17.8     2024-09-11 [1] Github (LieberInstitute/spatialLIBD@c82c789)
# spatialNAcUtils        * 0.99.0     2023-12-12 [1] Github (LieberInstitute/spatialNAcUtils@24c45ac)
# spatstat.data            3.0-3      2023-10-24 [1] CRAN (R 4.3.2)
# spatstat.explore         3.2-5      2023-10-22 [1] CRAN (R 4.3.2)
# spatstat.geom            3.2-7      2023-10-20 [1] CRAN (R 4.3.2)
# spatstat.random          3.2-2      2023-11-29 [1] CRAN (R 4.3.2)
# spatstat.sparse          3.0-3      2023-10-24 [1] CRAN (R 4.3.2)
# spatstat.utils           3.1-0      2024-08-17 [1] CRAN (R 4.3.2)
# statmod                  1.5.0      2023-01-06 [1] CRAN (R 4.3.2)
# stringi                  1.8.3      2023-12-11 [1] CRAN (R 4.3.2)
# stringr                * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
# SummarizedExperiment   * 1.32.0     2023-10-24 [1] Bioconductor
# survival                 3.5-7      2023-08-14 [1] CRAN (R 4.3.1)
# systemfonts              1.0.5      2023-10-09 [1] CRAN (R 4.3.1)
# tensor                   1.5        2012-05-05 [1] CRAN (R 4.3.2)
# textshaping              0.3.7      2023-10-09 [1] CRAN (R 4.3.1)
# tibble                 * 3.2.1      2023-03-20 [1] CRAN (R 4.3.0)
# tidyr                  * 1.3.0      2023-01-24 [1] CRAN (R 4.3.0)
# tidyselect               1.2.0      2022-10-10 [1] CRAN (R 4.3.0)
# tidyverse              * 2.0.0      2023-02-22 [1] CRA(R 4.3.0)
# timechange               0.2.0      2023-01-11 [1] CRAN (R 4.3.0)
# tzdb                     0.4.0      2023-05-12 [1] CRAN (R 4.3.0)
# utf8                     1.2.4      2023-10-22 [1] CRAN (R 4.3.1)
# uwot                     0.1.16     2023-06-29 [1] CRAN (R 4.3.2)
# vctrs                    0.6.5      2023-12-01 [1] CRAN (R 4.3.2)
# vipor                    0.4.5      2017-03-22 [1] CRAN (R 4.3.2)
# viridis                  0.6.4      2023-07-22 [1] CRAN (R 4.3.2)
# viridisLite              0.4.2      2023-05-02 [1] CRAN (R 4.3.0)
# withr                    2.5.2      2023-10-30 [1] CRAN (R 4.3.1)
# XML                      3.99-0.16  2023-11-29 [1] CRAN (R 4.3.2)
# xtable                   1.8-4      2019-04-21 [1] CRAN (R 4.3.0)
# XVector                  0.42.0     2023-10-24 [1] Bioconductor
# yaml                     2.3.8      2023-12-11 [1] CRAN (R 4.3.2)
# zlibbioc                 1.48.0     2023-10-24 [1] Bioconductor
# zoo                      1.8-12     2023-04-13 [1] CRAN (R 4.3.0)
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
