#Goal: Compute a gene-gene correlation matrix to identify sets of genes that provide similar information for 
#specific cell types. For example, GJA1 and AQP4 probably provide the same information for Astrocytes. 
#compute with SPEARMAN 
# cd ~/NAc_Xenium_Panel/
# module load conda_R/4.4

#Load libraries
library(SingleCellExperiment)
library(sessioninfo)
library(here)

#load the sce object to pull data from. 
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Inspect and double check that everything is in the correct order. 
sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(41): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

#Read in the 1vALL DEGs 
load(here("Tables","markers_1vAll_CellType_Final.rda"),
     verbose = TRUE)
# Loading objects:
#   markers_1vALL_df

#rename the table
snRNA_markers <- markers_1vALL_df
rm(markers_1vALL_df)

dim(snRNA_markers)
#[1] 699540      9

#subset for only enriched genes or a logFC > 0 
snRNA_markers <- subset(snRNA_markers,subset=(logFC > 0))

dim(snRNA_markers)
#[1] 264188      9

#subset for top 10
snRNA_top10 <- subset(snRNA_markers,subset=(rank_marker %in% 1:10))

dim(snRNA_top10)
#[1] 200   9
#200 genes because 10 genes x 20 clusters. 
#Removed Neuron_Ambig while performing DEG testing. 

length(unique(snRNA_top10$gene_name))
#[1] 197
#197 genes are unique, which are duplicated

#Find which genes are duplicated 
dup_genes <- snRNA_top10[which(duplicated(snRNA_top10$gene_name)),"gene_name"]

#Which genes and which cell types
subset(snRNA_top10,subset=(gene_name %in% dup_genes))
# gene_id    logFC log.p.value    log.FDR std.logFC
# 244840 ENSG00000171509 2.828411  -45064.830 -45054.368  4.891555
# 419731 ENSG00000157542 1.528172   -7374.659  -7366.142  2.495972
# 454707 ENSG00000171509 2.891297   -7481.184  -7472.513  3.475509
# 489684 ENSG00000166006 3.356360   -6805.685  -6797.014  4.586224
# 664564 ENSG00000166006 2.680402   -8232.843  -8222.380  3.713516
# 664571 ENSG00000157542 1.403954   -3041.303  -3032.920  2.198909
# cellType.target rank_marker         anno_logFC gene_name
# 244840      DRD1_MSN_B           1  std logFC = 4.892     RXFP1
# 419731      DRD1_MSN_C           7  std logFC = 2.496     KCNJ6
# 454707      DRD1_MSN_D           6  std logFC = 3.476     RXFP1
# 489684           Inh_B           6  std logFC = 4.586     KCNC2
# 664564           Inh_F           1  std logFC = 3.714     KCNC2
# 664571           Inh_F           8  std logFC = 2.199     KCNJ6

#Make cellType.target a factor in the order that we want the genes listed on the heatmap. 
snRNA_top10$cellType.target <- factor(x = snRNA_top10$cellType.target,
                                      levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                 "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                                                 "Inh_C","Inh_D","Inh_E","Inh_F",
                                                 "Excitatory","Astrocyte_A","Astrocyte_B","Ependymal",
                                                 "Oligo","OPC","Microglia","Endothelial" ))

#Rearrange the dataframe in order of celltype
snRNA_top10 <- snRNA_top10[order(snRNA_top10$cellType.target),]

#Pull the logcounts for the genes that we care about. 
expression_mat <- as.matrix(assay(sce,"logcounts")[snRNA_top10$gene_id,])

#Are the rownames of the expression matrix and the snRNA_top10 dataframe in the same order? 
identical(rownames(expression_mat),snRNA_top10$gene_id)
#[1] TRUE

#Okay, cool. Now let's go ahead and replace the rownames for the gene_name
rownames(expression_mat) <- snRNA_top10$gene_name

#Create gene-gene correlation matrix. 
gene_gene_cor <- cor(t(expression_mat),method = "spearman")

#Complex heatmap to visualize 
library(ComplexHeatmap)

#order of the rows
row_order <- c(rep("DRD1_MSN_A",10),
               rep("DRD1_MSN_B",10),
               rep("DRD1_MSN_C",10),
               rep("DRD1_MSN_D",10),
               rep("DRD2_MSN_A",10),
               rep("DRD2_MSN_B",10),
               rep("Inh_A",10),
               rep("Inh_B",10),
               rep("Inh_C",10),
               rep("Inh_D",10),
               rep("Inh_E",10),
               rep("Inh_F",10),
               rep("Excitatory",10),
               rep("Astrocyte_A",10),
               rep("Astrocyte_B",10),
               rep("Ependymal",10),
               rep("Oligo",10),
               rep("OPC",10),
               rep("Microglia",10),
               rep("Endothelial",10))

#load cluster colors 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)
# Loading objects:
#   cluster_cols

cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C Neuro_Ambig  DRD1_MSN_D       Inh_B       Inh_C       Inh_D 
# "#2A4BF9"      "grey"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281" 
# Inh_E  Excitatory       Inh_F 
# "#0DFBFA"   "#5C6300"     "black" 

#remove neuron_ambig
cluster_cols <- cluster_cols[-14]

cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C  DRD1_MSN_D       Inh_B       Inh_C       Inh_D       Inh_E 
# "#2A4BF9"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281"   "#0DFBFA" 
# Excitatory       Inh_F 
# "#5C6300"     "black" 

row_colors <- list(Identity = c("DRD1_MSN_A" = "#ECA31C",
                                "DRD1_MSN_B" =  "#D079AA",
                                "DRD1_MSN_C" = "#2A4BF9",  
                                "DRD1_MSN_D" = "#FB3DD9" ,
                                "DRD2_MSN_A" = "#58B6ED",
                                "DRD2_MSN_B" = "#F80091",
                                "Inh_A" = "#35FB00",
                                "Inh_B" = "#7A0096",
                                "Inh_C" = "#854222",
                                "Inh_D" = "#A7F281",     
                                "Inh_E" = "#0DFBFA",
                                "Inh_F" = "black",
                                "Excitatory" = "#5C6300",
                                "Astrocyte_A" =  "#D95F00",
                                "Astrocyte_B" = "#FF0016",
                                "Ependymal" = "#0077B9",  
                                "Oligo" = "#4F4753",
                                "OPC" = "#0D9F72",
                                "Microglia" =  "#F2E642",
                                "Endothelial" = "#D00DFF"))

#Make a secondary color scheme for neurons or non-neurons
neuron_pops <- ifelse(row_order %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D","DRD2_MSN_A",
                                       "DRD2_MSN_B","Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F","Excitatory"),
                      "Neuronal",
                      "Non-neuronal")

neuron_pops <- factor(x = neuron_pops,levels = c("Neuronal","Non-neuronal"))


colors_neurons <- list(class = c(Neuronal = "dodgerblue",
                                 `Non-neuronal` = "orange1"))


# Create row annotation based on the identity and broad cell class
row_anno <- rowAnnotation(Identity = row_order, 
                          class = neuron_pops,
                          col = c(row_colors,colors_neurons))

# Draw heatmap with row annotation
hm <- Heatmap(gene_gene_cor, 
              name = "Spearman's correlation\ncoefficient",
              cluster_columns=FALSE,
              cluster_rows=FALSE, #Leave the order the same.
              right_annotation = row_anno,
              row_names_gp = gpar(fontsize = 4.5),
              column_names_gp = gpar(fontsize = 4.5))

pdf(here("Plots","Top10_gene_gene_correlation_Spearman.pdf"),height = 12, width = 12)
draw(hm,
     column_title = "Gene-Gene Correlation Matrix\nTop 10 marker genes for each cluster\nSpearman",
     column_title_gp = gpar(fontsize = 16))
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.4.0 Patched (2024-05-22 r86590)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-10-25
# pandoc   3.1.13 @ /jhpce/shared/community/core/conda_R/4.4/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
# package              * version date (UTC) lib source
# abind                  1.4-5   2016-07-21 [2] CRAN (R 4.4.0)
# Biobase              * 2.64.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# BiocGenerics         * 0.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# Cairo                  1.6-2   2023-11-28 [2] CRAN (R 4.4.0)
# circlize               0.4.16  2024-02-20 [2] CRAN (R 4.4.0)
# cli                    3.6.3   2024-06-21 [2] CRAN (R 4.4.0)
# clue                   0.3-65  2023-09-23 [2] CRAN (R 4.4.0)
# cluster                2.1.6   2023-12-01 [3] CRAN (R 4.4.0)
# codetools              0.2-20  2024-03-31 [3] CRAN (R 4.4.0)
# colorspace             2.1-1   2024-07-26 [2] CRAN (R 4.4.0)
# ComplexHeatmap       * 2.20.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# crayon                 1.5.3   2024-06-20 [2] CRAN (R 4.4.0)
# DelayedArray           0.30.1  2024-05-07 [2] Bioconductor 3.19 (R 4.4.0)
# digest                 0.6.36  2024-06-23 [2] CRAN (R 4.4.0)
# doParallel             1.0.17  2022-02-07 [2] CRAN (R 4.4.0)
# foreach                1.5.2   2022-02-02 [2] CRAN (R 4.4.0)
# GenomeInfoDb         * 1.40.1  2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# GenomeInfoDbData       1.2.12  2024-05-23 [2] Bioconductor
# GenomicRanges        * 1.56.1  2024-06-12 [2] Bioconductor 3.19 (R 4.4.0)
# GetoptLong             1.0.5   2020-12-15 [2] CRAN (R 4.4.0)
# GlobalOptions          0.1.2   2020-06-10 [2] CRAN (R 4.4.0)
# here                 * 1.0.1   2020-12-13 [2] CRAN (R 4.4.0)
# httr                   1.4.7   2023-08-15 [2] CRAN (R 4.4.0)
# IRanges              * 2.38.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# iterators              1.0.14  2022-02-05 [2] CRAN (R 4.4.0)
# jsonlite               1.8.8   2023-12-04 [2] CRAN (R 4.4.0)
# lattice                0.22-6  2024-03-20 [3] CRAN (R 4.4.0)
# magick                 2.8.4   2024-07-14 [2] CRAN (R 4.4.0)
# magrittr               2.0.3   2022-03-30 [2] CRAN (R 4.4.0)
# Matrix                 1.7-0   2024-04-26 [3] CRAN (R 4.4.0)
# MatrixGenerics       * 1.16.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# matrixStats          * 1.3.0   2024-04-11 [2] CRAN (R 4.4.0)
# png                    0.1-8   2022-11-29 [2] CRAN (R 4.4.0)
# R6                     2.5.1   2021-08-19 [2] CRAN (R 4.4.0)
# RColorBrewer           1.1-3   2022-04-03 [2] CRAN (R 4.4.0)
# Rcpp                   1.0.13  2024-07-17 [2] CRAN (R 4.4.0)
# rjson                  0.2.21  2022-01-09 [2] CRAN (R 4.4.0)
# rprojroot              2.0.4   2023-11-05 [2] CRAN (R 4.4.0)
# S4Arrays               1.4.1   2024-05-20 [2] Bioconductor 3.19 (R 4.4.0)
# S4Vectors            * 0.42.1  2024-07-03 [2] Bioconductor 3.19 (R 4.4.0)
# sessioninfo          * 1.2.2   2021-12-06 [2] CRAN (R 4.4.0)
# shape                  1.4.6.1 2024-02-23 [2] CRAN (R 4.4.0)
# SingleCellExperiment * 1.26.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# SparseArray            1.4.8   2024-05-24 [2] Bioconductor 3.19 (R 4.4.0)
# SummarizedExperiment * 1.34.0  2024-05-01 [2] Bioconductor 3.19 (R 4.4.0)
# UCSC.utils             1.0.0   2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# XVector                0.44.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# zlibbioc               1.50.0  2024-04-30 [2] Bioconductor 3.19 (R 4.4.0)
# 
# [1] /users/rphillip/R/4.4
# [2] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/site-library
# [3] /jhpce/shared/community/core/conda_R/4.4/R/lib64/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────
