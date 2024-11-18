#Goal: Compute a gene-gene correlation matrices for top 25 genes for each cluster. Each heatmap
#will only contain top 25 genes for a single cluster. 
# cd ~/NAc_Xenium_Panel/
# module load conda_R/4.4

#Load libraries
library(SingleCellExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(circlize)
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
snRNA_top25 <- subset(snRNA_markers,subset=(rank_marker %in% 1:25))

dim(snRNA_top25)
#[1] 500   9
#500 genes because 25*20 

#How many genes are unique
length(unique(snRNA_top25$gene_name))
#[1] 455
#45 shared

#Find which genes are duplicated 
dup_genes <- snRNA_top25[which(duplicated(snRNA_top25$gene_name)),"gene_name"]
dup_genes 
# [1] "HTR2C"      "HTR4"       "COCH"       "AC008574.1" "ANO3"      
# [6] "CNTN5"      "AQP4"       "KIRREL1"    "RXFP1"      "CUX2"      
# [11] "PLCXD3"     "ZNF804B"    "TRHDE"      "SEMA5B"     "GRIK3"     
# [16] "FOXP2"      "TSHZ1"      "TLL1"       "KCNC2"      "ADARB2"    
# [21] "KCNJ3"      "CALB2"      "LHX6"       "SFTA3"      "CHRM2"     
# [26] "SAMD3"      "LHX6"       "NXPH2"      "CHRM2"      "AC007368.1"
# [31] "NELL1"      "CLSTN2"     "KCNC2"      "RAB3B"      "CHRM2"     
# [36] "KCNJ3"      "ELAVL2"     "KCNJ6"      "SIDT1"      "EPHA5"     
# [41] "SLC27A6"    "AFF2"       "SPHKAP"     "NXPH1"      "AL445623.2"

#Which genes and which cell types
subset(snRNA_top25,subset=(gene_name %in% dup_genes))[,c("gene_name","cellType.target")]
#         gene_name cellType.target
# 34987  AC008574.1      DRD1_MSN_A
# 34992       HTR2C      DRD1_MSN_A
# 34994        COCH      DRD1_MSN_A
# 34996        ANO3      DRD1_MSN_A
# 34998        HTR4      DRD1_MSN_A
# 139918      HTR2C      DRD2_MSN_A
# 139919       HTR4      DRD2_MSN_A
# 139923       COCH      DRD2_MSN_A
# 139924 AC008574.1      DRD2_MSN_A
# 139926       ANO3      DRD2_MSN_A
# 209881       AQP4     Astrocyte_A
# 244840      RXFP1      DRD1_MSN_B
# 244841      TRHDE      DRD1_MSN_B
# 244844      FOXP2      DRD1_MSN_B
# 244845     SEMA5B      DRD1_MSN_B
# 244848    ZNF804B      DRD1_MSN_B
# 244849      CNTN5      DRD1_MSN_B
# 244851      TSHZ1      DRD1_MSN_B
# 244852       TLL1      DRD1_MSN_B
# 244855     ADARB2      DRD1_MSN_B
# 244862    KIRREL1      DRD1_MSN_B
# 314803      NXPH1           Inh_A
# 314804      SFTA3           Inh_A
# 314806       CUX2           Inh_A
# 314808      KCNC2           Inh_A
# 314810       LHX6           Inh_A
# 349774      GRIK3      DRD2_MSN_B
# 349788      CNTN5      DRD2_MSN_B
# 349793     PLCXD3      DRD2_MSN_B
# 384772       AQP4     Astrocyte_B
# 419728      EPHA5      DRD1_MSN_C
# 419731      KCNJ6      DRD1_MSN_C
# 419739    KIRREL1      DRD1_MSN_C
# 454704     CLSTN2      DRD1_MSN_D
# 454707      RXFP1      DRD1_MSN_D
# 454708      SAMD3      DRD1_MSN_D
# 454710       CUX2      DRD1_MSN_D
# 454712     PLCXD3      DRD1_MSN_D
# 454713    ZNF804B      DRD1_MSN_D
# 454716      TRHDE      DRD1_MSN_D
# 454718     SEMA5B      DRD1_MSN_D
# 454719      GRIK3      DRD1_MSN_D
# 454721      FOXP2      DRD1_MSN_D
# 454722      TSHZ1      DRD1_MSN_D
# 454724       TLL1      DRD1_MSN_D
# 454726      KCNJ3      DRD1_MSN_D
# 489684      KCNC2           Inh_B
# 489687 AC007368.1           Inh_B
# 489694     ADARB2           Inh_B
# 489699      CALB2           Inh_B
# 489702      KCNJ3           Inh_B
# 524659      NXPH2           Inh_C
# 524664      CALB2           Inh_C
# 524666      CHRM2           Inh_C
# 524668     ELAVL2           Inh_C
# 524670       LHX6           Inh_C
# 524675      SFTA3           Inh_C
# 524678 AL445623.2           Inh_C
# 524680      NELL1           Inh_C
# 559647     SPHKAP           Inh_D
# 559648      CHRM2           Inh_D
# 559649      SAMD3           Inh_D
# 594614       AFF2           Inh_E
# 594616       LHX6           Inh_E
# 594621      NXPH2           Inh_E
# 594622    SLC27A6           Inh_E
# 594627      RAB3B           Inh_E
# 594628      CHRM2           Inh_E
# 594629 AC007368.1           Inh_E
# 594631      NELL1           Inh_E
# 629603     CLSTN2      Excitatory
# 629610      SIDT1      Excitatory
# 664564      KCNC2           Inh_F
# 664566      RAB3B           Inh_F
# 664567      CHRM2           Inh_F
# 664569      KCNJ3           Inh_F
# 664570     ELAVL2           Inh_F
# 664571      KCNJ6           Inh_F
# 664574      SIDT1           Inh_F
# 664576      EPHA5           Inh_F
# 664577    SLC27A6           Inh_F
# 664578       AFF2           Inh_F
# 664579     SPHKAP           Inh_F
# 664580      NXPH1           Inh_F
# 664582 AL445623.2           Inh_F

#Make cellType.target a factor in the order that we want the genes listed on the heatmap. 
snRNA_top25$cellType.target <- factor(x = snRNA_top25$cellType.target,
                                      levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                 "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                                                 "Inh_C","Inh_D","Inh_E","Inh_F",
                                                 "Excitatory","Astrocyte_A","Astrocyte_B","Ependymal",
                                                 "Oligo","OPC","Microglia","Endothelial" ))

#Rearrange the dataframe in order of celltype
snRNA_top25 <- snRNA_top25[order(snRNA_top25$cellType.target),]

#subset for neurons only 
snRNA_top25 <- subset(snRNA_top25,subset=(cellType.target %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                                 "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                                                                 "Inh_C","Inh_D","Inh_E","Inh_F","Excitatory")))


#Pull the logcounts for the genes that we care about. 
expression_mat <- as.matrix(assay(sce,"logcounts")[snRNA_top25$gene_id,])

dim(expression_mat)
#[1]    325 103785

#Are the rownames of the expression matrix and the snRNA_top10 dataframe in the same order? 
identical(rownames(expression_mat),snRNA_top25$gene_id)
#[1] TRUE

#Okay, cool. Now let's go ahead and replace the rownames for the gene_name
rownames(expression_mat) <- snRNA_top25$gene_name

#Create gene-gene correlation matrix. 
gene_gene_cor <- cor(t(expression_mat),method = "spearman")

#Complex heatmap to visualize 
library(ComplexHeatmap)

#order of the rows
row_order <- c(rep("DRD1_MSN_A",25),
               rep("DRD1_MSN_B",25),
               rep("DRD1_MSN_C",25),
               rep("DRD1_MSN_D",25),
               rep("DRD2_MSN_A",25),
               rep("DRD2_MSN_B",25),
               rep("Inh_A",25),
               rep("Inh_B",25),
               rep("Inh_C",25),
               rep("Inh_D",25),
               rep("Inh_E",25),
               rep("Inh_F",25),
               rep("Excitatory",25))

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

#Keep only neuronal colors
cluster_cols <- cluster_cols[c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                               "DRD2_MSN_A","DRD2_MSN_B","Inh_A","Inh_B",
                               "Inh_C","Inh_D","Inh_E","Inh_F","Excitatory")]

cluster_cols
# DRD1_MSN_A DRD1_MSN_B DRD1_MSN_C DRD1_MSN_D DRD2_MSN_A DRD2_MSN_B      Inh_A 
# "#ECA31C"  "#D079AA"  "#2A4BF9"  "#FB3DD9"  "#58B6ED"  "#F80091"  "#35FB00" 
# Inh_B      Inh_C      Inh_D      Inh_E      Inh_F Excitatory 
# "#7A0096"  "#854222"  "#A7F281"  "#0DFBFA"    "black"  "#5C6300" 

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
                                "Excitatory" = "#5C6300"))

#Make a secondary color scheme for neurons or non-neurons
neuron_pops <- ifelse(row_order %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D"),
                      "D1_MSN",
                      ifelse(row_order %in% c("DRD2_MSN_A","DRD2_MSN_B"),
                             "D2_MSN",
                             ifelse(row_order %in% c("Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F"),
                                    "Inhibitory",
                                    "Excitatory")
                      )
)

neuron_pops <- factor(x = neuron_pops,levels = c("D1_MSN","D2_MSN","Inhibitory","Excitatory"))


colors_neurons <- list(class = c(D1_MSN = "dodgerblue",
                                 D2_MSN = "orange1",
                                 Inhibitory = "tomato",
                                 Excitatory = "#5C6300"))

# Create row annotation based on the identity and broad cell class
row_anno <- rowAnnotation(Identity = row_order, 
                          class = neuron_pops,
                          col = c(row_colors,colors_neurons))

# Draw heatmap with row annotation
hm <- Heatmap(gene_gene_cor, 
              name = "Spearman's correlation\ncoefficient (r)",
              cluster_columns=FALSE,
              cluster_rows=FALSE, #Leave the order the same.
              right_annotation = row_anno,
              row_names_gp = gpar(fontsize = 3),
              column_names_gp = gpar(fontsize = 3))

pdf(here("Plots","Neurons_Top25_gene_gene_correlation_Spearman.pdf"),height = 12, width = 12)
draw(hm,
     column_title = "Gene-Gene Correlation Matrix\nTop 25 marker genes for each neuronal cluster\nSpearman",
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
# circlize             * 0.4.16  2024-02-20 [2] CRAN (R 4.4.0)
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
