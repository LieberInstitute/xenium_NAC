#Goal for this script: Generate MEX files and other essential files for the Xenium Panel Designer
#Scripts modified from https://www.10xgenomics.com/analysis-guides/creating-single-cell-references-for-xenium-custom-panel-design-from-seurat-or-anndata 
#and
#https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/09_xenium_panel/format_for%20_10x_submission.R
#cd ~/NAc_Xenium_Panel
#module load r_nac

#load libraries
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(readxl)
library(dplyr)
library(scran)
library(here)

#Load the v3 panel that currently contains the genes, ranks, and other statistics. 
V3_Panel <- as.data.frame(read_excel(path = here("Tables","PanelV3_Ranks.xlsx")))

colnames(V3_Panel)
# [1] "Gene"                "Rank"                "Class"              
# [4] "nnsvg_prop_sig_adj"  "nnsvg_prop_top_100"  "nnsvg_avg_rank"     
# [7] "nnsvg_avg_rank_rank" "Cluster_DE"          "Cluster_Marked"     
# [10] "Base_Panel"          "Notes"    

#Just pull the Gene name and Rank 
V3_Panel <- V3_Panel[,c("Gene","Rank")]

#Need to add the ensembl gene IDs to this. 

#load the sce object. 
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

identical(rownames(colData(sce)),colnames(sce))
#[1] TRUE

#Now add ensemble gene id to the V3_Panel dataframe
dim(V3_Panel)
#[1] 100   2

#Add the ensembl gene id
V3_Panel <- merge(x = V3_Panel,
                  y = rowData(sce)[,c("gene_id","gene_name")],
                  by.x = "Gene",
                  by.y = "gene_name")

dim(V3_Panel)
#[1] 100   3

#Rearrange and order the dataframe
V3_Panel <- V3_Panel[order(V3_Panel$Rank),c("gene_id","Rank","Gene")]

head(V3_Panel)
# gene_id      Rank        Gene
# <character> <numeric> <character>
# 1 ENSG00000184845         1        DRD1
# 2 ENSG00000149295         2        DRD2
# 3 ENSG00000101327         3        PDYN
# 4 ENSG00000128271         4     ADORA2A
# 5 ENSG00000181195         5        PENK
# 6 ENSG00000112038         6       OPRM1

#Write out this file. 
write.csv(x = V3_Panel,
          file = here("Tables","NAc_Xenium_Panel_Final.csv"),
          row.names = FALSE,quote = FALSE)

#Generate MEX files needed for Xenium Panel Designer upload. 
#Check CheckIntegerCounts.log for proof that all values are integers in the counts(sce) matrix. 

#Before moving to writing out the MEX files, check that there are no duplicated barcodes. 
#This is because there needs to be 1 unique barcode per column of the count matrix. 
any(duplicated(colnames(sce)))
#[1] FALSE

#Another way to check this. 
length(colnames(counts(sce))) == length(unique(colnames(sce)))
#[1] TRUE

#Unfortunately, the files will be too large for Xenium Panel Designer.
#the recommendation via 10x is ~50k cells. I have attempted several previous subsets and the only way this will work 
#is to do a 30% subset. 
table(sce$CellType.Final)
# Astrocyte_A  Astrocyte_B   DRD1_MSN_A   DRD1_MSN_B   DRD1_MSN_C   DRD1_MSN_D 
# 8682         1250        22801         6544         2653         1370 
# DRD2_MSN_A   DRD2_MSN_B  Endothelial    Ependymal   Excitatory        Inh_A 
# 22704         1841          879          957         1583         1988 
# Inh_B        Inh_C        Inh_D        Inh_E        Inh_F    Microglia 
# 698         1456          226          554         1316         4438 
# Neuron_Ambig        Oligo          OPC 
# 446        17636         3763 

#First go ahead and remove the Neuron_Ambig cluster. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

# Extract cell type information
cell_metadata <- colData(sce)
cell_types <- cell_metadata$CellType.Final

# Split cells by type
split_cells <- split(seq_along(cell_types), cell_types)
#Creates a list where each part of the list are cells within that celltype/cluster. 

# Sample 50% from each group
set.seed(445) # For reproducibility
sampled_cells <- unlist(lapply(split_cells, function(cells) {
  sample(cells, size = length(cells) * 0.30)
}))

length(sampled_cells)
#[1] 30992

#Subset to half. 
sce_subset <- sce[, sampled_cells]

#Check that the strategy worked. 
table(sce_subset$CellType.Final)
# Astrocyte_A Astrocyte_B  DRD1_MSN_A  DRD1_MSN_B  DRD1_MSN_C  DRD1_MSN_D 
# 2604         375        6840        1963         795         411 
# DRD2_MSN_A  DRD2_MSN_B Endothelial   Ependymal  Excitatory       Inh_A 
# 6811         552         263         287         474         596 
# Inh_B       Inh_C       Inh_D       Inh_E       Inh_F   Microglia 
# 209         436          67         166         394        1331 
# Oligo         OPC 
# 5290        1128 

#Check that all is good with the subsetted object before writing out the files. 
#Are rownames of colData in same order as colnames of the counts matrix? 
identical(rownames(colData(sce_subset)),colnames(sce_subset))
#[1] TRUE

#Does colnames of the object call the colnames of the counts matrix? 
all(colnames(sce_subset) == colnames(counts(sce_subset)))
#[1] TRUE

#Are the ensembl gene ids in the rowData in the same order as the rownames of the object? 
all(rowData(sce_subset)$gene_id == rownames(sce_subset))
#[1] TRUE

#Are the ensembl gene ids in the rowData in the same order as the rownames of the counts matrix? 
all(rowData(sce_subset)$gene_id == rownames(counts(sce_subset)))
#[1] TRUE

#are there any duplicated barcodes? 
any(duplicated(colnames(sce_subset)))
#[1] FALSE

#another way to check duplicated barcodes?
length(colnames(sce_subset)) == length(unique(colnames(sce_subset)))
#[1] TRUE

#Use Droplet.Utils write10xCounts() to make the mex files. 
write10xCounts(path = here("Panel_Design_Files"),
               x = counts(sce_subset),
               gene.id = rowData(sce_subset)$gene_id,
               gene.symbol = rowData(sce_subset)$gene_name,
               barcodes = colnames(sce_subset),
               type = "sparse",
               version = "3")

list.files(here("Panel_Design_Files"))
#[1] "barcodes.tsv.gz" "features.tsv.gz" "matrix.mtx.gz"  

#Below is directly from: https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/09_xenium_panel/format_for%20_10x_submission.R
# ====== Adding spatial domain annotations ======
# bundleOutputs is a function provided by 10x. Note that their code had a typo where
# the "data" argument was not actually called and instead directly used "seurat_obj".
# 
# I have corrected this below.

# Define function
bundleOutputs <- function(out_dir, data, barcodes = colnames(data), cell_type = "cell_type", subset = 1:length(barcodes)) {
  
  if (require("data.table", quietly = TRUE)) {
    data.table::fwrite(
      data.table::data.table(
        barcode = barcodes,
        annotation = unlist(data[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv")
    )
  } else {
    write.table(
      data.frame(
        barcode = barcodes,
        annotation = unlist(data[[cell_type]])
      )[subset, ],
      file.path(out_dir, "annotations.csv"),
      sep = ",", row.names = FALSE
    )
  }
  
  bundle <- file.path(out_dir, paste0(basename(out_dir), ".zip"))
  
  utils::zip(
    bundle,
    list.files(out_dir, full.names = TRUE),
    zip = "zip"
  )
  
  if (file.info(bundle)$size / 1e6 > 500) {
    warning("The output file is more than 500 MB and will need to be subset further.")
  }
}

bundleOutputs(out_dir = here("Panel_Design_Files"), 
              data = sce_subset, 
              barcodes = colnames(sce_subset),
              cell_type = "CellType.Final")

list.files(here("Panel_Design_Files"))
# [1] "annotations.csv"        "barcodes.tsv.gz"        "features.tsv.gz"       
# [4] "matrix.mtx.gz"          "Panel_Design_Files.zip"

# read annotations.csv to check
test <- read.csv(here("Panel_Design_Files","annotations.csv"))
head(test)
#                 barcode  annotation
# 1 10_CCGTAGGAGCTGAGTG-1 Astrocyte_A
# 2 11_GGGACCTCAAAGCACG-1 Astrocyte_A
# 3 17_TCATTCAGTTCCTAAG-1 Astrocyte_A
# 4  9_AGTGCCGGTTGGTGTT-1 Astrocyte_A
# 5 11_ATGTCCCTCCGATGCG-1 Astrocyte_A
# 6 15_CCTAAGACACCATTCC-1 Astrocyte_A

# check to see if there are any duplicates in barcode
any(duplicated(test$barcode))
#[1] FALSE

length(unique(test$barcode)) == nrow(test)
#[1] TRUE

sessioninfo::session_info()
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2024-11-18
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────
# package              * version   date (UTC) lib source
# abind                  1.4-5     2016-07-21 [1] CRAN (R 4.3.2)
# beachmat               2.18.0    2023-10-24 [1] Bioconductor
# beeswarm               0.4.0     2021-06-01 [1] CRAN (R 4.3.2)
# Biobase              * 2.62.0    2023-10-24 [1] Bioconductor
# BiocGenerics         * 0.48.1    2023-11-01 [1] Bioconductor
# BiocNeighbors          1.20.2    2024-01-07 [1] Bioconductor 3.18 (R 4.3.2)
# BiocParallel           1.36.0    2023-10-24 [1] Bioconductor
# BiocSingular           1.18.0    2023-10-24 [1] Bioconductor
# bitops                 1.0-7     2021-04-24 [1] CRAN (R 4.3.2)
# bluster                1.11.4    2024-02-02 [1] Github (LTLA/bluster@17dd9c8)
# cellranger             1.1.0     2016-07-27 [1] CRAN (R 4.3.0)
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# data.table           * 1.14.10   2023-12-08 [1] CRAN (R 4.3.2)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# dplyr                * 1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
# DropletUtils         * 1.22.0    2023-10-24 [1] Bioconductor
# edgeR                  4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# generics               0.1.3     2022-07-05 [1] CRAN (R 4.3.0)
# GenomeInfoDb         * 1.38.1    2023-11-08 [1] Bioconductor
# GenomeInfoDbData       1.2.11    2023-12-12 [1] Bioconductor
# GenomicRanges        * 1.54.1    2023-10-29 [1] Bioconductor
# ggbeeswarm             0.7.2     2023-04-29 [1] CRAN (R 4.3.2)
# ggplot2              * 3.5.1     2024-04-23 [1] CRAN (R 4.3.2)
# ggrepel                0.9.4     2023-10-13 [1] CRAN (R 4.3.2)
# glue                   1.7.0     2024-01-09 [1] CRAN (R 4.3.2)
# gridExtra              2.3       2017-09-09 [1] CRAN (R 4.3.2)
# gtable                 0.3.4     2023-08-21 [1] CRAN (R 4.3.1)
# HDF5Array              1.30.0    2023-10-24 [1] Bioconductor
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# igraph                 2.0.3     2024-03-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# lattice                0.22-5    2023-10-24 [1] CRAN (R 4.3.1)
# lifecycle              1.0.4     2023-11-07 [1] CRAN (R 4.3.2)
# limma                  3.58.1    2023-10-31 [1] Bioconductor
# locfit                 1.5-9.8   2023-06-11 [1] CRAN (R 4.3.2)
# magrittr               2.0.3     2022-03-30 [1] CRAN (R 4.3.0)
# Matrix                 1.6-4     2023-11-30 [1] CRAN (R 4.3.2)
# MatrixGenerics       * 1.14.0    2023-10-24 [1] Bioconductor
# matrixStats          * 1.2.0     2023-12-11 [1] CRAN (R 4.3.2)
# metapod                1.10.0    2023-10-24 [1] Bioconductor
# munsell                0.5.0     2018-06-12 [1] CRAN (R 4.3.0)
# pillar                 1.9.0     2023-03-22 [1] CRAN (R 4.3.0)
# pkgconfig              2.0.3     2019-09-22 [1] CRAN (R 4.3.0)
# R.methodsS3            1.8.2     2022-06-13 [1] CRAN (R 4.3.2)
# R.oo                   1.26.0    2024-01-24 [1] CRAN (R 4.3.2)
# R.utils                2.12.3    2023-11-18 [1] CRAN (R 4.3.2)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# readxl               * 1.4.3     2023-07-06 [1] CRAN (R 4.3.0)
# rhdf5                  2.46.1    2023-11-29 [1] Bioconductor 3.18 (R 4.3.2)
# rhdf5filters           1.14.1    2023-11-06 [1] Bioconductor
# Rhdf5lib               1.24.1    2023-12-11 [1] Bioconductor 3.18 (R 4.3.2)
# rlang                  1.1.3     2024-01-10 [1] CRAN (R 4.3.2)
# rprojroot              2.0.4     2023-11-05 [1] CRAN (R 4.3.2)
# rsvd                   1.0.5     2021-04-16 [1] CRAN (R 4.3.2)
# S4Arrays               1.2.0     2023-10-24 [1] Bioconductor
# S4Vectors            * 0.40.2    2023-11-23 [1] Bioconductor 3.18 (R 4.3.2)
# ScaledMatrix           1.10.0    2023-10-24 [1] Bioconductor
# scales                 1.3.0     2023-11-28 [1] CRAN (R 4.3.2)
# scater               * 1.30.1    2023-11-16 [1] Bioconductor
# scran                * 1.30.0    2023-10-24 [1] Bioconductor
# scuttle              * 1.12.0    2023-10-24 [1] Bioconductor
# sessioninfo            1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# tibble                 3.2.1     2023-03-20 [1] CRAN (R 4.3.0)
# tidyselect             1.2.0     2022-10-10 [1] CRAN (R 4.3.0)
# utf8                   1.2.4     2023-10-22 [1] CRAN (R 4.3.1)
# vctrs                  0.6.5     2023-12-01 [1] CRAN (R 4.3.2)
# vipor                  0.4.5     2017-03-22 [1] CRAN (R 4.3.2)
# viridis                0.6.4     2023-07-22 [1] CRAN (R 4.3.2)
# viridisLite            0.4.2     2023-05-02 [1] CRAN (R 4.3.0)
# withr                  2.5.2     2023-10-30 [1] CRAN (R 4.3.1)
# XVector                0.42.0    2023-10-24 [1] Bioconductor
# zlibbioc               1.48.0    2023-10-24 [1] Bioconductor
# 
# [1] /jhpce/shared/libd/core/r_nac/1.0/nac_env/lib/R/library
# 
# ──────────────────────────────────────────────────────────────────────────────────────────────────
# 
# 
# 
# 

