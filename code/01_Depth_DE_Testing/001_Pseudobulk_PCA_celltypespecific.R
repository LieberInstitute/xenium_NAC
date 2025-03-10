# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
#module load r_nac
#cd cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/




library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

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

#Remove the Neuron_Ambig group
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

dim(sce)
#[1]  36601 103339

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)
# FALSE  TRUE 
# 34977  1624 

#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

dim(sce)
#[1]  34977 103339

#For this each sample needs an anterior/middle/posterior designation. 
#Combine the Anterior and Posterior samples into a single group called "Anterior_Posterior"
#Make a dataframe. 
Ant_Mid_Post <- data.frame(Brain_ID = unique(sce$Brain_ID))

Ant_Mid_Post <- cbind(Ant_Mid_Post,c("Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior",
                                     "Anterior_Posterior","Anterior_Posterior",
                                     "Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior"))

colnames(Ant_Mid_Post)[2] <- "Depth"

Ant_Mid_Post
# Brain_ID              Depth
# 1    Br8325 Anterior_Posterior
# 2    Br8492             Middle
# 3    Br2720             Middle
# 4    Br6423 Anterior_Posterior
# 5    Br2743 Anterior_Posterior
# 6    Br3942 Anterior_Posterior
# 7    Br6432 Anterior_Posterior
# 8    Br6471             Middle
# 9    Br6522             Middle
# 10   Br8667 Anterior_Posterior

#Add depth to the sce object
sce$Depth <- Ant_Mid_Post[match(sce$Brain_ID,Ant_Mid_Post$Brain_ID),"Depth"]

as.data.frame(unique(colData(sce)[,c("Brain_ID","Depth")]))
# Brain_ID              Depth
# 1_AAACCCAAGACCAACG-1    Br8325 Anterior_Posterior
# 3_AAACCCAAGGTGAGCT-1    Br8492             Middle
# 5_AAACCCAGTAATTAGG-1    Br2720             Middle
# 7_AAACCCACACCCTTAC-1    Br6423 Anterior_Posterior
# 9_AAACCCACATTGTAGC-1    Br2743 Anterior_Posterior
# 11_AAACCCAAGACTCCGC-1   Br3942 Anterior_Posterior
# 13_AAACCCAAGGACTGGT-1   Br6432 Anterior_Posterior
# 15_AAACCCAAGGGAGATA-1   Br6471             Middle
# 17_AAACCCAAGATTGCGG-1   Br6522             Middle
# 19_AAACCCACAAGGCCTC-1   Br8667 Anterior_Posterior



#Pseudobulk across CellType and Brain_ID
sce_pb <- aggregateAcrossCells(sce,ids = colData(sce)[,c("CellType.Final","Brain_ID")])

sce_pb
# class: SingleCellExperiment 
# dim: 34977 199 
# metadata(1): Samples
# assays(1): counts
# rownames(34977): ENSG00000243485 ENSG00000186092 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames: NULL
# colData names(45): Sample Barcode ... Brain_ID ncells
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

table(sce_pb$CellType.Final)
# Astrocyte_A Astrocyte_B  DRD1_MSN_A  DRD1_MSN_B  DRD1_MSN_C  DRD1_MSN_D 
# 10          10          10          10          10          10 
# DRD2_MSN_A  DRD2_MSN_B Endothelial   Ependymal  Excitatory       Inh_A 
# 10          10          10           9          10          10 
# Inh_B       Inh_C       Inh_D       Inh_E       Inh_F   Microglia 
# 10          10          10          10          10          10 
# Oligo         OPC 
# 10          10 

table(sce_pb$Brain_ID)
# Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667 
# 20     20     20     20     20     20     20     20     19     20 

de.results <- pseudoBulkDGE(sce_pb,
                            label = sce_pb$CellType.Final,
                            design = ~Depth,
                            coef = "DepthMiddle",
                            condition = sce_pb$Depth)

saveRDS(de.results,file = here("processed-data","DE_Testing","Depth_DE_res.Rds"))

is.de <- decideTestsPerLabel(de.results, threshold=0.1)


x <- as.data.frame(summarizeTestsPerLabel(is.de))
#Remove the NA column
x <- x[,-4]
x$celltype <- rownames(x)
x_melt <- reshape2::melt(x)
#Using celltype as id variables

p <- ggplot(data = subset(x_melt,subset=(variable == -1 | variable == 1)),aes(x = celltype,y=value,fill = variable)) +
  geom_bar(stat="identity") +
  labs(x = "Cell Type",
       y = "# of Genes",
       fill = "logFC Direction") +
  ggtitle("# of DEGs (FDR<0.1)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(hjust = 0.5)) 
ggsave(p,filename = here("plots","01_Depth_DE_Testing","No_DEGs_0.1FDR_Bargraph.pdf"))
 
#Pull dataframe of enriched DEGs

#D1_B
D1_B_DEGs <- subset(de.results[["DRD1_MSN_B"]],subset=(FDR<=0.1 & logFC > 0))
D1_B_DEGs$gene_id <- rownames(D1_B_DEGs)
D1_B_DEGs <- merge(D1_B_DEGs,rowData(sce)[,c("gene_id","gene_name")],by = "gene_id")  
# DataFrame with 2 rows and 7 columns
# gene_id     logFC    logCPM         F      PValue       FDR
# <character> <numeric> <numeric> <numeric>   <numeric> <numeric>
#   1 ENSG00000064655   1.39426   5.95117   53.6043 1.05594e-05 0.0655246
# 2 ENSG00000184984   2.80088   4.84438   60.6905 5.72597e-06 0.0532973
# gene_name
# <character>
#   1        EYA2
# 2       CHRM5

#D1_D
D1_D_DEGs <- subset(de.results[["DRD1_MSN_D"]],subset=(FDR<=0.1 & logFC > 0 ))
D1_D_DEGs$gene_id <- rownames(D1_D_DEGs)
D1_D_DEGs <- merge(D1_D_DEGs,rowData(sce)[,c("gene_id","gene_name")],by = "gene_id")  
#DataFrame with 1 row and 7 columns
# gene_id     logFC    logCPM         F      PValue       FDR
# <character> <numeric> <numeric> <numeric>   <numeric> <numeric>
#   1 ENSG00000091831   1.66362   4.09871   29.5183 4.24559e-05 0.0704775
# gene_name
# <character>
#   1        ESR1


#D1_D
Inh_D_DEGs <- subset(de.results[["Inh_D"]],subset=(FDR<=0.1 & logFC > 0 ))
Inh_D_DEGs$gene_id <- rownames(Inh_D_DEGs)
Inh_D_DEGs <- merge(Inh_D_DEGs,rowData(sce)[,c("gene_id","gene_name")],by = "gene_id")  
Inh_D_DEGs$gene_name
# [1] "NEDD4L"     "ERICH1"     "OPRM1"      "ADRA1A"     "PDZRN3"    
# [6] "RASGEF1B"   "SLIT2"      "ADAM12"     "ADCY8"      "GRIK2"     
# [11] "PLEKHA7"    "SLCO3A1"    "COL25A1"    "RYR1"       "LPAR1"     
# [16] "AC023095.1"


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
# [1] "Reproducibility information:"
# [1] "2025-03-10 09:34:08 EDT"
# user   system  elapsed 
# 171.292    9.605 1387.027 
# ─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
# setting  value
# version  R version 4.3.2 (2023-10-31)
# os       Rocky Linux 9.4 (Blue Onyx)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# ctype    en_US.UTF-8
# tz       US/Eastern
# date     2025-03-10
# pandoc   3.1.3 @ /jhpce/shared/libd/core/r_nac/1.0/nac_env/bin/pandoc
# 
# ─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
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
# cli                    3.6.2     2023-12-11 [1] CRAN (R 4.3.2)
# cluster                2.1.6     2023-12-01 [1] CRAN (R 4.3.2)
# codetools              0.2-19    2023-02-01 [1] CRAN (R 4.3.0)
# colorspace             2.1-0     2023-01-23 [1] CRAN (R 4.3.0)
# crayon                 1.5.2     2022-09-29 [1] CRAN (R 4.3.0)
# DelayedArray           0.28.0    2023-10-24 [1] Bioconductor
# DelayedMatrixStats     1.24.0    2023-10-24 [1] Bioconductor
# dplyr                  1.1.4     2023-11-17 [1] CRAN (R 4.3.2)
# dqrng                  0.3.2     2023-11-29 [1] CRAN (R 4.3.2)
# edgeR                  4.0.3     2023-12-10 [1] Bioconductor 3.18 (R 4.3.2)
# fansi                  1.0.6     2023-12-08 [1] CRAN (R 4.3.2)
# farver                 2.1.1     2022-07-06 [1] CRAN (R 4.3.0)
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
# here                 * 1.0.1     2020-12-13 [1] CRAN (R 4.3.2)
# igraph                 2.0.3     2024-03-13 [1] CRAN (R 4.3.2)
# IRanges              * 2.36.0    2023-10-24 [1] Bioconductor
# irlba                  2.3.5.1   2022-10-03 [1] CRAN (R 4.3.2)
# labeling               0.4.3     2023-08-29 [1] CRAN (R 4.3.1)
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
# plyr                   1.8.9     2023-10-02 [1] CRAN (R 4.3.1)
# R6                     2.5.1     2021-08-19 [1] CRAN (R 4.3.0)
# ragg                   1.2.7     2023-12-11 [1] CRAN (R 4.3.2)
# Rcpp                   1.0.12    2024-01-09 [1] CRAN (R 4.3.2)
# RCurl                  1.98-1.13 2023-11-02 [1] CRAN (R 4.3.2)
# reshape2               1.4.4     2020-04-09 [1] CRAN (R 4.3.0)
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
# sessioninfo          * 1.2.2     2021-12-06 [1] CRAN (R 4.3.2)
# SingleCellExperiment * 1.24.0    2023-10-24 [1] Bioconductor
# SparseArray            1.2.2     2023-11-07 [1] Bioconductor
# sparseMatrixStats      1.14.0    2023-10-24 [1] Bioconductor
# statmod                1.5.0     2023-01-06 [1] CRAN (R 4.3.2)
# stringi                1.8.3     2023-12-11 [1] CRAN (R 4.3.2)
# stringr                1.5.1     2023-11-14 [1] CRAN (R 4.3.2)
# SummarizedExperiment * 1.32.0    2023-10-24 [1] Bioconductor
# systemfonts            1.0.5     2023-10-09 [1] CRAN (R 4.3.1)
# textshaping            0.3.7     2023-10-09 [1] CRAN (R 4.3.1)
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
# ──────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────────
