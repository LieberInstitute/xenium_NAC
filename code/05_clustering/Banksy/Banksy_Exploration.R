##Goal: Investigate Banksy louvain clustering. 
#lambda 0.8 + k_geom = 15
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe
# class: SpatialExperiment 
# dim: 541 4884175 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
# DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4884175): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Subset for donor Br6660
spe <- spe[,spe$Donor == "Br6660"]

spe
# class: SpatialExperiment 
# dim: 541 2444243 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
# DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2444243): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Read in the louvain clustering that was generatd with a resolution value of 1.0 
banksy_louvain_clusters <- read.csv(file = here("processed-data","05_Clustering","Banksy_louvain_clusters.csv"))

head(banksy_louvain_clusters)
# X V1                V2
# 1 1  1 Br6660_NAc1_580_1
# 2 2  1 Br6660_NAc1_580_2
# 3 3  2 Br6660_NAc1_580_3
# 4 4  1 Br6660_NAc1_580_4
# 5 5  3 Br6660_NAc1_580_5
# 6 6  1 Br6660_NAc1_580_6

banksy_louvain_clusters <- banksy_louvain_clusters[,-1]

colnames(banksy_louvain_clusters) <- c("Cluster","Key")

#How many clusters? How many cells in each cluster
table(banksy_louvain_clusters$Cluster)
# 1      2      3      4      5      6      7      8      9     10     11 
# 66237 155113  59329  96613  56579 157550  52003  83237 331771 316878  93658 
# 12     13     14     15     16     17     18     19     20     21     22 
# 53656  89748  43132 181037 113538  85597  34117  69190 123142  67299  69981 
# 23     24     25 
# 13130  31564    144 

dim(banksy_louvain_clusters)
#[1] 2444243       2

identical(banksy_louvain_clusters$Key,rownames(colData(spe)))
#[1] 2444243       2

identical(colnames(spe),banksy_louvain_clusters$Key)
#[1] TRUE

#Add the clusters to the object
spe$Banksy_louvain_res1.0 <- as.factor(banksy_louvain_clusters$Cluster)


#Create a color palette that makes sense. 
#Plot Banksy clusters on the tissue sections to help with annotation. 
cluster_cols <- Polychrome::createPalette(length(unique(spe$Banksy_louvain_res1.0)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(cluster_cols) <- unique(spe$Banksy_louvain_res1.0)

#Plot with escheR 
for (i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[, spe$Sample == i]
  
  p <-make_escheR(sub_spe) |>
    add_fill("Banksy_louvain_res1.0")+
    ggtitle(i) +
    scale_fill_manual(values = cluster_cols) +
    theme(element_text(hjust = 0.5))
  
  png(filename = here("plots","05_clustering","Banksy","Pre_Annotation","Louvain_res1.0",paste0(i,"_PreAnnotation_louvainres1.0.png")),
      units = "in",height = 20,width = 20,res = 300)
  print(p)
  dev.off()
  
}



table(spe$Sample,spe$Banksy_louvain_res1.0)
# 1     2     3     4     5     6     7     8     9    10
# Br6660_NAc1_580    6737   543  1305 14384  4784 21343 10240  7681 34172 26473
# Br6660_NAc2_1090   6760  3625  1547  8623  3707 13354  4263  5108 21751 17484
# Br6660_NAc3_1580  12709 21731  3547 10398  5110 17500  6627 11227 32598 28775
# Br6660_NAc4_2080  10344 12894  3595  9172  4789 14279  5080  7230 31571 28110
# Br6660_NAc5_2580  11450 14707  3385  8629  4866 16767  5448  6648 40774 38979
# Br6660_NAc6_3080   5526 10240  3553  8360  4441 12813  4176  6549 37494 36199
# Br6660_NAc7_3580   4488 18029  3681  7305  4279 10868  3350  6021 33041 35517
# Br6660_NAc8_4580   2410 21050 10551  6487  5735 11925  2175  8698 26877 28948
# Br6660_NAc9_5080   1675 19605 10975  8261  6584 13078  3850  8042 27054 26395
# Br6660_Nac10_4080  2910 11398  8657  7370  5873 11969  3360  8587 31216 36765
# Br6660_Nac11_5580  1228 21291  8533  7624  6411 13654  3434  7446 15223 13233
# 
# 11    12    13    14    15    16    17    18    19    20
# Br6660_NAc1_580    8197  4204  5760  7004  3302  7355  5043  2348  2692  5815
# Br6660_NAc2_1090   5297  2729  3863  2574  4026  4952  4659  1625  1590  3389
# Br6660_NAc3_1580   7356  4277  5562  4543  9788  9145  8523  2713  3194  6969
# Br6660_NAc4_2080   8119  4846  7554  2303  9442  9178  3603  3148  3139  8531
# Br6660_NAc5_2580  10213  6387  9430  4014  8014 11859  8295  3468  4597 12737
# Br6660_NAc6_3080   9202  6050  9067  1994 11744 11641 15155  2604  7371 14216
# Br6660_NAc7_3580   8936  5618  9503  1897 12601 11465 13501  2503  8477 12663
# Br6660_NAc8_4580   9290  4992 11199  4217 30674 12036  6225  3484  8934 14885
# Br6660_NAc9_5080   9727  5263  9357  4222 34302 12672  9191  4406  9759 15986
# Br6660_Nac10_4080 10147  5600  9811  7919 27679 13731  7980  3078 10893 16507
# Br6660_Nac11_5580  7174  3690  8642  2445 29465  9504  3422  4740  8544 11444
# 
#                       21    22    23    24    25
# Br6660_NAc1_580    2687  1347    93     0     0
# Br6660_NAc2_1090   1979  1504   134     7     0
# Br6660_NAc3_1580   5056  3905   721    22     0
# Br6660_NAc4_2080   5505 10788   125   658     0
# Br6660_NAc5_2580   8858  9364  1636  1553     0
# Br6660_NAc6_3080   8542  7608  1114  3671     0
# Br6660_NAc7_3580   7776  7132  8150  3129     0
# Br6660_NAc8_4580   6938  6787   192  6645     0
# Br6660_NAc9_5080   7197  8130   361  2489     0
# Br6660_Nac10_4080  8968  7859   398  9671     0
# Br6660_Nac11_5580  3793  5557   206  3719   144

#The foil fell on tissue section 7 of Br6660. There is a cluster that seems to overlap with where the foil covered. 
#Cluster 23? 
spe$cluster_23 <- ifelse(spe$Banksy_louvain_res1.0 == 23,
                         "cluster_23",
                         "Other")

spe$cluster_23 <- as.factor(spe$cluster_23)

sub_spe <- spe[, spe$Sample == "Br6660_NAc7_3580"]

p <-make_escheR(spe[, spe$Sample == "Br6660_NAc7_3580"]) |>
  add_fill("cluster_23")+
  ggtitle(i) +
  theme(element_text(hjust = 0.5))

png(filename = here("plots","05_clustering","Banksy","Cluster23_Br6660_NAc7_3580.png"),
       units = "in",height = 20,width = 20,res = 300)
print(p)
dev.off()

#Are clusters differing by # of detected features or # of total counts? 
library(scuttle)
spe <- addPerCellQC(spe)

#Plot detected and sum
library(scater)

detected_vln <- plotColData(object = spe,
            x = "Banksy_louvain_res1.0",y = "detected") + 
  scale_y_log10() + 
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) 
ggsave(filename = here("plots","05_clustering","Banksy","Detected_Banksy_louvain_res1.0.png"),plot = detected_vln)


sum_vln <- plotColData(object = spe,
                            x = "Banksy_louvain_res1.0",y = "sum") + 
  scale_y_log10() + 
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) 
ggsave(filename = here("plots","05_clustering","Banksy","sum_Banksy_louvain_res1.0.png"),plot = sum_vln)

#Cluster 23 has relatively normal levels of genes + UMIs 


#What genes are expresssed in this cluster? 
markers_all <- c("SNAP25","PPP1R1B","BCL11B","GAD1",
                 "DRD1","PDYN","DRD2","PENK","ADORA2A",
                 "FOXP2","GABRQ","SEMA5B","CPNE4","VIP",
                 "CHAT","SLC5A7",
                 "SST","NPY","PNOC","CHODL",
                 "SLC17A7","TBR1",
                 "MOBP","ST18","OLIG1","OLIG2","OPALIN",
                 "SLC1A2","GFAP","GJA1","AQP4",
                 "CFAP157",
                 "C3","ARHGAP15",
                 "DCN","EBF1")

#Make a dotplot with some specific genes. 
library(scDotPlot)

assay(spe, "logcounts") <- assay(spe, "nucleus_normcounts")
p <- scDotPlot(object = spe,
          features = rev(markers_all),
          group = "Banksy_louvain_res1.0",
          groupAnno = "Banksy_louvain_res1.0",
          clusterColumns = TRUE,
          clusterRows = TRUE,
          scale = TRUE)

pdf(file = here("plots","05_clustering","GeneExpression","DotPlot_Markers_Subset.pdf"))
print(p)
dev.off()


#Set marker genes to be included on the heatmap.  
Panel <- readxl::read_excel(here("processed-data","Files_For_Upload","NAc_Xenium_Panel_Final_withNotes.xlsx"))
markers_all <- Panel$Gene[order(Panel$Cell_Type)]


p <- scDotPlot(object = spe,
               features = rev(markers_all),
               group = "Banksy_louvain_res1.0",
               groupAnno = "Banksy_louvain_res1.0",
               clusterColumns = TRUE,
               clusterRows = TRUE,
               scale = TRUE)

pdf(file = here("plots","05_clustering","GeneExpression","DotPlot_Markers_custompanel_100genes.pdf"),
    height = 20, width = 14)
print(p)
dev.off()


#Now let's do all of the genes. 
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]


p <- scDotPlot(object = spe,
               features = rownames(spe[gene_expression_idx,]),
               group = "Banksy_louvain_res1.0",
               groupAnno = "Banksy_louvain_res1.0",
               clusterColumns = TRUE,
               clusterRows = TRUE,
               scale = TRUE)

pdf(file = here("plots","05_clustering","GeneExpression","DotPlot_Markers_custompanel_ALLgenes.pdf"),
    height = 54, width = 14)
print(p)
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
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] scDotPlot_1.2.1             here_1.0.1                 
# [3] escheR_1.8.0                ggplot2_3.5.2              
# [5] sessioninfo_1.2.3           SpatialExperiment_1.18.1   
# [7] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [9] Biobase_2.68.0              GenomicRanges_1.60.0       
# [11] GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [13] S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [15] generics_0.1.4              MatrixGenerics_1.20.0      
# [17] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] RcppAnnoy_0.0.22        splines_4.5.0           later_1.4.2            
# [4] ggplotify_0.1.2         cellranger_1.1.0        tibble_3.2.1           
# [7] polyclip_1.10-7         fastDummies_1.7.5       lifecycle_1.0.4        
# [10] rprojroot_2.0.4         globals_0.18.0          lattice_0.22-7         
# [13] MASS_7.3-65             magrittr_2.0.3          plotly_4.10.4          
# [16] httpuv_1.6.16           Seurat_5.3.0            sctransform_0.4.2      
# [19] spam_2.11-1             sp_2.2-0                spatstat.sparse_3.1-0  
# [22] reticulate_1.42.0       cowplot_1.1.3           pbapply_1.7-2          
# [25] RColorBrewer_1.1-3      abind_1.4-8             Rtsne_0.17             
# [28] purrr_1.0.4             yulab.utils_0.2.0       GenomeInfoDbData_1.2.14
# [31] ggrepel_0.9.6           irlba_2.3.5.1           listenv_0.9.1          
# [34] spatstat.utils_3.1-4    tidytree_0.4.6          goftest_1.2-3          
# [37] RSpectra_0.16-2         spatstat.random_3.4-1   fitdistrplus_1.2-2     
# [40] parallelly_1.44.0       codetools_0.2-20        DelayedArray_0.34.1    
# [43] scuttle_1.18.0          tidyselect_1.2.1        aplot_0.2.7            
# [46] UCSC.utils_1.4.0        farver_2.1.2            viridis_0.6.5          
# [49] ScaledMatrix_1.16.0     spatstat.explore_3.4-3  jsonlite_2.0.0         
# [52] BiocNeighbors_2.2.0     progressr_0.15.1        ggridges_0.5.6         
# [55] survival_3.8-3          scater_1.36.0           tools_4.5.0            
# [58] treeio_1.32.0           ica_1.0-3               Rcpp_1.0.14            
# [61] glue_1.8.0              gridExtra_2.3           SparseArray_1.8.0      
# [64] dplyr_1.1.4             withr_3.0.2             fastmap_1.2.0          
# [67] digest_0.6.37           rsvd_1.0.5              R6_2.6.1               
# [70] mime_0.13               gridGraphics_0.5-1      colorspace_2.1-1       
# [73] scattermore_1.2         tensor_1.5              dichromat_2.0-0.1      
# [76] spatstat.data_3.1-6     utf8_1.2.5              tidyr_1.3.1            
# [79] ggsci_3.2.0             data.table_1.17.2       httr_1.4.7             
# [82] htmlwidgets_1.6.4       S4Arrays_1.8.0          scatterplot3d_0.3-44   
# [85] uwot_0.2.3              pkgconfig_2.0.3         gtable_0.3.6           
# [88] lmtest_0.9-40           XVector_0.48.0          htmltools_0.5.8.1      
# [91] dotCall64_1.2           SeuratObject_5.1.0      scales_1.4.0           
# [94] png_0.1-8               spatstat.univar_3.1-3   ggfun_0.1.8            
# [97] reshape2_1.4.4          rjson_0.2.23            nlme_3.1-168           
# [100] zoo_1.8-14              Polychrome_1.5.4        stringr_1.5.1          
# [103] KernSmooth_2.23-26      vipor_0.4.7             parallel_4.5.0         
# [106] miniUI_0.1.2            pillar_1.10.2           grid_4.5.0             
# [109] vctrs_0.6.5             RANN_2.6.2              promises_1.3.2         
# [112] BiocSingular_1.24.0     beachmat_2.24.0         xtable_1.8-4           
# [115] cluster_2.1.8.1         beeswarm_0.4.0          magick_2.8.6           
# [118] cli_3.6.5               compiler_4.5.0          rlang_1.1.6            
# [121] crayon_1.5.3            future.apply_1.11.3     labeling_0.4.3         
# [124] plyr_1.8.9              fs_1.6.6                ggbeeswarm_0.7.2       
# [127] stringi_1.8.7           viridisLite_0.4.2       deldir_2.0-4           
# [130] BiocParallel_1.42.0     lazyeval_0.2.2          spatstat.geom_3.4-1    
# [133] Matrix_1.7-3            RcppHNSW_0.6.0          patchwork_1.3.0        
# [136] future_1.49.0           shiny_1.10.0            ROCR_1.0-11            
# [139] igraph_2.1.4            ggtree_3.16.0           readxl_1.4.5           
# [142] ape_5.8-1     
# 
# 
