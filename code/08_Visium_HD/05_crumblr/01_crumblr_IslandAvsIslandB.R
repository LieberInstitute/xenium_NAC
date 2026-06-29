# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(dreamlet)
library(crumblr)
library(here)

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


sfe <- sfe[,sfe$snRNA_label %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D")]

#Factorize columns for use 
sfe$spatial_0.4 <- as.factor(sfe$spatial_0.4)
table(is.na(sfe$spatial_0.4))
# FALSE 
# 80760

sfe$Spatial_Domain <- as.factor(sfe$Spatial_Domain)
table(is.na(sfe$Spatial_Domain))
# FALSE 
# 80760

sfe$Sample <- paste0(sfe$sample_id, "_", sfe$Spatial_Domain)
sfe$Sample <- as.factor(sfe$Sample)
table(is.na(sfe$Sample))
# FALSE 
# 80760

sfe$sample_id <- as.factor(sfe$sample_id)
table(is.na(sfe$sample_id))
# FALSE 
# 80760 

sfe$snRNA_label <- as.factor(sfe$snRNA_label)
table(is.na(sfe$snRNA_label))
# FALSE 
# 80760 

# Pseudobulk: cluster by cell type, aggregate by sample
# Warning messages coming from internal package calls. Nothing to worry about. 
pb <- aggregateToPseudoBulk(sfe,
                            assay    = "counts",
                            cluster_id = "snRNA_label",
                            sample_id  = "Sample",
                            verbose  = TRUE)
# Processing block [[2/2, 8/8]] ... OKWarning messages:
#   1: In S4Vectors:::anyMissing(runValue(x_seqnames)) :
#   'S4Vectors:::anyMissing()' is deprecated.
# Use 'anyNA()' instead.
# See help("Deprecated")
# 2: In S4Vectors:::anyMissing(runValue(strand(x))) :
#   'S4Vectors:::anyMissing()' is deprecated.
# Use 'anyNA()' instead.
# See help("Deprecated")

cobj <- crumblr(cellCounts(pb))

################################################
# 0 + removes that global intercept and get 12 standalone coefficients, 
# one per domain, each equal to that domain's own mean.

# Spatial_Domain is the fixed effect --> 12-level domain factor. 
# Mean CLR value for that domain, which you subtract from one another to compare domains. 

# (1 | sample_id) is  random intercept grouped by sample. 
#  Each of your 8 samples gets its own baseline offset, modeled as a draw from a normal distribution whose variance is estimated. 
# This soaks up each donor's overall composition quirks so the domain comparison is effectively within-sample

form <- ~ 0 + `Spatial_Domain` + (1 | sample_id)

#Test between the putative island domains
L <- makeContrastsDream(form, colData(pb),
                        contrasts = c(IslandAvsIsland_B = "Spatial_DomainD1_Island_A - Spatial_DomainD1_Island_B"))

fit <- dream(cobj, form, colData(pb),L)

fit <- eBayes(fit)

topTable(fit, coef = "IslandAvsIsland_B", number = Inf)
# logFC     AveExpr          t      P.Value    adj.P.Val
# DRD1_MSN_C -1.40164950  0.30932004 -6.8565511 4.039689e-10 1.615876e-09
# DRD1_MSN_D  1.69172482  0.04632587  6.6719327 1.003356e-09 2.006712e-09
# DRD1_MSN_B -0.25927114  0.06187440 -1.1912520 2.360676e-01 3.147568e-01
# DRD1_MSN_A -0.09946615 -0.41752031 -0.4458111 6.565915e-01 6.565915e-01
# B     z.std
# DRD1_MSN_C 12.671905 -6.252486
# DRD1_MSN_D 11.795391  6.108875
# DRD1_MSN_B -6.193438 -1.184873
# DRD1_MSN_A -6.764217 -0.444624


#Build the cluster tree
hcl <- buildClusterTreeFromPB(pb)

res <- treeTest(fit, cobj, hcl, coef = "IslandAvsIsland_B")

forest_plot <- plotForest(res)


pdf(file = here("plots","HD_Full_Analysis","crumblr","Forest_plot_IslandAvsIsland_B.pdf"),height =8, width = 8)
forest_plot
dev.off()

pdf(file = here("plots","HD_Full_Analysis","crumblr","plotTreeTest_IslandAvsIsland_B.pdf"),height =8, width = 8)
plotTreeTest(res)
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
#   [1] here_1.0.1                      crumblr_1.0.0                  
# [3] dreamlet_1.9.2                  variancePartition_1.40.2       
# [5] BiocParallel_1.42.0             limma_3.64.3                   
# [7] ggplot2_3.5.2                   HDF5Array_1.36.0               
# [9] h5mread_1.0.0                   rhdf5_2.52.0                   
# [11] DelayedArray_0.34.1             SparseArray_1.8.0              
# [13] S4Arrays_1.8.0                  abind_1.4-8                    
# [15] Matrix_1.7-3                    sessioninfo_1.2.3              
# [17] SpatialExperiment_1.18.1        SingleCellExperiment_1.30.1    
# [19] SummarizedExperiment_1.38.1     Biobase_2.68.0                 
# [21] GenomicRanges_1.60.0            GenomeInfoDb_1.44.0            
# [23] IRanges_2.42.0                  S4Vectors_0.48.1               
# [25] BiocGenerics_0.54.0             generics_0.1.4                 
# [27] MatrixGenerics_1.20.0           matrixStats_1.5.0              
# [29] SpatialFeatureExperiment_1.10.1
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.6                  spatialreg_1.3-6         
# [3] bitops_1.0-9              sf_1.0-21                
# [5] EBImage_4.50.0            httr_1.4.7               
# [7] RColorBrewer_1.1-3        Rgraphviz_2.52.0         
# [9] numDeriv_2016.8-1.1       tools_4.5.0              
# [11] backports_1.5.0           utf8_1.2.5               
# [13] R6_2.6.1                  metafor_5.0-1            
# [15] lazyeval_0.2.2            rhdf5filters_1.20.0      
# [17] withr_3.0.2               sp_2.2-0                 
# [19] gridExtra_2.3             prettyunits_1.2.0        
# [21] cli_3.6.5                 sandwich_3.1-1           
# [23] labeling_0.4.3            KEGGgraph_1.70.0         
# [25] SQUAREM_2021.1            mvtnorm_1.3-3            
# [27] proxy_0.4-27              mixsqp_0.3-54            
# [29] yulab.utils_0.2.4         R.utils_2.13.0           
# [31] zenith_1.12.0             dichromat_2.0-0.1        
# [33] invgamma_1.1              RSQLite_2.3.11           
# [35] gridGraphics_0.5-1        gtools_3.9.5             
# [37] spdep_1.3-11              dplyr_1.1.4              
# [39] metadat_1.6-0             ggbeeswarm_0.7.2         
# [41] R.methodsS3_1.8.2         terra_1.8-50             
# [43] lifecycle_1.0.4           multcomp_1.4-28          
# [45] edgeR_4.6.2               mathjaxr_2.0-0           
# [47] gplots_3.2.0              grid_4.5.0               
# [49] blob_1.2.4                dqrng_0.4.1              
# [51] crayon_1.5.3              lattice_0.22-7           
# [53] beachmat_2.24.0           msigdbr_25.1.1           
# [55] annotate_1.86.0           KEGGREST_1.48.0          
# [57] magick_2.8.6              zeallot_0.1.0            
# [59] pillar_1.10.2             rjson_0.2.23             
# [61] boot_1.3-31               corpcor_1.6.10           
# [63] codetools_0.2-20          wk_0.9.4                 
# [65] glue_1.8.0                ggfun_0.1.8              
# [67] data.table_1.17.2         treeio_1.32.0            
# [69] vctrs_0.6.5               png_0.1-8                
# [71] Rdpack_2.6.4              gtable_0.3.6             
# [73] assertthat_0.2.1          cachem_1.1.0             
# [75] zigg_0.0.2                rbibutils_2.3            
# [77] DropletUtils_1.28.0       Rfast_2.1.5.2            
# [79] coda_0.19-4.1             reformulas_0.4.4         
# [81] survival_3.8-3            sfheaders_0.4.4          
# [83] iterators_1.0.14          units_0.8-7              
# [85] statmod_1.5.0             TH.data_1.1-3            
# [87] nlme_3.1-168              pbkrtest_0.5.4           
# [89] ggtree_3.16.0             bit64_4.6.0-1            
# [91] progress_1.2.3            EnvStats_3.1.0           
# [93] rprojroot_2.0.4           irlba_2.3.5.1            
# [95] vipor_0.4.7               KernSmooth_2.23-26       
# [97] spData_2.3.4              rmeta_3.0                
# [99] DBI_1.2.3                 tidyselect_1.2.1         
# [101] bit_4.6.0                 compiler_4.5.0           
# [103] curl_6.2.2                graph_1.86.0             
# [105] BiocNeighbors_2.2.0       scales_1.4.0             
# [107] caTools_1.18.3            classInt_0.4-11          
# [109] remaCor_0.0.20            rappdirs_0.3.3           
# [111] tiff_0.1-12               stringr_1.5.1            
# [113] digest_0.6.37             fftwtools_0.9-11         
# [115] minqa_1.2.8               aod_1.3.3                
# [117] XVector_0.48.0            RhpcBLASctl_0.23-42      
# [119] htmltools_0.5.8.1         pkgconfig_2.0.3          
# [121] jpeg_0.1-11               lme4_2.0-1               
# [123] sparseMatrixStats_1.20.0  mashr_0.2.79             
# [125] fastmap_1.2.0             rlang_1.1.6              
# [127] htmlwidgets_1.6.4         UCSC.utils_1.4.0         
# [129] DelayedMatrixStats_1.30.0 farver_2.1.2             
# [131] zoo_1.8-14                jsonlite_2.0.0           
# [133] R.oo_1.27.1               RCurl_1.98-1.17          
# [135] magrittr_2.0.3            ggplotify_0.1.2          
# [137] scuttle_1.18.0            GenomeInfoDbData_1.2.14  
# [139] s2_1.1.8                  patchwork_1.3.0          
# [141] Rhdf5lib_1.30.0           Rcpp_1.1.1-1.1           
# [143] viridis_0.6.5             ape_5.8-1                
# [145] babelgene_22.9            EnrichmentBrowser_2.40.0 
# [147] stringi_1.8.7             MASS_7.3-65              
# [149] plyr_1.8.9                parallel_4.5.0           
# [151] ggrepel_0.9.6             deldir_2.0-4             
# [153] Biostrings_2.76.0         splines_4.5.0            
# [155] hms_1.1.3                 locfit_1.5-9.12          
# [157] reshape2_1.4.4            LearnBayes_2.15.1        
# [159] XML_3.99-0.18             RcppParallel_5.1.10      
# [161] nloptr_2.2.1              tidyr_1.3.1              
# [163] purrr_1.0.4               scattermore_1.2          
# [165] ashr_2.2-63               broom_1.0.8              
# [167] xtable_1.8-4              tidytree_0.4.6           
# [169] fANCOVA_0.6-1             e1071_1.7-16             
# [171] viridisLite_0.4.2         class_7.3-23             
# [173] truncnorm_1.0-9           tibble_3.2.1             
# [175] aplot_0.2.7               lmerTest_3.2-1           
# [177] memoise_2.0.1             beeswarm_0.4.0           
# [179] AnnotationDbi_1.70.0      GSEABase_1.70.0          
