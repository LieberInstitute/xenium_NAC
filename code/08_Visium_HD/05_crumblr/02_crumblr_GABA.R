# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(HDF5Array)
library(dreamlet)
library(crumblr)
library(here)

# Load object 
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_annotated_labels")

sfe <- loadHDF5SummarizedExperiment(sfe_dir)

#Ependymal and Excitatory cells are primarily located within their own domains. The tests will not be meaningful for
# those cell types. Given that, let's remove them. 
sfe <- sfe[,sfe$snRNA_label %in% c("Astrocyte_A","Astrocyte_B",
                                   "DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                   "DRD2_MSN_A","DRD2_MSN_B","Endothelial","Inh_A",
                                   "Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                   "Microglia","Oligo","OPC")]

#Factorize columns for use 
sfe$spatial_domain <- as.factor(sfe$spatial_domain)
table(is.na(sfe$spatial_domain))
# FALSE 
# 426497 

sfe$Sample <- paste0(sfe$sample_id, "_", sfe$spatial_0.4)
sfe$Sample <- as.factor(sfe$Sample)
table(is.na(sfe$Sample))
# FALSE 
# 426497 

sfe$sample_id <- as.factor(sfe$sample_id)
table(is.na(sfe$sample_id))
# FALSE 
# 426497 

sfe$snRNA_label <- as.factor(sfe$snRNA_label)
table(is.na(sfe$snRNA_label))
# FALSE 
# 426497 

# Pseudobulk: cluster by cell type, aggregate by sample
# Warning messages coming from internal package calls. Nothing to worry about. 
pb <- aggregateToPseudoBulk(sfe,
                            assay    = "counts",
                            cluster_id = "snRNA_label",
                            sample_id  = "Sample",
                            verbose  = TRUE)
# Processing block [[2/2, 39/39]] ... OKWarning messages:
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

# spatial_0.4 is the fixed effect --> 12-level domain factor. 
# Mean CLR value for that domain, which you subtract from one another to compare domains. 

# (1 | sample_id) is  random intercept grouped by sample. 
#  Each of your 8 samples gets its own baseline offset, modeled as a draw from a normal distribution whose variance is estimated. 
# This soaks up each donor's overall composition quirks so the domain comparison is effectively within-sample

form <- ~ 0 + spatial_domain + (1 | sample_id)

#Test between the putative island domains
L <- makeContrastsDream(form, colData(pb),
                        contrasts = c(AvsB = "spatial_domainGABA_A- spatial_domainGABA_B",
                                      AvsC = "spatial_domainGABA_A- spatial_domainGABA_C",
                                      BvsC = "spatial_domainGABA_B- spatial_domainGABA_C"))

fit <- dream(cobj, form, colData(pb),L)

fit <- eBayes(fit)


#Build the cluster tree
hcl <- buildClusterTreeFromPB(pb)

resAvB <- treeTest(fit, cobj, hcl, coef = "AvsB")
resAvC <- treeTest(fit, cobj, hcl, coef = "AvsC")
resBvC <- treeTest(fit, cobj, hcl, coef = "BvsC")

pdf(file = here("plots","HD_Full_Analysis","crumblr","Forest_plot_GABAAvsGABAB.pdf"),height =8, width = 8)
plotForest(resAvB)
dev.off()

pdf(file = here("plots","HD_Full_Analysis","crumblr","Forest_plot_GABAAvsGABAC.pdf"),height =8, width = 8)
plotForest(resAvC)
dev.off()


pdf(file = here("plots","HD_Full_Analysis","crumblr","Forest_plot_GABABvsGABAC.pdf"),height =8, width = 8)
plotForest(resBvC)
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
#   [1] here_1.0.1                      crumblr_1.0.0                  
# [3] dreamlet_1.9.2                  variancePartition_1.40.2       
# [5] BiocParallel_1.42.0             limma_3.64.3                   
# [7] ggplot2_3.5.2                   HDF5Array_1.36.0               
# [9] h5mread_1.0.0                   rhdf5_2.52.0                   
# [11] DelayedArray_0.34.1             SparseArray_1.8.0              
# [13] S4Arrays_1.8.0                  abind_1.4-8                    
# [15] Matrix_1.7-3                    sessioninfo_1.2.3              
# [17] ComplexHeatmap_2.24.0           SpatialExperiment_1.18.1       
# [19] SingleCellExperiment_1.30.1     SummarizedExperiment_1.38.1    
# [21] Biobase_2.68.0                  GenomicRanges_1.60.0           
# [23] GenomeInfoDb_1.44.0             IRanges_2.42.0                 
# [25] S4Vectors_0.48.1                BiocGenerics_0.54.0            
# [27] generics_0.1.4                  MatrixGenerics_1.20.0          
# [29] matrixStats_1.5.0               SpatialFeatureExperiment_1.10.1
# 
# loaded via a namespace (and not attached):
#   [1] fs_1.6.6                  spatialreg_1.3-6         
# [3] bitops_1.0-9              sf_1.0-21                
# [5] EBImage_4.50.0            httr_1.4.7               
# [7] RColorBrewer_1.1-3        doParallel_1.0.17        
# [9] Rgraphviz_2.52.0          numDeriv_2016.8-1.1      
# [11] tools_4.5.0               backports_1.5.0          
# [13] R6_2.6.1                  metafor_5.0-1            
# [15] lazyeval_0.2.2            rhdf5filters_1.20.0      
# [17] GetoptLong_1.0.5          withr_3.0.2              
# [19] sp_2.2-0                  gridExtra_2.3            
# [21] prettyunits_1.2.0         cli_3.6.5                
# [23] sandwich_3.1-1            labeling_0.4.3           
# [25] KEGGgraph_1.70.0          SQUAREM_2021.1           
# [27] mvtnorm_1.3-3             proxy_0.4-27             
# [29] mixsqp_0.3-54             yulab.utils_0.2.4        
# [31] R.utils_2.13.0            zenith_1.12.0            
# [33] dichromat_2.0-0.1         invgamma_1.1             
# [35] RSQLite_2.3.11            gridGraphics_0.5-1       
# [37] shape_1.4.6.1             gtools_3.9.5             
# [39] spdep_1.3-11              dplyr_1.1.4              
# [41] metadat_1.6-0             ggbeeswarm_0.7.2         
# [43] R.methodsS3_1.8.2         terra_1.8-50             
# [45] lifecycle_1.0.4           multcomp_1.4-28          
# [47] edgeR_4.6.2               mathjaxr_2.0-0           
# [49] gplots_3.2.0              blob_1.2.4               
# [51] dqrng_0.4.1               crayon_1.5.3             
# [53] lattice_0.22-7            msigdbr_25.1.1           
# [55] beachmat_2.24.0           annotate_1.86.0          
# [57] KEGGREST_1.48.0           magick_2.8.6             
# [59] zeallot_0.1.0             pillar_1.10.2            
# [61] rjson_0.2.23              boot_1.3-31              
# [63] corpcor_1.6.10            codetools_0.2-20         
# [65] wk_0.9.4                  glue_1.8.0               
# [67] ggfun_0.1.8               data.table_1.17.2        
# [69] treeio_1.32.0             vctrs_0.6.5              
# [71] png_0.1-8                 Rdpack_2.6.4             
# [73] gtable_0.3.6              assertthat_0.2.1         
# [75] cachem_1.1.0              zigg_0.0.2               
# [77] rbibutils_2.3             DropletUtils_1.28.0      
# [79] Rfast_2.1.5.2             coda_0.19-4.1            
# [81] reformulas_0.4.4          survival_3.8-3           
# [83] sfheaders_0.4.4           iterators_1.0.14         
# [85] units_0.8-7               statmod_1.5.0            
# [87] TH.data_1.1-3             nlme_3.1-168             
# [89] pbkrtest_0.5.4            ggtree_3.16.0            
# [91] bit64_4.6.0-1             progress_1.2.3           
# [93] EnvStats_3.1.0            rprojroot_2.0.4          
# [95] irlba_2.3.5.1             vipor_0.4.7              
# [97] KernSmooth_2.23-26        colorspace_2.1-1         
# [99] spData_2.3.4              rmeta_3.0                
# [101] DBI_1.2.3                 tidyselect_1.2.1         
# [103] curl_6.2.2                bit_4.6.0                
# [105] compiler_4.5.0            graph_1.86.0             
# [107] BiocNeighbors_2.2.0       scales_1.4.0             
# [109] caTools_1.18.3            classInt_0.4-11          
# [111] remaCor_0.0.20            rappdirs_0.3.3           
# [113] tiff_0.1-12               stringr_1.5.1            
# [115] digest_0.6.37             fftwtools_0.9-11         
# [117] minqa_1.2.8               aod_1.3.3                
# [119] XVector_0.48.0            RhpcBLASctl_0.23-42      
# [121] htmltools_0.5.8.1         pkgconfig_2.0.3          
# [123] jpeg_0.1-11               lme4_2.0-1               
# [125] sparseMatrixStats_1.20.0  mashr_0.2.79             
# [127] fastmap_1.2.0             rlang_1.1.6              
# [129] GlobalOptions_0.1.2       htmlwidgets_1.6.4        
# [131] UCSC.utils_1.4.0          DelayedMatrixStats_1.30.0
# [133] farver_2.1.2              zoo_1.8-14               
# [135] jsonlite_2.0.0            R.oo_1.27.1              
# [137] RCurl_1.98-1.17           magrittr_2.0.3           
# [139] ggplotify_0.1.2           scuttle_1.18.0           
# [141] GenomeInfoDbData_1.2.14   s2_1.1.8                 
# [143] patchwork_1.3.0           Rhdf5lib_1.30.0          
# [145] Rcpp_1.1.1-1.1            viridis_0.6.5            
# [147] ape_5.8-1                 babelgene_22.9           
# [149] EnrichmentBrowser_2.40.0  stringi_1.8.7            
# [151] MASS_7.3-65               plyr_1.8.9               
# [153] parallel_4.5.0            ggrepel_0.9.6            
# [155] deldir_2.0-4              Biostrings_2.76.0        
# [157] splines_4.5.0             hms_1.1.3                
# [159] circlize_0.4.16           locfit_1.5-9.12          
# [161] reshape2_1.4.4            LearnBayes_2.15.1        
# [163] XML_3.99-0.18             RcppParallel_5.1.10      
# [165] nloptr_2.2.1              foreach_1.5.2            
# [167] tidyr_1.3.1               purrr_1.0.4              
# [169] clue_0.3-66               scattermore_1.2          
# [171] ashr_2.2-63               broom_1.0.8              
# [173] xtable_1.8-4              tidytree_0.4.6           
# [175] fANCOVA_0.6-1             e1071_1.7-16             
# [177] viridisLite_0.4.2         class_7.3-23             
# [179] truncnorm_1.0-9           tibble_3.2.1             
# [181] aplot_0.2.7               lmerTest_3.2-1           
# [183] memoise_2.0.1             beeswarm_0.4.0           
# [185] AnnotationDbi_1.70.0      cluster_2.1.8.1          
# [187] GSEABase_1.70.0          
