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

#Ependymal and Excitatory cells are primarily located within their own domains. The tests will not be meaningful for
# those cell types. Given that, let's remove them. 
sfe <- sfe[,sfe$snRNA_label %in% c("Astrocyte_A","Astrocyte_B",
                                   "DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                   "DRD2_MSN_A","DRD2_MSN_B","Endothelial","Inh_A",
                                   "Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                   "Microglia","Oligo","OPC")]

#Factorize columns for use 
sfe$spatial_0.4 <- as.factor(sfe$spatial_0.4)
table(is.na(sfe$spatial_0.4))
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

form <- ~ 0 + `spatial_0.4` + (1 | sample_id)

#Test between the putative island domains
L <- makeContrastsDream(form, colData(pb),
                        contrasts = c(D8vsD6 = "spatial_0.48- spatial_0.46"))

fit <- dream(cobj, form, colData(pb),L)

fit <- eBayes(fit)

topTable(fit, coef = "D8vsD6", number = Inf)
# logFC      AveExpr            t      P.Value    adj.P.Val
# DRD1_MSN_D   1.7650450284 -0.020548445  6.498282998 3.652484e-09 6.574472e-08
# DRD1_MSN_C  -1.3140024266  0.214166777 -4.743041003 7.349023e-06 6.614121e-05
# DRD2_MSN_B  -1.3357483468 -1.113272934 -4.534640554 1.670740e-05 1.002444e-04
# Inh_D        1.1330285834 -1.104409929  3.662159441 4.093768e-04 1.842195e-03
# Inh_A        0.5650033755 -0.818994918  2.578039527 1.145177e-02 3.606053e-02
# Inh_B       -0.7314466578 -0.792724558 -2.544820770 1.252891e-02 3.606053e-02
# Inh_F       -0.8062477493 -0.036026425 -2.502695698 1.402354e-02 3.606053e-02
# Microglia    0.5337309431  0.311144550  2.141425392 3.477879e-02 7.825228e-02
# OPC          0.5110499048 -0.176872972  2.074461931 4.072102e-02 8.144203e-02
# Oligo        0.4316603113  1.720100493  1.670344268 9.811647e-02 1.766097e-01
# Astrocyte_B -0.5151911175  1.179226178 -1.538728805 1.271668e-01 2.080912e-01
# DRD1_MSN_B  -0.1977067721 -0.006671803 -0.659926449 5.108860e-01 7.663289e-01
# Endothelial  0.1576385664  1.619274253  0.544026006 5.876883e-01 8.137223e-01
# Inh_C       -0.1131851857 -1.071534183 -0.426369836 6.707910e-01 8.624456e-01
# DRD1_MSN_A  -0.0421160478 -0.434518694 -0.121796494 9.033148e-01 9.991389e-01
# Astrocyte_A -0.0118592497  1.457860249 -0.047776493 9.619939e-01 9.991389e-01
# DRD2_MSN_A   0.0131886686 -0.542271450  0.034887154 9.722422e-01 9.991389e-01
# Inh_E        0.0002921975 -0.383926189  0.001082047 9.991389e-01 9.991389e-01
# B        z.std
# DRD1_MSN_D  10.5968331  5.899208953
# DRD1_MSN_C   3.2505940 -4.483321057
# DRD2_MSN_B   2.4617767 -4.304883073
# Inh_D       -0.5739607  3.533964042
# Inh_A       -3.6015746  2.528602440
# Inh_B       -3.7289612 -2.496886331
# Inh_F       -3.8230424 -2.456659906
# Microglia   -4.5664221  2.110924614
# OPC         -4.7315557  2.046359483
# Oligo       -5.4232608  1.654054335
# Astrocyte_B -5.6128912 -1.525370069
# DRD1_MSN_B  -6.6195004 -0.657458806
# Endothelial -6.6488725  0.542188996
# Inh_C       -6.7266808 -0.425062609
# DRD1_MSN_A  -6.7875115 -0.121475012
# Astrocyte_A -6.7839604 -0.047651659
# DRD2_MSN_A  -6.8049749  0.034796228
# Inh_E       -6.8314270  0.001079235


#Build the cluster tree
hcl <- buildClusterTreeFromPB(pb)

res <- treeTest(fit, cobj, hcl, coef = "D8vsD6")

forest_plot <- plotForest(res)


pdf(file = here("plots","HD_Full_Analysis","crumblr","Forest_plot_D8vsD6.pdf"),height =8, width = 8)
forest_plot
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
# [11] backports_1.5.0           R6_2.6.1                 
# [13] metafor_5.0-1             lazyeval_0.2.2           
# [15] rhdf5filters_1.20.0       withr_3.0.2              
# [17] sp_2.2-0                  gridExtra_2.3            
# [19] prettyunits_1.2.0         cli_3.6.5                
# [21] sandwich_3.1-1            labeling_0.4.3           
# [23] KEGGgraph_1.70.0          SQUAREM_2021.1           
# [25] mvtnorm_1.3-3             proxy_0.4-27             
# [27] mixsqp_0.3-54             yulab.utils_0.2.4        
# [29] R.utils_2.13.0            zenith_1.12.0            
# [31] dichromat_2.0-0.1         invgamma_1.1             
# [33] RSQLite_2.3.11            gridGraphics_0.5-1       
# [35] gtools_3.9.5              spdep_1.3-11             
# [37] dplyr_1.1.4               metadat_1.6-0            
# [39] ggbeeswarm_0.7.2          R.methodsS3_1.8.2        
# [41] terra_1.8-50              lifecycle_1.0.4          
# [43] multcomp_1.4-28           edgeR_4.6.2              
# [45] mathjaxr_2.0-0            gplots_3.2.0             
# [47] grid_4.5.0                blob_1.2.4               
# [49] dqrng_0.4.1               crayon_1.5.3             
# [51] lattice_0.22-7            beachmat_2.24.0          
# [53] msigdbr_25.1.1            annotate_1.86.0          
# [55] KEGGREST_1.48.0           magick_2.8.6             
# [57] zeallot_0.1.0             pillar_1.10.2            
# [59] rjson_0.2.23              boot_1.3-31              
# [61] corpcor_1.6.10            codetools_0.2-20         
# [63] wk_0.9.4                  glue_1.8.0               
# [65] ggfun_0.1.8               data.table_1.17.2        
# [67] treeio_1.32.0             vctrs_0.6.5              
# [69] png_0.1-8                 Rdpack_2.6.4             
# [71] gtable_0.3.6              assertthat_0.2.1         
# [73] cachem_1.1.0              zigg_0.0.2               
# [75] rbibutils_2.3             DropletUtils_1.28.0      
# [77] Rfast_2.1.5.2             coda_0.19-4.1            
# [79] reformulas_0.4.4          survival_3.8-3           
# [81] sfheaders_0.4.4           iterators_1.0.14         
# [83] units_0.8-7               statmod_1.5.0            
# [85] TH.data_1.1-3             nlme_3.1-168             
# [87] pbkrtest_0.5.4            ggtree_3.16.0            
# [89] bit64_4.6.0-1             progress_1.2.3           
# [91] EnvStats_3.1.0            rprojroot_2.0.4          
# [93] irlba_2.3.5.1             vipor_0.4.7              
# [95] KernSmooth_2.23-26        spData_2.3.4             
# [97] rmeta_3.0                 DBI_1.2.3                
# [99] tidyselect_1.2.1          bit_4.6.0                
# [101] compiler_4.5.0            curl_6.2.2               
# [103] graph_1.86.0              BiocNeighbors_2.2.0      
# [105] scales_1.4.0              caTools_1.18.3           
# [107] classInt_0.4-11           remaCor_0.0.20           
# [109] rappdirs_0.3.3            tiff_0.1-12              
# [111] stringr_1.5.1             digest_0.6.37            
# [113] fftwtools_0.9-11          minqa_1.2.8              
# [115] aod_1.3.3                 XVector_0.48.0           
# [117] RhpcBLASctl_0.23-42       htmltools_0.5.8.1        
# [119] pkgconfig_2.0.3           jpeg_0.1-11              
# [121] lme4_2.0-1                sparseMatrixStats_1.20.0 
# [123] mashr_0.2.79              fastmap_1.2.0            
# [125] rlang_1.1.6               htmlwidgets_1.6.4        
# [127] UCSC.utils_1.4.0          DelayedMatrixStats_1.30.0
# [129] farver_2.1.2              zoo_1.8-14               
# [131] jsonlite_2.0.0            R.oo_1.27.1              
# [133] RCurl_1.98-1.17           magrittr_2.0.3           
# [135] ggplotify_0.1.2           scuttle_1.18.0           
# [137] GenomeInfoDbData_1.2.14   s2_1.1.8                 
# [139] patchwork_1.3.0           Rhdf5lib_1.30.0          
# [141] Rcpp_1.1.1-1.1            viridis_0.6.5            
# [143] ape_5.8-1                 babelgene_22.9           
# [145] EnrichmentBrowser_2.40.0  stringi_1.8.7            
# [147] MASS_7.3-65               plyr_1.8.9               
# [149] parallel_4.5.0            ggrepel_0.9.6            
# [151] deldir_2.0-4              Biostrings_2.76.0        
# [153] splines_4.5.0             hms_1.1.3                
# [155] locfit_1.5-9.12           reshape2_1.4.4           
# [157] LearnBayes_2.15.1         XML_3.99-0.18            
# [159] RcppParallel_5.1.10       nloptr_2.2.1             
# [161] tidyr_1.3.1               purrr_1.0.4              
# [163] scattermore_1.2           ashr_2.2-63              
# [165] broom_1.0.8               xtable_1.8-4             
# [167] tidytree_0.4.6            fANCOVA_0.6-1            
# [169] e1071_1.7-16              viridisLite_0.4.2        
# [171] class_7.3-23              truncnorm_1.0-9          
# [173] tibble_3.2.1              aplot_0.2.7              
# [175] lmerTest_3.2-1            memoise_2.0.1            
# [177] beeswarm_0.4.0            AnnotationDbi_1.70.0     
# [179] GSEABase_1.70.0          
