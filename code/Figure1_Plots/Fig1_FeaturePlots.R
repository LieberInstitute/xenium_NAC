# /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# module load conda_R/4.5

library(SpatialExperiment)
library(escheR)
library(here)

#Read in the spe object
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))


#Read in the colors
CellType_cols <- readRDS(here("processed-data","05_Clustering","CellType_cols.Rds"))

CellType_cols
# Excitatory  Microglia_A         WM_A         WM_B         WM_C  D1_Island_B 
# "#DB1C5F"    "#0D87E4"    "#FDCB22"    "#F3700D"    "#385A51"    "#00FF0D" 
# Astro_A Fibroblast_A Fibroblast_B     DRD2_MSN     DRD1_MSN      Astro_B 
# "#FD00FD"    "#FFA7E2"    "#2AFECA"    "#7AAA16"    "#9400FF"    "#823526" 
# MSN_Oligo    Inh_PVALB          OPC  Microglia_B      Inh_SST         WM_D 
# "#F5DEC0"    "#93E5FF"    "#7A1699"    "#FF0DBC"    "#C4B3FB"    "#F80D2A" 
# D1_Island_A         CHAT    Ependymal         WM_E 
# "#224B82"    "#FBA475"    "#73690D"    "#B62A7A" 

#Rename colors 
names(CellType_cols)[18] <- "Astrocyte_Oligo"
names(CellType_cols)[22] <- "Microglia_Oligo"

#Remove WM_A-C colors
CellType_cols <- CellType_cols[-c(3:5)]

#Add a color for WM
CellType_cols <- c(CellType_cols,"orange")
names(CellType_cols)[20] <- "WM"

saveRDS(object = CellType_cols,file = here("processed-data","CellType_Cols_v2.Rds"))

#Plot a single depth so that we can pull the legend fror banksy colors
spe_sub <- spe[,spe$Sample == "Br6660_Nac10_4080"]

p1 <- make_escheR(spe_sub) |>
  add_fill("CellTypes") +
  scale_fill_manual(values = CellType_cols)

p1_legend <- cowplot::get_legend(p1)
ggsave(p1_legend,file = here("plots","05_clustering","Banksy_cellTypes_v2_legend.pdf"))



#Subset to just Br6660
spe_sub2 <- spe[,spe$Donor == "Br6660"]

spe_sub2$MOBP <- assay(spe_sub2,"nucleus_normcounts")["MOBP",]
spe_sub2$PPP1R1B <- assay(spe_sub2,"nucleus_normcounts")["PPP1R1B",]
spe_sub2$DRD1 <- assay(spe_sub2,"nucleus_normcounts")["DRD1",]
spe_sub2$DRD2 <- assay(spe_sub2,"nucleus_normcounts")["DRD2",]

#Make feature plots of MOBP,PPP1R1B, DRD1,DRD2
genes <- c("MOBP","PPP1R1B","DRD1","DRD2")
for(sample in unique(spe_sub2$Sample)){
  for(gene in genes){
    spe_sub3 <- spe_sub2[,spe_sub2$Sample == sample]
    p <- make_escheR(spe_sub3) |>
      add_fill(gene) +
      scale_fill_gradientn(colors = c("lightgrey","orange","red"))
    ggsave(p,
           filename = here("plots","Figure1","FeaturePlots",paste0(sample,"_",gene,".png")),
           height= 18, width = 20)
  }
}
  
#Save the legend for each gene
for(gene in genes){
  print(gene)
  p2 <- make_escheR(spe_sub3) |>
    add_fill(gene) +
    scale_fill_gradientn(colors = c("lightgrey","orange","red"))
  
  p2_legend <- cowplot::get_legend(p2)
  ggsave(p2_legend,file = here("plots","Figure1",paste0(gene,"_","FeaturePlotLegend.pdf")))
}

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
#   [1] here_1.0.1                  escheR_1.8.0               
# [3] ggplot2_3.5.2               SpatialExperiment_1.18.1   
# [5] SingleCellExperiment_1.30.1 SummarizedExperiment_1.38.1
# [7] Biobase_2.68.0              GenomicRanges_1.60.0       
# [9] GenomeInfoDb_1.44.0         IRanges_2.42.0             
# [11] S4Vectors_0.46.0            BiocGenerics_0.54.0        
# [13] generics_0.1.4              MatrixGenerics_1.20.0      
# [15] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] SparseArray_1.8.0       lattice_0.22-7          magrittr_2.0.3         
# [4] grid_4.5.0              RColorBrewer_1.1-3      rprojroot_2.0.4        
# [7] jsonlite_2.0.0          Matrix_1.7-3            httr_1.4.7             
# [10] viridisLite_0.4.2       UCSC.utils_1.4.0        scales_1.4.0           
# [13] textshaping_1.0.1       abind_1.4-8             cli_3.6.5              
# [16] rlang_1.1.6             crayon_1.5.3            XVector_0.48.0         
# [19] cowplot_1.1.3           withr_3.0.2             DelayedArray_0.34.1    
# [22] S4Arrays_1.8.0          tools_4.5.0             dplyr_1.1.4            
# [25] GenomeInfoDbData_1.2.14 vctrs_0.6.5             R6_2.6.1               
# [28] lifecycle_1.0.4         magick_2.8.6            ragg_1.4.0             
# [31] pkgconfig_2.0.3         pillar_1.10.2           gtable_0.3.6           
# [34] glue_1.8.0              Rcpp_1.0.14             systemfonts_1.2.3      
# [37] tibble_3.2.1            tidyselect_1.2.1        dichromat_2.0-0.1      
# [40] farver_2.1.2            rjson_0.2.23            labeling_0.4.3         
# [43] compiler_4.5.0     
  
  
  
