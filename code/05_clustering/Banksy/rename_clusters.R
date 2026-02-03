# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(escheR)
library(here)

#Load the xenium object 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

#Add the Banksy cell types to the domain. 
#Load the annotated banksy cell types
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

stopifnot(identical(Banksy_celltypes$cell_id,colnames(spe)))

spe$Banksy_celltypes <- Banksy_celltypes$CellType


spe
# class: SpatialExperiment 
# dim: 366 4849373 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4849373): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(60): Sample Barcode ... nucleus_area.sf Banksy_celltypes
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Subset for just white matter cell types
spe_wm <- spe[,spe$Banksy_celltypes %in% c("WM_A","WM_B","WM_C",
                                           "WM_D","WM_E")]

#Set the logcounts
logcounts(spe_wm) <- assay(spe_wm,"nucleus_normcounts")

mod <- with(colData(spe_wm), model.matrix(~ Sample))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- scran::findMarkers(spe_wm, 
                                       groups=spe_wm$Banksy_celltypes,
                                       assay.type="logcounts", 
                                       design=mod, 
                                       test="t",
                                       direction="up", 
                                       pval.type="all", 
                                       full.stats=T)

saveRDS(object = markers_pairwise,
        file = here("processed-data","WM_pairwise_DEGs.Rds"))



#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})
# $WM_A
# 
# FALSE  TRUE 
# 363     3 
# 
# $WM_B
# 
# FALSE 
# 366 
# 
# $WM_C
# 
# FALSE  TRUE 
# 287    79 
# 
# $WM_D
# 
# FALSE  TRUE 
# 221   145 
# 
# $WM_E
# 
# FALSE  TRUE 
# 275    91 

#WM_A,B,C seems to be WM based on expression and anatomy
#WM_D is Astrocyte_Oligo
#WM_E is Microglia_Oligo

#rename the cell types
spe$CellTypes <- spe$Banksy_celltypes

#Rename the WM celltypes 
spe[,spe$Banksy_celltypes %in% c("WM_A","WM_B","WM_C")]$CellTypes <- "WM"
spe[,spe$Banksy_celltypes %in% c("WM_D")]$CellTypes <- "Astrocyte_Oligo"
spe[,spe$Banksy_celltypes %in% c("WM_E")]$CellTypes <- "Microglia_Oligo"

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

for(sample in unique(spe$Sample)){
  print(sample)
  sub_spe <- spe[,spe$Sample == sample]
  p <- make_escheR(sub_spe) |>
    add_fill("CellTypes") +
    scale_fill_manual(values = CellType_cols)
  ggsave(p,filename = here("plots","05_clustering","Banksy","Post_Annotation","CellType_v2",paste0(sample,".png")),
         height= 18, width = 20)
}


#Remake the gene expression heatmap. 
#Make a complexheatmap of genes
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
library(ComplexHeatmap)
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$CellTypes)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SNAP25","GAD1","PPP1R1B", #Broad neurons
                 "DRD1","RXFP1","TSHZ1",#D1/D1 islands
                 "SEMA5B","TRHDE","GABRQ","VWC2L","CPNE4", #D1_Island_A
                 "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP",#D1_IslandB
                 "RELN","TAC1","PDYN", #D1_MSN
                 "DRD2","ADORA2A","PENK","GPR6", #D2_MSN
                 "SLC18A3","SLC5A7","CHAT", #CHAT
                 "PVALB","GFRA2","KIT", #Inh_PVALB
                 "NPY","CORT","SST","CHODL",#Inh_SST 
                 "SLC17A7","TBR1", #EXCITATORY
                 "GJA1","AQP4","GFAP","TNC","WIF1", #AStrocyte
                 "CFAP157","ANXA1", #Ependymal
                 "CLDN5","NR2F2",
                 "OPALIN","MOBP","MOG","FGFR2","PROX1", #Oligo
                 "PDGFRA","VCAN",#OPC
                 "C3","P2RY13","MS4A6A" ) #Microglia

#marker labels
marker_labels <- c(rep("Neuron",3),
                   rep("D1_Island_A",8),
                   rep("D1_Island_B",6),
                   rep("D1_MSN",3),
                   rep("D2_MSN",4),
                   rep("CHAT",3),
                   rep("PVALB",3),
                   rep("SST",4),
                   rep("Excitatory",2),
                   rep("Astrocyte",5),
                   rep("Ependymal",2),
                   rep("Fibro",2),
                   rep("WM",5),
                   rep("OPC",2),
                   rep("Microglia",3))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))



dat <- assay(spe,"nucleus_normcounts")
dim(dat)

dat <- dat[markers_all,]
dim(dat)

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

min(hm_mat)
#[1] -1.363948

max(hm_mat)
#[1] 4.247616

#Reorder the matrix
hm_mat2 <- hm_mat[c("D1_Island_A","D1_Island_B","DRD1_MSN","DRD2_MSN",
                    "CHAT","Inh_PVALB","Inh_SST","Excitatory",
                    "Astro_A","Astro_B","Ependymal","Fibroblast_A","Fibroblast_B",
                    "WM","MSN_Oligo","Astrocyte_Oligo","Microglia_Oligo",
                    "OPC","Microglia_A","Microglia_B"),]

col_fun <- circlize::colorRamp2(c(-1.5,0,4.5),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat2,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
                              #column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              #bottom_annotation = col_ha, 
                              # right_annotation = row_ha,
                              column_split = marker_labels,
                              # row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = grid::gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)

pdf(file = here("plots","05_clustering","Banksy","Expression_heatmap_v2_celltypes.pdf"),height = 12, width = 18)
draw(hm)
dev.off()


#Save the object as celltype_v2
saveRDS(object = spe,file = here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))

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
#   [1] ComplexHeatmap_2.24.0       here_1.0.1                 
# [3] escheR_1.8.0                ggplot2_3.5.2              
# [5] SpatialExperiment_1.18.1    SingleCellExperiment_1.30.1
# [7] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [11] IRanges_2.42.0              S4Vectors_0.46.0           
# [13] BiocGenerics_0.54.0         generics_0.1.4             
# [15] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        dplyr_1.1.4             farver_2.1.2           
# [4] bluster_1.18.0          rsvd_1.0.5              digest_0.6.37          
# [7] lifecycle_1.0.4         cluster_2.1.8.1         Cairo_1.6-2            
# [10] statmod_1.5.0           magrittr_2.0.3          compiler_4.5.0         
# [13] rlang_1.1.6             tools_4.5.0             igraph_2.1.4           
# [16] S4Arrays_1.8.0          labeling_0.4.3          dqrng_0.4.1            
# [19] DelayedArray_0.34.1     RColorBrewer_1.1-3      abind_1.4-8            
# [22] BiocParallel_1.42.0     withr_3.0.2             beachmat_2.24.0        
# [25] colorspace_2.1-1        edgeR_4.6.2             scales_1.4.0           
# [28] iterators_1.0.14        dichromat_2.0-0.1       cli_3.6.5              
# [31] crayon_1.5.3            ragg_1.4.0              metapod_1.16.0         
# [34] httr_1.4.7              rjson_0.2.23            scuttle_1.18.0         
# [37] parallel_4.5.0          XVector_0.48.0          vctrs_0.6.5            
# [40] Matrix_1.7-3            jsonlite_2.0.0          BiocSingular_1.24.0    
# [43] GetoptLong_1.0.5        BiocNeighbors_2.2.0     irlba_2.3.5.1          
# [46] clue_0.3-66             systemfonts_1.2.3       magick_2.8.6           
# [49] locfit_1.5-9.12         foreach_1.5.2           limma_3.64.3           
# [52] glue_1.8.0              codetools_0.2-20        gtable_0.3.6           
# [55] shape_1.4.6.1           UCSC.utils_1.4.0        ScaledMatrix_1.16.0    
# [58] tibble_3.2.1            pillar_1.10.2           GenomeInfoDbData_1.2.14
# [61] circlize_0.4.16         R6_2.6.1                textshaping_1.0.1      
# [64] doParallel_1.0.17       rprojroot_2.0.4         lattice_0.22-7         
# [67] png_0.1-8               scran_1.36.0            Rcpp_1.0.14            
# [70] SparseArray_1.8.0       pkgconfig_2.0.3         GlobalOptions_0.1.2    
