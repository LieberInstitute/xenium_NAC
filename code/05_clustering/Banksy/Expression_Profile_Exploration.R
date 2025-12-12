#cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
library(ComplexHeatmap)
library(escheR)
library(here)

#Load in the annotated object --> this contains annotations from Banksy clustering. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

dim(spe)
#[1]     366 4849373

spe
# class: SpatialExperiment 
# dim: 366 4849373 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4849373): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

genes <- c("DRD1","RXFP1","TSHZ1","OPRM1","FOXP2",#D1/D1 islands
           "SEMA5B","TRHDE","GABRQ", #D1_Island_A
           "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP",#D1_IslandB
           "RELN","TAC1","PDYN")
spe_sub <- spe[,spe$Donor == "Br6660"]
logcounts(spe_sub) <- assay(spe_sub,"nucleus_normcounts")
for(sample in unique(spe_sub$Sample)){
  print(sample)
  spe_sub2 <- spe_sub[,spe_sub$Sample == sample]
  for(gene in genes){
    print(gene)
    spe_sub2$Gene_Expression <- logcounts(spe_sub2)[gene,]
    p <- make_escheR(spe_sub2) |>
      add_fill("Gene_Expression") +
      scale_fill_gradientn(colors = c("lightgrey","orange","red"))
    ggsave(filename = here("plots","05_clustering","GeneExpression",
                           "D1_Island_Markers",paste0(gene,"_",sample,".png")),
           height = 22, width = 22)
  }
}

#Load the annotated banksy cell types
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

head(Banksy_celltypes)
# cell_id    CellType
# 1 Br6660_NAc1_580_1  Excitatory
# 2 Br6660_NAc1_580_2  Excitatory
# 3 Br6660_NAc1_580_3  Excitatory
# 4 Br6660_NAc1_580_4 Microglia_A
# 5 Br6660_NAc1_580_5  Excitatory
# 6 Br6660_NAc1_580_6        WM_A

identical(Banksy_celltypes$cell_id,colnames(spe))
#[1] TRUE

spe$Banksy_celltypes <- Banksy_celltypes$CellType

#Co-expression matrices. Are D1_Island_B in fact D1_Islands? 
# ---- user inputs ----
genes <- c("SNAP25","PPP1R1B","BCL11B",
           "DRD1","RXFP1","TSHZ1","FOXP2","GABRQ",
           "SEMA5B","TRHDE","SEMA3E","PROK2",
           "NPY1R","RXFP2","KCNH5","VIP",
           "RELN","TAC1","PDYN")  
cluster_col <- "Banksy_celltypes"                                 
assay_name <- "nucleus_normcounts"                          
expr_threshold <- 0                                       

# ---- checks ----
cluster_col %in% colnames(colData(spe))
#[1] TRUE
genes_present <- genes[genes %in% rownames(spe)]

genes_present
# [1] "SNAP25"  "PPP1R1B" "BCL11B"  "DRD1"    "RXFP1"   "TSHZ1"   "FOXP2"  
# [8] "GABRQ"   "SEMA5B"  "TRHDE"   "SEMA3E"  "PROK2"   "NPY1R"   "RXFP2"  
# [15] "KCNH5"   "VIP"     "RELN"    "TAC1"    "PDYN" 

# ---- helper: percent co-expression matrix within ONE cluster ----
pct_coexp_one_cluster <- function(Xg) {
  # Xg: dgCMatrix genes x cells (subset to one cluster)
  # Convert to logical expressed matrix
  E <- Xg > expr_threshold
  
  # If E is logical sparse, crossprod is fast and memory-friendly
  both <- as.matrix(E %*% t(E))  # Matrix multiplication to create a logical matrix consisting of TRUE/FALSE depending on whether value is >0.
  #Each non-diagnoal value in both is the number of cells in the cluster that have non-zero expression values for two genes. 
  #However, the diagonal is the number of cells in the cluster that have non-zero values for that one specific gene. 
  cell_totals <- diag(both) #Save for calculation. 
  
  # percent of cells expressing column gene that also express row gene:
  pct <- sweep(both, 2, cell_totals, FUN = "/") * 100
  pct
}

# ---- compute per-cluster heatmaps ----
clusters <- as.character(colData(spe)[[cluster_col]]) #Pull cluster identities
X <- assay(spe, assay_name)
X <- X[genes_present, , drop = FALSE]

# split cell indices by cluster
idx_list <- split(seq_len(ncol(spe)), clusters)

# build a list of percent matrices (one per cluster)
pct_list <- lapply(idx_list, function(idx) {
  Xg <- X[, idx, drop = FALSE]
  pct_coexp_one_cluster(Xg)
})


for(i in names(pct_list)){
  print(i)
  mat <- pct_list[[i]]
  diag(mat) <- NA
  pdf(here("plots","05_clustering","GeneExpression","D1_Island_Markers","CoExpression",paste0(i,".pdf")), width = 8, height = 8)
  pheatmap::pheatmap(
    mat,
    main = paste0(i, "\n% of Cells Expressing Column Gene that Also Express Row Genes"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "grey95",
    border_color = NA
  )
  dev.off()
}



#Factorize the Banksy_CellType
spe$Banksy_celltypes <- factor(x = spe$Banksy_celltypes,
                              levels = c("D1_Island_A","D1_Island_B","DRD1_MSN","DRD2_MSN",
                                         "CHAT","Inh_PVALB","Inh_SST","Excitatory","MSN_Oligo",
                                         "Astro_A","Astro_B","Ependymal",
                                         "Fibroblast_A","Fibroblast_B",
                                         "WM_A","WM_B","WM_C","WM_D","WM_E","OPC",
                                         "Microglia_A","Microglia_B"))


logcounts(spe) <- assay(spe,"nucleus_normcounts")

#Make a complexheatmap of genes
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$Banksy_celltypes)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SNAP25","GAD1","PPP1R1B", #Broad neurons
                 "DRD1","RXFP1","TSHZ1",#D1/D1 islands
                 "SEMA5B","TRHDE","GABRQ", #D1_Island_A
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
                   rep("D1_Island_A",6),
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



dat <- assay(spe,"logcounts")
dim(dat)
#[1]     366 4849373

dat <- dat[markers_all,]
dim(dat)
#[1]      53 4849373

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)



min(hm_mat)
#[1] -1.255829

max(hm_mat)
#[1] 4.47619

col_fun <- circlize::colorRamp2(c(-1.5,0,4.5),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
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

pdf(file = here("plots","05_clustering","Banksy","Expression_heatmap.pdf"),height = 12, width = 12)
draw(hm)
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
#   [1] here_1.0.1                  escheR_1.8.0               
# [3] ggplot2_3.5.2               ComplexHeatmap_2.24.0      
# [5] SpatialExperiment_1.18.1    SingleCellExperiment_1.30.1
# [7] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [11] IRanges_2.42.0              S4Vectors_0.46.0           
# [13] BiocGenerics_0.54.0         generics_0.1.4             
# [15] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6            circlize_0.4.16         shape_1.4.6.1          
# [4] rjson_0.2.23            GlobalOptions_0.1.2     lattice_0.22-7         
# [7] Cairo_1.6-2             vctrs_0.6.5             tools_4.5.0            
# [10] parallel_4.5.0          tibble_3.2.1            cluster_2.1.8.1        
# [13] pkgconfig_2.0.3         pheatmap_1.0.12         Matrix_1.7-3           
# [16] RColorBrewer_1.1-3      lifecycle_1.0.4         GenomeInfoDbData_1.2.14
# [19] compiler_4.5.0          farver_2.1.2            textshaping_1.0.1      
# [22] codetools_0.2-20        clue_0.3-66             pillar_1.10.2          
# [25] crayon_1.5.3            DelayedArray_0.34.1     magick_2.8.6           
# [28] iterators_1.0.14        abind_1.4-8             foreach_1.5.2          
# [31] tidyselect_1.2.1        digest_0.6.37           dplyr_1.1.4            
# [34] labeling_0.4.3          rprojroot_2.0.4         colorspace_2.1-1       
# [37] cli_3.6.5               SparseArray_1.8.0       magrittr_2.0.3         
# [40] S4Arrays_1.8.0          dichromat_2.0-0.1       withr_3.0.2            
# [43] scales_1.4.0            UCSC.utils_1.4.0        XVector_0.48.0         
# [46] httr_1.4.7              ragg_1.4.0              png_0.1-8              
# [49] GetoptLong_1.0.5        doParallel_1.0.17       viridisLite_0.4.2      
# [52] rlang_1.1.6             Rcpp_1.0.14             glue_1.8.0             
# [55] jsonlite_2.0.0          R6_2.6.1                systemfonts_1.2.3      
# 
# 
# 
# 
# 
