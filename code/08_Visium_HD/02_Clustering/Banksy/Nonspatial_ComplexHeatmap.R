#Investigate 0.2 Lambda clusters
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(escheR)
library(here)


################################################################################
#   Read in object, says 0 but was run with lambda 0.2
################################################################################
# x will only have 2000 genes because took top 2000 HVGs for 
x <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_banksy_0_cell_level"
))

x

#Read in the filtered spe object
spe_filtered_path <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_path)

#Make sure cells are all in same order
identical(colnames(spe),colnames(x))

#Add cluster to object
spe$clust_M0_lam0.2_k50_res0.5 <- x$clust_M0_lam0.2_k50_res0.5
spe$clust_M0_lam0.2_k50_res0.8 <- x$clust_M0_lam0.2_k50_res0.8
spe$clust_M0_lam0.2_k50_res1   <- x$clust_M0_lam0.2_k50_res1

#  Work with resolution of 1.0 
#Make a complexheatmap of genes
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$clust_M0_lam0.2_k50_res1)

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


dat <- assay(spe, "logcounts")[markers_all, ]

mean_mat <- do.call(cbind, lapply(cell_idx, function(i) {
  Matrix::rowMeans(dat[, i, drop = FALSE])
}))

hm_mat <- scale(t(mean_mat), center = TRUE, scale = TRUE)

min(hm_mat)

max(hm_mat)

col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              column_split = marker_labels,
                              row_title = NULL,
                              rect_gp = grid::gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)

pdf(file = here("plots","HD_Full_Analysis","Banksy",
                "Cell_Level","Nonspatial","Expression",
                "Nonspatial_ComplexHeatmap.pdf"),height = 12, width = 18)
draw(hm)
dev.off()

################################################################################
#   Save clustered BANKSY object
################################################################################

message(Sys.time(), " | Saving BANKSY object")

saveHDF5SummarizedExperiment(spe, here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_banksy_0.2_cell_level"
))


sessionInfo()

