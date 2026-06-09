# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# moudle load ocnda

library(SpatialExperiment)
library(DeconvoBuddies)
library(ComplexHeatmap)
library(circlize)
library(HDF5Array)
library(ggplot2)
library(here)

###### Prep spatial data
#Read in the filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

spe

#Add spatial clusters and non-spatial clusters
spatial_clust <- read.csv(here("processed-data", "HD_Full_Analysis",
                               "spatial_banksy_clusters_res_0.4.csv"))

#Add spatial data to spe 
stopifnot(identical(spatial_clust$V2,colnames(spe)))

spe$spatial_0.4 <- spatial_clust$V1

#Load colors for the clusters 
spatial_colors <- readRDS(here("processed-data", "HD_Full_Analysis", "Cluster_colors",
                               "clust_M1_lam0.8_k50_res0.4_cell_level_spatial_colors.Rds"))

#Make a complexheatmap of genes
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

#Factorize the Banksy_CellType
set.seed(123)

# Create a grouping variable combining donor and cluster
group <- paste0(spe$sample_id, "_", spe$spatial_0.4)

# Sample 25% from each group
keep_idx <- unlist(lapply(splitit(group), function(i) {
  n <- max(1, round(length(i) * 0.25))  # at least 1 cell per group
  sample(i, n)
}))

# Subset the spe
spe_sub <- spe[, keep_idx]

cell_idx <- splitit(spe_sub$spatial_0.4)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SNAP25","GAD1","PPP1R1B", #Broad neurons
                 "DRD1","RXFP1","TSHZ1","OPRM1",#D1/D1 islands
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

dat <- assay(spe_sub,"logcounts")
dim(dat)

dat <- dat[markers_all,]
dim(dat)

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))

# Get the cluster labels from the heatmap rows
cluster_labels <- rownames(hm_mat)

# Build row annotation with your spatial colors
row_ha <- ComplexHeatmap::rowAnnotation(
  Cluster = cluster_labels,
  col = list(Cluster = spatial_colors),
  show_legend = TRUE
)

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              row_title = NULL,
                              rect_gp = grid::gpar(col = "gray50", lwd = 0.5),
                              col = col_fun,
                              left_annotation = row_ha)

pdf(file = here("plots","HD_Full_Analysis",
                "Banksy","spatial_0.4_Expression_heatmap.pdf"),height = 12, width = 18)
draw(hm)
dev.off()


################################################################################
#  Find Markers of each cluster to help with annotation
################################################################################
message("Calculate mean ratio - ",Sys.time())
#Run mean ratio
gmr <- get_mean_ratio(spe,
                       cellType_col = "spatial_0.4",
                       assay_name = "logcounts")

save(gmr,
     file = here("processed-data","HD_Full_Analysis","Banksy","spatial_0.4_getmeanratio.rda"))

###############
message("Starting DEG testing",Sys.time())
################################################
message("Pairwise DEG testing",Sys.time())
#Pairwise DEG testing
mod <- with(colData(spe), model.matrix(~ sample_id))
mod <- mod[ , -1, drop=F] # intercept otherwise automatically dropped by `findMarkers()`

# Run pairwise t-tests
markers_pairwise <- findMarkers(spe, 
                                groups=spe$spatial_0.4,
                                assay.type="logcounts", 
                                design=mod, 
                                test="t",
                                direction="up", 
                                pval.type="all", 
                                full.stats=T)

#How many DEGs for each cluster? 
sapply(markers_pairwise, function(x){table(x$FDR<0.05)})

save(markers_pairwise,
     file = here("processed-data","HD_Full_Analysis","Banksy","spatial_0.4_pairwise.rda"))


################################################
message("1vALL DEG testing",Sys.time())
markers_1vALL_enrich <- findMarkers_1vAll(spe, 
                                          assay_name = "logcounts", 
                                          cellType_col = "spatial_0.4", 
                                          mod = "~sample_id")


save(markers_1vALL_enrich,file = here("processed-data","HD_Full_Analysis","Banksy","spatial_0.4_1vALL.rda"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
