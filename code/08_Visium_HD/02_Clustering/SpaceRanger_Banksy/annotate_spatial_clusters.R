# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#  module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(ComplexHeatmap)
library(HDF5Array)
library(escheR)
library(here)

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

#Make a complexheatmap of genes
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

#Factorize the Banksy_CellType
set.seed(123)

# Create a grouping variable combining donor and cluster
group <- paste0(sfe$sample_id, "_", sfe$Spatial_Domain)

# Sample 25% from each group
keep_idx <- unlist(lapply(splitit(group), function(i) {
  n <- max(1, round(length(i) * 0.25))  # at least 1 cell per group
  sample(i, n)
}))

# Subset the sfe
sfe_sub <- sfe[, keep_idx]

cell_idx <- splitit(sfe_sub$Spatial_Domain)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SNAP25","GAD1","PPP1R1B","BCL11B","ISL1",#Broad neurons
                 "DRD1","RELN","TAC1","PDYN", "CALB1",#D1_MSN
                 "DRD2","ADORA2A","PENK","GPR6", #D2_MSN
                 "NPY","CORT","SST","CHODL",#Inh_SST 
                 "RXFP1","TSHZ1","OPRM1","FOXP2",#D1/D1 islands
                 "SEMA5B","TRHDE","GABRQ","VWC2L", #D1_Island_A
                 "DRD3","CPNE4","SEMA3E","PROK2","NPY1R","KCNH5","VIP",#D1_IslandB
                 "SLC18A3","SLC5A7","CHAT", #CHAT
                 "PVALB","GFRA2","KIT", #Inh_PVALB
                 "SLC17A7","TBR1", #EXCITATORY
                 "AVP","OXT","SIM1", 
                 "CLDN5","NR2F2",
                 "GJA1","AQP4","GFAP","TNC", #AStrocyte
                 "CFAP157","ANXA1", #Ependymal
                 "C3","P2RY13","MS4A6A", #Microglia
                 "OPALIN","MOBP","MOG") #WM

#marker labels
marker_labels <- c(rep("Neuron",2),
                   rep("MSN",3),
                   rep("D1_MSN",5),
                   rep("D2_MSN",4),
                   rep("Inhibitory",4),
                   rep("D1_Island_A",8),
                   rep("D1_Island_B",7),
                   rep("CHAT",3),
                   rep("PVALB",3),
                   rep("Excitatory",2),
                   rep("Hypo",3),
                   rep("Endothelial",2),
                   rep("Astrocyte",4),
                   rep("Ependymal",2),
                   rep("Microglia",3),
                   rep("WM",3))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE)


dat <- assay(sfe_sub,"logcounts")
dim(dat)

dat <- dat[markers_all,]
dim(dat)

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[c("MSN_Inh_A","MSN_Inh_B","MSN_Inh_C",
                   "D1_Island_A","D1_Island_B",
                   "CHAT_PVALB","Excitatory","Hypo",
                   "Endothelial","Endo_Micro_Astro","Ependymal","WM"),]

col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))

# Get the cluster labels from the heatmap rows
cluster_labels <- rownames(hm_mat)


#Load colors for the clusters 
spatial_colors <- readRDS(here("processed-data", "HD_Full_Analysis", "Cluster_colors",
                               "clust_M1_lam0.8_k50_res0.4_cell_level_spatial_colors_sr.Rds"))

names(spatial_colors) <- anno_df[match(names(spatial_colors),anno_df$spatial_0.4),"Annotation"]

# Build row annotation with your spatial colors
row_ha <- ComplexHeatmap::rowAnnotation(
  Spatial_Domain = cluster_labels,
  col = list(Spatial_Domain = spatial_colors),
  show_legend = TRUE
)

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              row_title = NULL,
                              bottom_annotation = col_ha,
                              column_title = NULL,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              rect_gp = grid::gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)

pdf(file = here("plots","HD_Full_Analysis",
                "Banksy_sr","spatial_0.4_Annotated_Expression_heatmap.pdf"),height = 12, width = 18)
draw(hm)
dev.off()


saveHDF5SummarizedExperiment(x = sfe,
                             dir = here("processed-data", 
                                        "HD_Full_Analysis", 
                                        "sfe_spatial_annotated"),
                             as.sparse = TRUE,
                             verbose = TRUE)
sessionInfo()
