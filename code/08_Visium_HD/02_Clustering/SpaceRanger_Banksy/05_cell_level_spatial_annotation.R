# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(HDF5Array)
library(escheR)
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

#Annotate the spatial domains
anno_df <- read.csv(here("processed-data","HD_Full_Analysis","spatial_0.4_Annotations.csv"))
sfe$spatial_domain <- anno_df[match(sfe$spatial_0.4,anno_df$spatial_0.4),"Annotation"]

#Load the spatial colors
cluster_cols <-  readRDS(here("processed-data", "HD_Full_Analysis", "Cluster_colors",
                      "clust_M1_lam0.8_k50_res0.4_cell_level_spatial_colors_sr.Rds"))

names(cluster_cols) <- anno_df[match(names(cluster_cols),anno_df$spatial_0.4),"Annotation"]

# Plot per sample
for (sample in unique(sfe$sample_id)) {
  message(Sys.time(), " | Plotting: ", sample)
  sfe_sub <- sfe[, sfe$sample_id == sample]
  p <- make_escheR(sfe_sub) |>
    add_fill("spatial_domain") +
    scale_fill_manual(values = cluster_cols)
  ggsave(
    here("plots", "HD_Full_Analysis", "Banksy_sr",
         "Cell_Level", "Spatial_Annotated",
         paste0(sample,".png")),
    p, width = 20, height = 14, dpi = 200
  )
}


#Complex Heatmap 
###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(sfe$spatial_domain)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SNAP25","GAD1","GAD2","SLC32A1", # GABA neurons
                 "PPP1R1B","FOXP2","BCL11B", #Broad neurons
                 "NPY","CORT","SST","CHODL", #Inhibitory subset
                 "DRD1","RELN","TAC1","PDYN", #D1_MSN
                 "DRD2","ADORA2A","PENK","GPR6", #D2_MSN
                 "RXFP1","TSHZ1",#D1/D1 islands
                 "SEMA5B","TRHDE","GABRQ","VWC2L","CPNE4", #D1_Island_A
                 "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP",#D1_IslandB
                 "SLC18A3","SLC5A7","CHAT","PVALB","GFRA2","KIT", #CHAT
                 "SLC17A7","TBR1", #EXCITATORY
                 "CLDN5","NR2F2", #Endo
                 "C3","P2RY13","MS4A6A", #Microglia
                 "GJA1","AQP4","GFAP", #AStrocyte
                 "CFAP157","CAPS", #Ependymal
                 "OPALIN","MOBP","MOG", #Oligo
                 "PDGFRA","VCAN") #OPC)

#marker labels
marker_labels <- c(rep("GABA",4),
                   rep("MSN",3),
                   rep("SST",4),
                   rep("DRD1_MSN",4),
                   rep("DRD2_MSN",4),
                   rep("D1_Islands",13),
                   rep("CHAT_PVALB",6),
                   rep("Excitatory",2),
                   rep("Endo",2),
                   rep("Microglia",3),
                   rep("Astrocyte",3),
                   rep("WM",3),
                   rep("OPC",2),
                   rep("Ependymal",2))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

# After creating colors_markers, drop any NAs and manually add Neuron
col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(GABA_A = "GABA_A",
                     GABA_B = "GABA_B",
                     GABA_C = "GABA_C",
                     D1_Island_A = "D1_Island_A",
                     D1_Island_B = "D1_Island_B",
                     CHAT_PVALB = "CHAT_PVALB",
                     Excitatory = "Excitatory",
                     Endothelial = "Endothelial",
                     Endo_Micro_Astro = "Endo_Micro_Astro",
                     Endo_Astro = "Endo_Astro",
                     Ependymal = "Ependymal",
                     WM = "WM")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))

cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = cluster_cols[names(cluster_pops)])

row_ha <- rowAnnotation(population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = FALSE,
                        col = pop_markers)



dat <- assay(sfe,"logcounts")[markers_all,]
dim(dat)

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]

col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))



hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "Gene expression across spatial domains",
                              column_title_gp = gpar(fontface = "bold"),
                              cluster_rows = FALSE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              column_split = marker_labels,
                              row_split = cluster_pops_rev,
                              row_title = NULL,
                              rect_gp = gpar(col = "gray50", lwd = 0.5),
                              col = col_fun)


pdf(file = here("plots","HD_Full_Analysis","Banksy_sr",
                "Cell_Level", "Spatial_Annotated","Annotated_Spatial_Domains_Heatmap.pdf"),
    width = 16,
    height = 12)
draw(hm)
dev.off()


saveHDF5SummarizedExperiment(x = sfe,dir = here("processed-data", "HD_Full_Analysis", "sfe_annotated_labels"))

sessionInfo()
