#Goal: Generate ComplexHeatmap for figure 1 of the NAc A-P manuscript. Also save legend of escheR for figure 1. '
# /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# module load conda_R/4.5

library(SpatialExperiment)
library(ComplexHeatmap)
library(cowplot)
library(escheR)
library(here)

#Read in the spe object
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))

spe

#colors
cluster_cols <- readRDS(here("processed-data","CellType_Cols_v2.Rds"))

###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
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
                   rep("DRD1_MSN",3),
                   rep("DRD2_MSN",4),
                   rep("CHAT",3),
                   rep("Inh_PVALB",3),
                   rep("Inh_SST",4),
                   rep("Excitatory",2),
                   rep("Astro_A",4),
                   rep("Astro_B",1),
                   rep("Ependymal",2),
                   "Fibroblast_A",
                   "Fibroblast_B",
                   rep("WM",5),
                   rep("OPC",2),
                   rep("Microglia_A",2),
                   "Microglia_B")

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))


colors_markers <- list(marker = cluster_cols[match(unique(marker_labels),names(cluster_cols))])

# After creating colors_markers, drop any NAs and manually add Neuron
colors_markers$marker <- colors_markers$marker[!is.na(colors_markers$marker)]
colors_markers$marker["Neuron"] <- "gray60"
col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(D1_Island_A = "D1_Island_A",
                     D1_Island_B = "D1_Island_B",
                     DRD1_MSN = "DRD1_MSN",
                     DRD2_MSN = "DRD2_MSN",
                     CHAT = "CHAT",
                     Inh_PVALB = "Inh_PVALB",
                     Inh_SST = "Inh_SST",
                     Excitatory = "Excitatory",
                     Astro_A = "Astro_A",
                     Astro_B = "Astro_B",
                     Ependymal = "Ependymal",
                     Fibroblast_A = "Fibroblast_A",
                     Fibroblast_B = "Fibroblast_B",
                     WM = "WM",
                     MSN_Oligo = "MSN_Oligo",
                     Astrocyte_Oligo = "Astrocyte_Oligo",
                     Microglia_Oligo = "Microglia_Oligo",
                     OPC = "OPC",
                     Microglia_A = "Microglia_A",
                     Microglia_B = "Microglia_B")

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



dat <- assay(spe,"nucleus_normcounts")
dim(dat)


dat <- dat[markers_all,]
dim(dat)


dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]

max(hm_mat)

min(hm_mat)


col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))



hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "Gene expression across clusters",
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


pdf(file = here("plots","Figure1","Cluster_Heatmap.pdf"),
    width = 16,
    height = 12)
draw(hm)
dev.off()

sessionInfo()
