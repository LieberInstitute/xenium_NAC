# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(scuttle)
library(escheR)
library(here)


#Load sfe with nmf info
sfe <-  loadHDF5SummarizedExperiment(here("processed-data", "HD_Full_Analysis", "sfe_with_NMF"))

sfe


# Add RCTD results to the sfe object
rctd_coldata <- readRDS(here("processed-data", "HD_Full_Analysis","sfe_cell_RCTD_doublet_colData.Rds"))

stopifnot(identical(rownames(rctd_coldata),colnames(sfe)))
stopifnot(identical(rownames(rctd_coldata),rownames(colData(sfe))))

rctd_cols <- colnames(rctd_coldata)[grep("rctd",colnames(rctd_coldata))]

colData(sfe) <- cbind(colData(sfe),rctd_coldata[,rctd_cols])

# Add singleR results to the Sfe
files <- list.files(here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger","SingleR"),full.names = TRUE)

singler_res <- lapply(X = files,FUN = readRDS)
singler_res <- do.call(what = rbind,singler_res)

scores_mat <- singler_res$scores
labels     <- singler_res$labels   

label_scores <- scores_mat[cbind(seq_len(nrow(scores_mat)),
                                 match(labels, colnames(scores_mat)))]

singler_res$label_score <- label_scores

singler_res <- as.data.frame(singler_res)

#Force same order
singler_res <- singler_res[colnames(sfe),]
stopifnot(identical(rownames(singler_res),colnames(sfe)))

#Add single r to sfe
colData(sfe) <- cbind(colData(sfe),singler_res)


sfe$snRNA_label <- ifelse(sfe$rctd_spot_class == "singlet",
                          sfe$rctd_first_type,
                          sfe$pruned.labels)

#11 cells don't have labels. Move on without them. 
sfe <- sfe[,!is.na(sfe$snRNA_label)]

sfe

#Load colors
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE) 

plot_dir <- here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger","snRNA_labels")

#Plot all cell types on each sample
for(sample in unique(sfe$sample_id)){
  print(sample)
  sfe_sub <- sfe[,sfe$sample_id == sample]
  p <- make_escheR(sfe_sub) |>
    add_fill("snRNA_label") +
    scale_fill_manual(values = cluster_cols)
  ggsave(plot = p,filename = file.path(plot_dir,paste0(sample,".png")),
         height = 10, width = 12,dpi = 200)
}

#Plot each cell type by itself on the sample 
for(ct in unique(sfe$pruned.labels)){
  print(ct)
  sfe$ct <- ifelse(sfe$snRNA_label == ct,
                   ct,
                   "Other")
  dir.create(path = file.path(plot_dir,ct))
  for(sample in unique(sfe$sample_id)){
    print(sample)
    sfe_sub <- sfe[,sfe$sample_id == sample]
    p <- make_escheR(sfe_sub) |>
      add_fill("ct") +
      scale_fill_manual(values = cluster_cols)
    ggsave(plot = p,filename = file.path(plot_dir,ct,paste0(sample,".png")),
           height = 12, width = 16,dpi = 200)
  }
}


###########Complext heatmap of basic markers
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

set.seed(123)

# Create a grouping variable combining donor and cluster
group <- paste0(sfe$sample_id, "_", sfe$snRNA_label)

# Sample 20% from each group
keep_idx <- unlist(lapply(splitit(group), function(i) {
  n <- max(1, round(length(i) * 0.20))  # at least 1 cell per group
  sample(i, n)
}))

# Subset the sfe
sfe_sub <- sfe[, keep_idx]

cell_idx <- splitit(sfe_sub$snRNA_label)


############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
markers_all <- c("SYT1","SNAP25","GAD1","PPP1R1B", #Broad neurons
                 "DRD1","RXFP1", #D1/D1 islands
                 "RELN","CNTNAP3B", #D1_A
                 "TRHDE","CPNE4",#D1_B
                 "RXFP2","SEMA3E", #D1_C
                 "VWC2L","CLSTN2", #D1_D
                 "DRD2","ADORA2A", #General D2
                 "PENK","PTPRM",#D2_A
                 "CLMP","GRIK3",#D2_B
                 "IL1RAPL2","PDGFD", #Inh_A
                 "VIP","CCK", #Inh_B
                 "GLP1R","TAC3", #Inh_C
                 "CHAT","SLC5A7", #Inh_D
                 "NPY","SST", #Inh_E
                 "KCNC2","ANK1", #Inh_F
                 "SLC17A7","TBR1", #Excitatory
                 "AQP4","GFAP", #Astrocyte
                 "CAPS","FOXJ1", #Ependymal
                 "MBP","MOBP", #Oligo
                 "PDGFRA","VCAN", #OPC
                 "C3","DOCK8", #Microglia
                 "DCN","CLDN5") #Endothelial  

#marker labels
marker_labels <- c(rep("Neuron",4),
                   rep("D1_MSN",10),
                   rep("D2_MSN",6),
                   rep("Inhibitory",12),
                   rep("Excitatory",2),
                   rep("Astrocyte",2),
                   rep("Ependymal",2),
                   rep("Oligo",2),
                   rep("OPC",2),
                   rep("Microglia",2),
                   rep("Endothelial",2))

marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c(Neuron = "black",
                                  D1_MSN = "#332288",
                                  D2_MSN = "#D81B60",
                                  Inhibitory = "#44AA99",
                                  Excitatory = as.character(cluster_cols["Excitatory"]),
                                  Astrocyte = "#DDCC77",
                                  Ependymal = as.character(cluster_cols["Ependymal"]),
                                  Oligo = as.character(cluster_cols["Oligo"]),
                                  OPC = as.character(cluster_cols["OPC"]),
                                  Microglia = as.character(cluster_cols["Microglia"]),
                                  Endothelial = as.character(cluster_cols["Endothelial"])))

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
# cluster labels
cluster_pops <- list(D1_MSN = c("DRD1_MSN_A","DRD1_MSN_B",
                                "DRD1_MSN_C","DRD1_MSN_D"),
                     D2_MSN = c("DRD2_MSN_A","DRD2_MSN_B"),
                     Inhibitory = c("Inh_A","Inh_B","Inh_C",
                                    "Inh_D","Inh_E","Inh_F"),
                     Excitatory = "Excitatory",
                     Astrocyte = c("Astrocyte_A","Astrocyte_B"),
                     Ependymal = "Ependymal",
                     Oligo = "Oligo",
                     OPC = "OPC",
                     Microglia = "Microglia",
                     Endothelial = "Endothelial")

# cluster labels order
# # cluster labels order
cluster_pops_order <- unname(unlist(cluster_pops))

# swap values and names of list
cluster_pops_rev <- rep(names(cluster_pops),
                        times = sapply(cluster_pops, length))
names(cluster_pops_rev) <- unname(unlist(cluster_pops))
#cluster_pops_rev <- cluster_pops_rev[as.character(sort(cluster_pops_order))]
cluster_pops_rev <- factor(cluster_pops_rev, levels = names(cluster_pops))

# second set of cluster labels
neuron_pops <- ifelse(cluster_pops_rev %in% c("Inhibitory","D1_MSN",
                                              "D2_MSN","Excitatory"),
                      "Neuron",
                      "Non-neuron")

neuron_pops <- factor(x = neuron_pops,levels = c("Neuron","Non-neuron"))

colors_neurons <- list(class = c(Neuron = "black",
                                 `Non-neuron` = "gray65"))

#n <- table(sce$CellType.Final)

#row annotation dataframe. 
# row annotation
pop_markers <- list(population = c(D1_MSN = "#332288",
                                   D2_MSN = "#D81B60",
                                   Inhibitory = "#44AA99",
                                   Excitatory = as.character(cluster_cols["Excitatory"]),
                                   Astrocyte = "#DDCC77",
                                   Ependymal = as.character(cluster_cols["Ependymal"]),
                                   Oligo = as.character(cluster_cols["Oligo"]),
                                   OPC = as.character(cluster_cols["OPC"]),
                                   Microglia = as.character(cluster_cols["Microglia"]),
                                   Endothelial = as.character(cluster_cols["Endothelial"])))

row_ha <- rowAnnotation(class = neuron_pops,
                        population = cluster_pops_rev,
                        show_annotation_name = FALSE,
                        show_legend = FALSE,
                        col = c(pop_markers,colors_neurons))


dat <- assay(sfe_sub,"logcounts")[markers_all,]
#rownames(dat) <- rowData(sce)$gene_name
dim(dat)


dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)

hm_mat <- hm_mat[names(cluster_pops_rev),]




col_fun <- circlize::colorRamp2(c(min(hm_mat),0,max(hm_mat)),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              column_title = "General cell class marker \ngene expression across clusters",
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


pdf(file = file.path(plot_dir,"VHD_snRNA_label_Heatmap.pdf"),
    width = 16,
    height = 12)
draw(hm)
dev.off()


#Save object (now contains nmf{i} and nmf{i}_scaled columns)
saveHDF5SummarizedExperiment(sfe, 
                             here("processed-data", "HD_Full_Analysis", "sfe_with_labels"),
                             replace = TRUE)


#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
