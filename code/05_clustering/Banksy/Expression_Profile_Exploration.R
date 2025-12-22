#cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
library(ComplexHeatmap)
library(escheR)
library(here)

#Load in the annotated object --> this contains annotations from Banksy clustering. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

dim(spe)


spe


#genes <- c("DRD1","RXFP1","TSHZ1","OPRM1","FOXP2",#D1/D1 islands
#           "SEMA5B","TRHDE","GABRQ","VWC2L",#D1_Island_A
#           "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP",#D1_IslandB
#           "RELN","TAC1","PDYN")
#spe_sub <- spe[,spe$Donor == "Br6660"]
#logcounts(spe_sub) <- assay(spe_sub,"nucleus_normcounts")
#for(sample in unique(spe_sub$Sample)){
#  print(sample)
#  spe_sub2 <- spe_sub[,spe_sub$Sample == sample]
#  for(gene in genes){
#    print(gene)
#    spe_sub2$Gene_Expression <- logcounts(spe_sub2)[gene,]
#    p <- make_escheR(spe_sub2) |>
#      add_fill("Gene_Expression") +
#      scale_fill_gradientn(colors = c("lightgrey","orange","red"))
#    ggsave(filename = here("plots","05_clustering","GeneExpression",
#                           "D1_Island_Markers",paste0(gene,"_",sample,".png")),
#           height = 22, width = 22)
#  }
#}

#Load the annotated banksy cell types
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

head(Banksy_celltypes)

stopifnot(identical(Banksy_celltypes$cell_id,colnames(spe)))

spe$Banksy_celltypes <- Banksy_celltypes$CellType

#Co-expression matrices. Are D1_Island_B in fact D1_Islands? 
# ---- user inputs ----
genes <- c("DRD1","RXFP1","TSHZ1","OPRM1","FOXP2",#D1/D1 islands
           "SEMA5B","TRHDE","GABRQ","VWC2L",#D1_Island_A
           "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP",#D1_IslandB
           "RELN","TAC1","PDYN")  
cluster_col <- "Banksy_celltypes"                                 
assay_name <- "nucleus_normcounts"                          
expr_threshold <- 0                                       

# ---- checks ----
cluster_col %in% colnames(colData(spe))
#[1] TRUE
genes_present <- genes[genes %in% rownames(spe)]

genes_present


# ---- helper: percent co-expression matrix within ONE cluster ----
pct_coexp_one_cluster <- function(Xg) {
  # Xg: dgCMatrix genes x cells (subset to one cluster)
  # Convert to logical expressed matrix
  E <- Xg > expr_threshold
  
  # If E is logical sparse, crossprod is fast and memory-friendly
  both <- as.matrix(E %*% t(E))  # Matrix multiplication to create a logical matrix consisting of TRUE/FALSE depending on whether value is >0.
  #Each non-diagnoal value in both is the number of cells in the cluster that have non-zero expression values for two genes. 
  n_cells <- ncol(E)
  
  # percent of cells in cluster that co-express two genes: %(genei & genej)
  pct  <- both / n_cells * 100 
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
  #Pull matrix. Diagonal is percent of cells in the cluster that express the gene/ 
  mat <- pct_list[[i]]
  #The matrix is symmetric so remove the upper half. 
  mat[upper.tri(mat)] <- NA
  #diag(mat) <- NA
  pdf(here("plots","05_clustering","GeneExpression","D1_Island_Markers","CoExpression",paste0(i,".pdf")), width = 8, height = 8)
  pheatmap::pheatmap(
    mat,
    main = paste0(i, "\n% of Cells in Cluster Co-Expressing Two Genes"),
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    na_col = "grey95",
    border_color = "black"
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



dat <- assay(spe,"logcounts")
dim(dat)

dat <- dat[markers_all,]
dim(dat)

dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)



min(hm_mat)

max(hm_mat)

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

pdf(file = here("plots","05_clustering","Banksy","Expression_heatmap.pdf"),height = 12, width = 18)
draw(hm)
dev.off()


sessionInfo()
