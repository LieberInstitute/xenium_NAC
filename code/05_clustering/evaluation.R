#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5

library(cluster)
library(here)

# DONOR <- args[1]
DONOR <- "Br6660"  # or "Br6436" or "All"
NUM_CLUST = 15


emb_file <- sprintf(
  here("processed-data", "05_Clustering", "STAligner_r0_%s_%d_embeddings.tsv"),
  DONOR, NUM_CLUST
)
cluster_file <- sprintf(
  here("processed-data", "05_Clustering", "STAligner_r0_%s_%d_clusters.csv"),
  DONOR, NUM_CLUST
)

emb <- read.table(emb_file,
  header      = FALSE,
  sep         = "\t",
  check.names = FALSE
)

clusters_df <- read.csv(
  cluster_file,
  header    = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

cluster_lab <- clusters_df$mclust
# sanity‑check dimensions
n_cells <- nrow(emb)
if (length(cluster_lab) != n_cells) {
  stop(sprintf(
    "Embedding has %d rows but cluster vector has length %d!",
    n_cells, length(cluster_lab)
  ))
}

set.seed(42)
sample_size <- min(50000, n_cells)
idx <- sample(n_cells, size = sample_size, replace = FALSE)
emb_sub     <- emb[idx, , drop = FALSE]
cluster_sub <- cluster_lab[idx]
cluster_int <- as.integer(as.factor(cluster_sub))
dist_mat <- dist(emb_sub, method = "euclidean")
sil <- silhouette(cluster_int, dist_mat)
sil_score <- mean(sil[ , "sil_width"])
cat(sprintf(
  "DONOR=%s    n=%d    silhouette=%.4f\n",
  DONOR, sample_size, sil_score
))

# DONOR=Br6660    n=50000    silhouette=0.1519
# DONOR=Br6436    n=50000    silhouette=0.1498
# DONOR=All    n=50000    silhouette=0.1394






library(scDotPlot)
library(SpatialExperiment)
library(ggplot2)
library(ggplotify)
library(ragg)

spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))
assay(spe, "logcounts") <- assay(spe, "nucleus_normcounts")
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]
spe_br6660 <- spe[, spe$Donor == "Br6660"]
colData(spe_br6660) <- S4Vectors::droplevels(colData(spe_br6660))
spe_br6660$STAligner_cluster <- clusters_df[colnames(spe_br6660), "mclust"]
spe_br6660$STAligner_cluster <- as.character(spe_br6660$STAligner_cluster)


markers_all <- c("SNAP25","PPP1R1B","BCL11B","GAD1",
                 "DRD1","PDYN","DRD2","PENK","ADORA2A",
                 "FOXP2","GABRQ","SEMA5B","CPNE4","VIP",
                 "CHAT","SLC5A7",
                 "SST","NPY","PNOC","CHODL",
                 "SLC17A7","TBR1",
                 "MOBP","ST18","OLIG1","OLIG2","OPALIN",
                 "SLC1A2","GFAP","GJA1","AQP4",
                 "CFAP157",
                 "C3","ARHGAP15",
                 "DCN","EBF1",
                 "TSHZ1", "CASZ1","RXFP1", "CORT", "GLP1R", "ONECUT2")

p <- scDotPlot(object = spe_br6660,
                features = rev(markers_all),
                group = "STAligner_cluster",
                groupAnno = "STAligner_cluster",
                clusterColumns = TRUE,
                clusterRows = TRUE,
                scale = TRUE)

outfile <- here(
  "plots","05_clustering","STAligner","Pre_Annotation",
   DONOR, paste0("k", NUM_CLUST), 
  sprintf("STAligner_%s_mclust%d_heatmap.png", DONOR, NUM_CLUST)
)

ragg::agg_png(outfile, width = 12, height = 10, units = "in", res = 300)
grid::grid.newpage()
grid::grid.draw(ggplotify::as.grob(p))   # or: print(ggplotify::as.ggplot(p))
dev.off()



# scDotPlot(object = spe,
#           features = rev(markers_all),
#           group = "Annotation",
#           groupAnno = "Annotation",
#           clusterColumns = TRUE,
#           clusterRows = TRUE,
#           scale = TRUE,
#           annoColors = list("Annotation" = cluster_cols))



#Generate a pseudobulked object
#Use AggregateAcrossCells to Pseudobulk the data and generate some boxplots. 
ids_df <- S4Vectors::DataFrame(
  Sample         = colData(spe)$Sample,
  Annotation = colData(spe)$Annotation
)

spe_pb <- scuttle::aggregateAcrossCells(spe, ids = ids_df, use.assay.type = "counts")

#Add the column names back. 
colnames(spe_pb) <- paste(colData(spe_pb)$Annotation,colData(spe_pb)$Sample,sep = ".")
rownames(colData(spe_pb)) <- paste(colData(spe_pb)$Annotation,colData(spe_pb)$Sample,sep = ".")
identical(colnames(spe_pb),rownames(colData(spe_pb)))
#[1] TRUE

# log-normalize
spe_pb <- scuttle::computeLibraryFactors(spe_pb)
spe_pb <- scuttle::logNormCounts(spe_pb)


spe_pb
# class: SpatialExperiment 
# dim: 366 308 
# metadata(22): Samples Samples ... Samples Samples
# assays(2): counts logcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(308): MSN_A.Br6660_NAc1_580 D1_Island_A.Br6660_NAc1_580 ...
# Microglia.Br6436_Nac_11_5650 Endothelial.Br6436_Nac_11_5650
# colData names(65): Sample Barcode ... ncells sizeFactor
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id


lib_sizes <- colSums(assay(spe_pb, "counts"))

df <- data.frame(
  Pseudobulk = colnames(spe_pb),
  Sample     = spe_pb$Sample,
  Annotation = spe_pb$Annotation,
  TotalCounts = lib_sizes,
  SizeFactor  = sizeFactors(spe_pb)   # from computeLibraryFactors()
)

subset(df,subset=(Sample == "Br6436_Nac_6_3150"))

spe

ggplot(df, aes(x = TotalCounts, y = SizeFactor, label = Annotation)) +
  geom_point() +
  geom_text(vjust = -0.5, size = 3, check_overlap = TRUE) +
  scale_x_log10() + scale_y_log10() +   # log scales help if ranges are large
  labs(title = "Library size vs Size factor",
       x = "Total counts (per pseudobulk sample)",
       y = "Size factor (from computeLibraryFactors)") +
  theme_minimal()

#Go ahead and add all of the gene names to the 
for(i in rownames(spe_pb)){
  print(i)
  colData(spe_pb)[[i]] <- assay(spe_pb,"logcounts")[i,]
}

#Make boxplots for some specific genes
genes <- c("NPY","CHAT","SLC18A3","OLIG1","SLC1A2")

#make boxplots to help with annotation 
for(i in genes){
  print(i)
  p <- ggplot(as.data.frame(colData(spe_pb)),aes(x = Annotation,y = .data[[i]],fill = Banksy_Cluster)) +
    geom_boxplot(outlier.shape = NA) + #Don't plot outlier twice
    geom_jitter(alpha = 0.6) +
    ggtitle(i) +
    labs(y = "logcounts") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 45, hjust =1 ))
  ggsave(p,filename = here("plots","05_clustering","GeneExpression",
                           "RetreatPoster",paste0(i,".pdf")),
         height = 5, width = 6.5)
}

