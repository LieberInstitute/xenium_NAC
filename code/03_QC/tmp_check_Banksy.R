library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)
library(scater)


#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts_R_4.4.Rds"))

spe
# class: SpatialExperiment 
# dim: 541 2375058 
# metadata(11): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
# DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2375058): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(56): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe
# class: SpatialExperiment 
# dim: 366 2375058 
# metadata(11): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2375058): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(56): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

#Add Banksy clusters to the object. 
Banksy_clusters <- read.csv(here("processed-data","05_Clustering","Banksy_clusters.csv"))
head(Banksy_clusters)
# X V1                V2
# 1 1  5 Br6660_NAc1_580_1
# 2 2  5 Br6660_NAc1_580_2
# 3 3  5 Br6660_NAc1_580_3
# 4 4  5 Br6660_NAc1_580_4
# 5 5  5 Br6660_NAc1_580_5
# 6 6  5 Br6660_NAc1_580_6

table(Banksy_clusters$V1)
# 1      2      3      4      5      6      7      8      9     10     11 
# 564628 417918 323302 317806 237014 169281  89446  60439  55478  51784  37109 
# 12     13 
# 31097  19756 

#X is just rownames, remove it. 
Banksy_clusters <- Banksy_clusters[,-1]

#Rename the columns
colnames(Banksy_clusters) <- c("Banksy_Cluster","key")

identical(rownames(colData(spe)),Banksy_clusters$key)
#[1] TRUE

#Just cbind the clusters into the colData
colData(spe)$Banksy_Cluster <- as.factor(Banksy_clusters$Banksy_Cluster)

###########Complext heatmap of basic markers to classify Banksy_Cluster
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$Banksy_Cluster)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
Panel <- readxl::read_excel(here("processed-data","Files_For_Upload","NAc_Xenium_Panel_Final_withNotes.xlsx"))
markers_all <- Panel$Gene[order(Panel$Cell_Type)]

marker_labels <-  Panel$Cell_Type[order(Panel$Cell_Type)]
marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c(`D1-MSN` = "#332288",
                                  `D2-MSN` = "#117733",
                                  `D1-Islands` = "#44AA99",
                                  MSNs = "#88CCEE",
                                  `Inhibitory-Neuron` = "#DDCC77",
                                  `Spatial-Marker` = "#8C564B",
                                  Astrocytes = "#D81B60",
                                  Endothelial = "#CC6677",
                                  Ependymal = "#AA4499",
                                  Microglia = "#882255",
                                  Other = "#FE6100",
                                  Excitatory = "#FFB000",
                                  Neurons = "#1E88E5",
                                  OT = "gray",
                                  `Lateral-Septum` = "black",
                                  `Diagonal-Band` = "white",
                                  `OUD-GWAS` = "green",
                                  `SCZ-GWAS` = "blue"))

correct_order <- unique(Panel$Cell_Type[order(Panel$Cell_Type)])

colors_markers$marker <- colors_markers$marker[correct_order]

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
dat <- assay(spe,"nucleus_normcounts")
#rownames(dat) <- rowData(sce)$gene_name
dim(dat)

dat <- dat[markers_all,]
dim(dat)
dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)
max(hm_mat)
#[1] 3.321265

min(hm_mat)
#[1] -1.944437


#Name the clusters
cluster_annotation_df <- data.frame(Banksy_cluster = levels(spe$Banksy_Cluster),
                                    Annotation = c(
                                                "Astrocyte-A",
                                                "Astrocyte-B",
                                                "MSN-1",
                                                "Other_A",
                                                "WM",
                                                "Excitatory",
                                                "Endothelial",
                                                "Microglia",
                                                "Inh_A",
                                                "D1_Island_A",
                                                "Other_B",
                                                "Ependymal",
                                                "D1_Island_B",
                                                "Chat_Inh"
                                                ))

spe$Annotation <- cluster_annotation_df$Annotation[match(spe$Banksy_Cluster,
                                                         cluster_annotation_df$Banksy_cluster)]


#Read in the singleR results. Maybe that will provide some information as to what "Other is"
#Load all of the SingleR results
filenames <- list.files(here("processed-data","06_label_transfer"),full.names = TRUE)
file_list <- lapply(filenames, function(x) get(load(x)))

res <- do.call(what = rbind,file_list)
dim(res)
#[1] 2375058       4

#Create a key column 
res$key <- rownames(res)
colData(spe)$key <- rownames(colData(spe))

#Put res in the order of the spe object
res <- res[match(colData(spe)$key,res$key),]

#Some sanity checks that it worked. 
identical(res$key,spe$key)
#[1] TRUE

identical(res$key,rownames(colData(spe)))
#[1] TRUE

#Add the labels and pruned labels to the object
#Just use the labels column 
colData(spe) <- cbind(colData(spe),res[,c("labels","pruned.labels")])

labels_annotation <- as.data.frame.matrix(table(spe$pruned.labels,spe$Annotation))

#Calculate percentages 
labels_annotation_pct <- sweep(labels_annotation,MARGIN = 2,colSums(labels_annotation),"/") * 100
labels_annotation_pct$snRNAseq_cluster <- rownames(labels_annotation_pct)
labels_annotation_pct_melt <- reshape2::melt(labels_annotation_pct)
#Using snRNAseq_cluster as id variables

load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)
# Loading objects:
#   cluster_cols

cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C Neuro_Ambig  DRD1_MSN_D       Inh_B       Inh_C       Inh_D 
# "#2A4BF9"      "grey"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281" 
# Inh_E  Excitatory       Inh_F 
# "#0DFBFA"   "#5C6300"     "black" 

cluster_cols <- cluster_cols[-14]


p1 <- ggplot(data = labels_annotation_pct_melt,aes(x = variable,y=value,fill =snRNAseq_cluster)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Annotation",
       y = "% of Annotation",
       fill = "snRNA-seq\nCluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = here("plots","05_clustering","tmp","Annotation_percentage_snRNA-seq_cluster_bar.pdf"),plot = p1)


#Change other to Inh_B
# spe$Annotation <- ifelse(spe$Annotation == "Other_A",
#                          "Inh_B",
#                          spe$Annotation)

#Create a color palette that makes sense. 
#make some new brain colors
cluster_cols <- Polychrome::createPalette(length(unique(spe$Annotation)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(cluster_cols) <- unique(spe$Annotation)

#Switch colors for Excitatory and D1_Island_B
# cluster_cols[1] <- "#B2E0FF" 
# cluster_cols[13] <- "#D70D5D"

#Plot with escheR 
for (i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[, spe$Sample == i]
  
  p <-make_escheR(sub_spe) |>
    add_fill("Annotation")+
    ggtitle(i) +
    scale_fill_manual(values = cluster_cols) +
    theme(element_text(hjust = 0.5))

  png(filename = here("plots","05_clustering","tmp",paste0(i,".png")),
      units = "in",height = 20,width = 20,res = 300)
  print(p)
  dev.off()
  
}

p2 <- plotColData(spe, x = "Annotation", y = "sum") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p2, filename = here("plots","05_clustering","tmp","umi_Banksy_violin.png"))   


p2 <- plotColData(spe, x = "Annotation", y = "sum") +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p2, filename = here("plots","05_clustering","tmp","umi_log_Banksy_violin.png"))   




