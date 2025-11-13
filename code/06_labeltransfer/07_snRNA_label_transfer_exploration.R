# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(escheR)
library(Seurat)
library(here)

#Load in the count normalized spe object
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

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

#load in the snRNA-seq predictions
Br6660_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_anno_predictions.Rds"))
Br6436_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_6436_snRNA_predictions.Rds"))
predictions <- rbind(Br6660_predictions,Br6436_predictions)


#Sanity checks before adding 
identical(nrow(predictions),ncol(spe))
#[1] TRUE

identical(rownames(predictions),rownames(colData(spe)))
#[1] TRUE

identical(rownames(predictions),colnames(spe))
#[1] TRUE

#Prediction dataframes have the predicted id and the weighted probability score for each designation. 
#Add that information to the colData of the spe object. prediction.score.max
colData(spe)$snRNA_predicted_CellType <- predictions[,"predicted.id"]
colData(spe)$predicted_max_score <- predictions[,"prediction.score.max"]

table(colData(spe)$snRNA_predicted_CellType)
# Astrocyte_A Astrocyte_B  DRD1_MSN_A  DRD1_MSN_B  DRD1_MSN_C  DRD1_MSN_D 
# 556353      150080     1109211       78429       55237       26038 
# DRD2_MSN_A  DRD2_MSN_B Endothelial   Ependymal  Excitatory       Inh_A 
# 306335       19920      426859       25396      144365       27353 
# Inh_B       Inh_C       Inh_D       Inh_E       Inh_F   Microglia 
# 13465       27044       21104       26348      146711      297154 
# Oligo         OPC 
# 1198472      193499 

#load cluster colors from snRNA-seq paper. 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE)
# Loading objects:
#   cluster_cols

#Remove neuron ambig
cluster_cols <- cluster_cols[-14]
cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C  DRD1_MSN_D       Inh_B       Inh_C       Inh_D       Inh_E 
# "#2A4BF9"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281"   "#0DFBFA" 
# Excitatory       Inh_F 
# "#5C6300"     "black" 


spe$snRNA_predicted_CellType <- factor(x = spe$snRNA_predicted_CellType,
                                       levels = c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                                  "DRD2_MSN_A","DRD2_MSN_B",
                                                  "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                                  "Excitatory",
                                                  "Oligo","OPC",
                                                  "Astrocyte_A","Astrocyte_B","Ependymal",
                                                  "Microglia",
                                                  "Endothelial"))


library(ggplot2)

max_score_box <- ggplot(colData(spe),aes(x = snRNA_predicted_CellType, 
                                         y = predicted_max_score,
                                         fill = snRNA_predicted_CellType)) +
  geom_boxplot() +
  geom_hline(yintercept = 0.50) +
  theme_bw() +
  scale_fill_manual(values = cluster_cols) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  labs(x = "Predicted Cell Type",
       y = "Max Score")
ggsave(max_score_box,filename = here("plots","06_label_transfer","snRNA_LabelTransfer",
                                     "Max_score_boxplot.pdf"))
