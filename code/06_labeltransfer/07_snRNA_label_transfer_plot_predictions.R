#Add snRNA
#cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
library(escheR)
library(Seurat)
library(here)

#Load in the annotated object --> this contains annotations from Banksy clustering. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

dim(spe)

spe

#load in the snRNA-seq predictions
Br6660_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_anno_predictions.Rds"))
Br6436_predictions <- readRDS(here("processed-data","06_label_transfer","seurat_6436_snRNA_predictions.Rds"))
predictions <- rbind(Br6660_predictions,Br6436_predictions)

dim(predictions)

#Sanity checks before adding 
stopifnot(identical(nrow(predictions),ncol(spe)))


stopifnot(identical(rownames(predictions),rownames(colData(spe))))


stopifnot(identical(rownames(predictions),colnames(spe)))


#Prediction dataframes have the predicted id and the weighted probability score for each designation. 
#Add that information to the colData of the spe object. prediction.score.max
colData(spe)$snRNA_predicted_CellType <- predictions[,"predicted.id"]
colData(spe)$predicted_max_score <- predictions[,"prediction.score.max"]

table(colData(spe)$snRNA_predicted_CellType)


#load cluster colors from snRNA-seq paper. 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)

#Remove neuron ambig
cluster_cols <- cluster_cols[-14]
cluster_cols


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

#Plot the predicted IDs on the tissue sections .
for(i in unique(spe$Sample)){
  print(i)
  spe_sub <- spe[,spe$Sample == i]
  p <- make_escheR(spe_sub) |>
    add_fill("snRNA_predicted_CellType") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here("plots","06_label_transfer","snRNA_LabelTransfer",paste0(i,".png")),
         height = 22, width = 22)
}


###Reproduciblity
sessionInfo()
