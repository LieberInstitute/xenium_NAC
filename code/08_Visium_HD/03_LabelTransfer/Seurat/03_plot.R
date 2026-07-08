#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(Seurat)
library(escheR)
library(here)

#Read in the filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

spe

#Read in the Seurat predictions. 
preds <- readRDS(here("processed-data","HD_Full_Analysis","Seurat","snRNA_Visium_HD_predictions.Rds"))

stopifnot(identical(rownames(preds),colnames(spe)))

spe$Seurat_Prediction <- preds$predicted.id
spe$prediction.score.max <- preds$prediction.score.max

#Load colors
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE) 

#Make a boxplot fo prediction scores
p2 <- ggplot(data = preds,aes(x = predicted.id,y = prediction.score.max,fill = predicted.id)) +
  geom_boxplot() +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none") +
  geom_hline(yintercept = .50,lty = 2)
ggsave(filename = here("plots","HD_Full_Analysis","LabelTransfer","Seurat","prediction_score_boxplot.pdf"),
       plot = p2)

for(sample in unique(spe$sample_id)){
  print(sample)
  spe_sub <- spe[,spe$sample_id == sample]
  p <- make_escheR(spe_sub) |>
    add_fill("Seurat_Prediction") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here("plots","HD_Full_Analysis","LabelTransfer","Seurat",paste0(sample,"_Seurat_Predictions.png")),
         height = 20, width = 14, dpi = 200)
  
}

for(ct in unique(spe$Seurat_Prediction)){
  print(ct)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","Seurat",ct)
  dir.create(path = dir_path)
  for(sample in unique(spe$sample_id)){
    print(sample)
    spe_sub <- spe[,spe$sample_id == sample]
    spe_sub$CellType <- ifelse(spe_sub$Seurat_Prediction == ct,
                               ct,
                               "Other")
    ct_colors <- c("lightgrey",cluster_cols[ct])
    names(ct_colors) <- c("Other",ct)
    p <- make_escheR(spe_sub) |>
      add_fill("CellType") +
      scale_fill_manual(values = ct_colors)
    ggsave(filename = here(dir_path,paste0(sample,".png")),
           plot = p,
           width = 20, height = 14, dpi = 200)
    
  }
}

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
