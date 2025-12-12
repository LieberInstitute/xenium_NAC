#Add snRNA
#cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
library(Seurat)
library(here)

#Read in the annotated spe containing just Br6660
message(paste0("Reading in spe_Br6660_anno - ",Sys.time()))
spe_Br6660_anno <- readRDS(here("processed-data", "05_Clustering", "SPEs", "spe_annotated.Rds"))

#Load suerat object for Br6436
message(paste0("Reading in seurat_6436 - ",Sys.time()))
seurat_6436 <- readRDS(here("processed-data", "06_label_transfer","Objects", "seurat_6436_predictions_object.Rds"))

#Make a dataframe consisting of cell_id and banksy cluster
message(paste0("Making the banksy cell type dataframe - ",Sys.time()))
Banksy_CellTypes <- data.frame(cell_id = c(colnames(spe_Br6660_anno),colnames(seurat_6436)),
                               CellType = c(spe_Br6660_anno$CellType,seurat_6436$predicted.id))

#Save as RDS file
message(paste0("Save the RDS - ",Sys.time()))
saveRDS(Banksy_CellTypes, file = here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))


#Reproducibility
sessionInfo()
