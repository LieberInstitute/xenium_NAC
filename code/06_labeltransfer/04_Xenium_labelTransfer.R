#Perform label transfer with Seurat
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5

library(Seurat)
library(here)

#Raise limit
options(future.globals.maxSize = 16 * 1024^3) #16 GB


#Load the seurat objects
seurat_6436 <- readRDS(here("processed-data", "06_label_transfer", 
                            "Objects", "seurat_6436.Rds"))

seurat_anno <- readRDS(here(here("processed-data", "06_label_transfer", 
                                 "Objects", "seurat_anno.Rds")))

#Perform label transfer
#Will use top 30 PCs for initial trial of Seurat label transfer
#Identify the transfer anchors
message(paste0("Finding transfer anchors - ", Sys.time()))
anchors <- FindTransferAnchors(reference = seurat_anno, 
                               query = seurat_6436, 
                               dims = 1:30,
                               reference.reduction = "seurat_pca")

#Perform the label transfer
message(paste0("Transferring Data - ", Sys.time()))
predictions <- TransferData(anchorset = anchors, 
                            refdata = seurat_anno$CellType, 
                            dims = 1:30)

#Add predictions to the object
message(paste0("Adding metadata - ", Sys.time()))
seurat_6436 <- AddMetaData(seurat_6436, 
                           metadata = predictions)

#Save the predictions  and updated seurat objectas an RDS file. 
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions,file = here("processed-data","06_label_transfer","seurat_6436_predictions.RDS"))
saveRDS(object = seurat_6436,file = saveRDS(seurat_anno, here("processed-data", "06_label_transfer", 
                                                              "Objects", "seurat_6436_predictions.Rds")))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
