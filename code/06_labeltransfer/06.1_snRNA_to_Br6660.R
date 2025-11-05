#Perform label transfer between the snRNA-seq data and the seurat objects for each donor
library(Seurat)
library(here)

#Raise limit
options(future.globals.maxSize = 16 * 1024^3) #16 GB

#load snRNA-seq seurat object
snRNA <- readRDS(here("processed-data", "06_label_transfer", 
                      "Objects", "seurat_snrna.Rds"))

snRNA

#load the annotated seurat object containing banksy clusters first. 
seurat_anno <- readRDS(here("processed-data", "06_label_transfer", 
                            "Objects", "seurat_anno.Rds"))

#Perform label transfer
#Will use top 30 PCs for initial trial of Seurat label transfer
#Identify the transfer anchors
message(paste0("Finding transfer anchors - ", Sys.time()))
anchors <- FindTransferAnchors(reference = snRNA, 
                               query = seurat_anno, 
                               dims = 1:30,reduction = "cca",
                               reference.reduction = "seurat_pca")

#Perform the label transfer
message(paste0("Transferring Data - ", Sys.time()))
predictions <- TransferData(anchorset = anchors, 
                            refdata = snRNA$CellType.Final, 
                            dims = 1:30)

#Add predictions to the object
message(paste0("Adding metadata - ", Sys.time()))
seurat_anno <- AddMetaData(seurat_anno, 
                           metadata = predictions)

#Save the predictions  and updated seurat objectas an RDS file. 
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions,file = here("processed-data","06_label_transfer","seurat_anno_predictions.Rds"))
saveRDS(object = seurat_anno,file = here("processed-data", "06_label_transfer","Objects", "seurat_anno_predictions_object.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
