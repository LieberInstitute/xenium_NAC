# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(Seurat)
library(here)

#Raise limit
options(future.globals.maxSize = 16 * 1024^3) #16 GB

#Load the visium seurat object from script 08
seurat_visium <- readRDS(here("processed-data", "06_label_transfer", 
                              "Objects", "seurat_visium.Rds"))

seurat_visium

#load the annotated seurat object containing banksy clusters first. 
seurat_6436 <- readRDS(here("processed-data", "06_label_transfer", 
                            "Objects", "seurat_6436.Rds"))

#Perform label transfer
#Identify the transfer anchors
message(paste0("Finding transfer anchors - ", Sys.time()))
anchors <- FindTransferAnchors(reference = seurat_visium, 
                               query = seurat_6436, 
                               dims = 1:50,
                               reference.reduction = "seurat_pca",
                               reduction = "rpca")

#Perform the label transfer
message(paste0("Transferring Data - ", Sys.time()))
predictions <- TransferData(anchorset = anchors, 
                            refdata = seurat_visium$spatial_domains, 
                            dims = 1:50)

#Add predictions to the object
message(paste0("Adding metadata - ", Sys.time()))
seurat_6436 <- AddMetaData(seurat_6436, 
                           metadata = predictions)

#Save the predictions  and updated seurat objectas an RDS file. 
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions,
        file = here("processed-data",
                    "06_label_transfer",
                    "seurat_6436_spatial_domains_predictions.Rds"))
saveRDS(object = seurat_6436,
        file = here("processed-data", 
                    "06_label_transfer",
                    "Objects", 
                    "seurat_6436_spatial_domains_predictions_object.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
