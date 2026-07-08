#Perform label transfer between the snRNA-seq data and the seurat objects for each donor
library(Seurat)
library(here)

#Raise limit
options(future.globals.maxSize = 16 * 1024^3) #16 GB

#load snRNA-seq seurat object
snRNA <- readRDS(here("processed-data", "06_label_transfer", 
                      "Objects", "seurat_snrna.Rds"))

snRNA

spe <- readRDS(here("processed-data", "HD_Full_Analysis", 
                          "SPEs", "seurat_spe.Rds"))

spe

#Identify top 2500 integration features
features <- SelectIntegrationFeatures(object.list = list(snRNA, spe),
                                       nfeatures = 2500)

features <- c(features,"DRD1","RXFP1","OPRM1","SEMA5B","TRHDE","GABRQ","VWC2L","CPNE4","SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP")
features <- unique(features)
length(features)
saveRDS(features[order(features)],file = here("processed-data","HD_Full_Analysis","Seurat","integration_features.Rds"))

#Perform label transfer
#Identify the transfer anchors
message(paste0("Finding transfer anchors - ", Sys.time()))
anchors <- FindTransferAnchors(reference = snRNA,
                               query = spe,
                               dims = 1:30,
			       features = features,
                               reference.reduction = "seurat_pca",
                               reduction = "rpca")

#Perform the label transfer
message(paste0("Transferring Data - ", Sys.time()))
predictions <- TransferData(anchorset = anchors,
                            refdata = snRNA$CellType.Final,
                            dims = 1:30)

#Add predictions to the object
message(paste0("Adding metadata - ", Sys.time()))
spe <- AddMetaData(spe,
                   metadata = predictions)

#Save the predictions  and updated seurat objectas an RDS file. 
message(paste0("Saving objects - ", Sys.time()))
saveRDS(object = predictions,file = here("processed-data","HD_Full_Analysis","Seurat","snRNA_Visium_HD_predictions.Rds"))
saveRDS(object = spe,file = here("processed-data", "HD_Full_Analysis","Seurat","snRNA_Visium_HD_predictions_object.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
