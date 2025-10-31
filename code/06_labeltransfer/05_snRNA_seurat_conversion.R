library(SingleCellExperiment)
library(sessioninfo)
library(Seurat)
library(here)

#Load spe object containing normalized coutns
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Convert to seurat
seurat_snrna <- as.Seurat(sce,counts = "counts",data = "logcounts")

#convert the HARMONY reduced dimensions to pca for the sake of label transfer. 
seurat_snrna[["pca"]] <- seurat_snrna[["HARMONY"]]
seurat_snrna@reductions[["HARMONY"]] <- NULL

#Save the object
saveRDS(seurat_snrna, here("processed-data", "06_label_transfer", 
                           "Objects", "seurat_snrna.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
