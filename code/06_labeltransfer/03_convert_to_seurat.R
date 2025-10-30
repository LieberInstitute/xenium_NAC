library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(Seurat)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

#Subset to Br6436
spe_6436 <- spe[,spe$Donor == "Br6436"]
logcounts(spe_6436) <- assay(spe_6436,"nucleus_normcounts")

#Convert to seurat
seurat_6436 <- as.Seurat(spe_6436,counts = "counts",data = "logcounts")

#Save the object
saveRDS(seurat_6436, here("processed-data", "06_label_transfer", 
                          "Objects", "seurat_6436.Rds"))

#Load spe object containing normalized coutns
spe_anno <- readRDS(here("processed-data", "05_Clustering", "SPEs", "spe_annotated.Rds"))
logcounts(spe_anno) <- assay(spe_anno,"nucleus_normcounts")

#Convert to seurat
seurat_anno <- as.Seurat(spe_anno,counts = "counts",data = "logcounts")

#Save the object
saveRDS(seurat_anno, here("processed-data", "06_label_transfer", 
                          "Objects", "seurat_anno.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
