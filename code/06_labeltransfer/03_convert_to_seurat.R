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

#Run PCA 
seurat_6436 <- FindVariableFeatures(seurat_6436)
seurat_6436 <- ScaleData(seurat_6436)
seurat_6436 <- RunPCA(seurat_6436,seed.use = 2051,npcs = 50,reduction.name = "seurat_pca")
seurat_6436 <- RunUMAP(seurat_6436,seed.use = 2051, dims= 1:50,reduction = "seurat_pca") 

#Save the object
saveRDS(seurat_6436, here("processed-data", "06_label_transfer", 
                          "Objects", "seurat_6436.Rds"))

#Load spe object containing normalized coutns
spe_anno <- readRDS(here("processed-data", "05_Clustering", "SPEs", "spe_annotated.Rds"))
logcounts(spe_anno) <- assay(spe_anno,"nucleus_normcounts")

#Convert to seurat
seurat_anno <- as.Seurat(spe_anno,counts = "counts",data = "logcounts")

#Run PCA
seurat_anno <- FindVariableFeatures(seurat_anno)
seurat_anno <- ScaleData(seurat_anno)
seurat_anno <- RunPCA(seurat_anno,seed.use = 2051,npcs = 50,reduction.name = "seurat_pca")
seurat_anno <- RunUMAP(seurat_anno,seed.use = 2051, dims= 1:50,reduction = "seurat_pca")

#Save the object
saveRDS(seurat_anno, here("processed-data", "06_label_transfer", 
                          "Objects", "seurat_anno.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
