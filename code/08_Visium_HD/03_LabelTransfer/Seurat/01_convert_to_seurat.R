library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(Seurat)
library(escheR)
library(Matrix)
library(here)

###### Prep spatial data
##Read in the filtered spe object
#Read in the filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

spe

#Convert counts to sparse matrices
counts(spe) <- as(counts(spe), "dgCMatrix")
logcounts(spe) <- as(logcounts(spe), "dgCMatrix")

#Convert to seurat
message("Converting spe to Seurat object - ", Sys.time())
seurat_spe <- as.Seurat(spe,counts = "counts",data = "logcounts")

#Run PCA 
message("Finding variable features - ", Sys.time())
seurat_spe <- FindVariableFeatures(seurat_spe)

message("Scaling Data - ", Sys.time())
seurat_spe <- ScaleData(seurat_spe)


message("Running PCA - ", Sys.time())
seurat_spe <- RunPCA(seurat_spe,seed.use = 1744,npcs = 50,reduction.name = "seurat_pca")

#Save the object
saveRDS(seurat_spe, here("processed-data", "HD_Full_Analysis", 
                          "SPEs", "seurat_spe.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
