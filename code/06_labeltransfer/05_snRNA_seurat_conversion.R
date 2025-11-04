#Convert single cell object to seurat
library(SingleCellExperiment)
library(sessioninfo)
library(Seurat)
library(here)

#Load spe object containing normalized coutns
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")
sce

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq


#Remove Neuron_Ambig and 0 count genes
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

sce

#Convert to seurat
seurat_snrna <- as.Seurat(sce,counts = "counts",data = "logcounts")

#Run PCA
seurat_snrna <- FindVariableFeatures(seurat_snrna)
seurat_snrna <- ScaleData(seurat_snrna)
seurat_snrna <- RunPCA(seurat_snrna,seed.use = 1318,npcs = 50,reduction.name = "seurat_pca")

#print the object
seurat_snrna

#Save the object
saveRDS(seurat_snrna, here("processed-data", "06_label_transfer", 
                           "Objects", "seurat_snrna.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
