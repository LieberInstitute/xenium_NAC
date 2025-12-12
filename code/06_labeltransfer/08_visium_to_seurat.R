library(SpatialExperiment)
library(HDF5Array)
library(Seurat)
library(dplyr)
library(here)

#Pull the genes
xen_genes <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds")) |> rownames()

# #Read in spe that contains the raw counts matrix. 
spe <- loadHDF5SummarizedExperiment("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/03-filter_normalize_spe/spe_filtered_hdf5")
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe),unique = TRUE)

#Are all of the xenium genes in the 
table(rownames(spe) %in% xen_genes)

#Which genes aren't in the spe object
setdiff(xen_genes,rownames(spe))

#subset for xenium genes
spe <- spe[rownames(spe) %in% xen_genes,]

spe

# Add the final spatial domains
clusters_file <- "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data//07_spatial_domains/01_precast/nnSVG_precast/final_clusters/precast_clusters.csv"
spe[["spatial_domains"]] = colData(spe) |>
  as_tibble() |>
  left_join(read.csv(clusters_file), by = 'key') |>
  pull(cluster) |>
  as.factor()

#Remove any spot with no spatial domains
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]

#Convert the matrices
message(paste0("Converting the counts matrix -",Sys.time()))
counts(spe) <- as(counts(spe),"dgCMatrix")
message(paste0("Converting the logcounts matrix -",Sys.time()))
logcounts(spe) <- as(logcounts(spe),"dgCMatrix")

spe


#Convert to seurat
seurat_visium <- as.Seurat(spe,counts = "counts",data = "logcounts")

#Run PCA 
seurat_visium <- FindVariableFeatures(seurat_visium)
seurat_visium <- ScaleData(seurat_visium)
seurat_visium <- RunPCA(seurat_visium,seed.use = 1652,npcs = 50,reduction.name = "seurat_pca")
seurat_visium <- RunUMAP(seurat_visium,seed.use = 1652, dims= 1:50,reduction = "seurat_pca") 

#Save the object
saveRDS(seurat_visium, here("processed-data", "06_label_transfer", 
                          "Objects", "seurat_visium.Rds"))

#Reproducibility
sessionInfo()
