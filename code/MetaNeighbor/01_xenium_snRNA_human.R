# Use MetaNeighbor to compare human snRNA-seq and Xenium
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(MetaNeighbor)
library(here)

# Load human spe object (Xenium)
spe <- readRDS(here("processed-data", "02_build_spe", "SPEs", "spe_celltype_v2.Rds"))

# Load human sce object (snRNA-seq)
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

# Remove ambiguous neurons and genes with zero counts
sce <- sce[, sce$CellType.Final != "Neuron_Ambig"]
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

# Uniquify gene names to avoid duplicates
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(
  rowData(sce)$gene_id, rowData(sce)$gene_name
)
rownames(sce) <- rowData(sce)$Symbol.uniq

##### Find intersection of genes and subset both objects #####
shared_genes <- intersect(rownames(sce), rownames(spe))
sce <- sce[shared_genes, ]
spe <- spe[shared_genes, ]

# Verify gene order is identical
stopifnot(identical(rownames(sce), rownames(spe)))

##### Add dataset labels #####
sce$dataset <- "snRNA"
spe$dataset <- "xenium"

##### Harmonize colData — use the SAME column names in both #####
# snRNA-seq
colData(sce) <- colData(sce)[, c("Sample", "CellType.Final", "dataset")]
colnames(colData(sce)) <- c("Sample", "CellType", "dataset")

# Xenium (spatial)
colData(spe) <- colData(spe)[, c("Sample", "CellTypes", "dataset")]
colnames(colData(spe)) <- c("Sample", "CellType", "dataset","sample_id") #SpatialExperiment forces a sample_id column


##### Harmonize assays #####
# Keep only counts and logcounts in snRNA-seq
assays(sce) <- assays(sce)[c("counts", "logcounts")]

# Keep counts and rename nucleus_normcounts -> logcounts in Xenium
assays(spe) <- assays(spe)[c("counts", "nucleus_normcounts")]
assayNames(spe)[2] <- "logcounts"

##### Rebuild as clean SingleCellExperiment objects #####
# This avoids rowRanges / seqnames incompatibility issues
# and the cbind problems between SCE and SPE classes
sce_clean <- SingleCellExperiment(
  assays = assays(sce),
  colData = colData(sce)
)
rownames(sce_clean) <- rownames(sce_clean)

spe_clean <- SingleCellExperiment(
  assays = assays(spe),
  colData = colData(spe)[,c("Sample","CellType","dataset")]
)
rownames(spe_clean) <- rownames(spe)

rm(sce,spe)

##### Combine into a single object #####
stopifnot(identical(nrow(sce_clean), nrow(spe_clean)))
stopifnot(identical(rownames(sce_clean), rownames(spe_clean)))

combo <- cbind(sce_clean, spe_clean)
rownames(combo) <- shared_genes
combo

#   Run unsupervised MetaNeighbor
aurocs <- MetaNeighborUS(
  var_genes = rownames(combo),
  dat = combo,
  study_id = combo$dataset,
  cell_type = combo$CellType,
  fast_version = TRUE,
  one_vs_best = TRUE, symmetric_output = FALSE
)

#Save the output 
saveRDS(aurocs,file = here("processed-data","MetaNeighbor","human_xenium_snRNA_aurocs.Rds"))


pdf(file = here("plots","MetaNeighbor",
                "Xenium_snRNA_HumanOnly_MetaNeighbor_Heatmap.pdf"),
    width = 12, height = 12)
plotHeatmap(aurocs, 
            cex = 0.75)

dev.off()

sessionInfo()
