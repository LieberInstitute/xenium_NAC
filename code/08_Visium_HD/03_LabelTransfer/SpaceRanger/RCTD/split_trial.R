# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#Goal: Add RCTD weights and predictions to the sfe
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(Seurat)
library(SPLIT)
library(here)

#Read in the filtered sfe object
#Now save.
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_with_labels"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

sfe


## Get cell type names from first available result
rctd_dir <- here("processed-data", "HD_Full_Analysis", "LabelTransfer","SpaceRanger","RCTD")

samples <- unique(sfe$sample_id)

## Get cell type names AND prediction column names from first available resultS
rctd_files <- vector(mode = "list",length = length(samples))
names(rctd_files) <- samples
for(s in samples){
  print(s)
  rctd_files[[s]] <- readRDS(file.path(rctd_dir, paste0("doublet_snRNA_VisiumHD_RCTD_", s, ".Rds")))
}

rctd_res <- do.call(what = cbind,rctd_files)
sfe <- sfe[,colnames(sfe) %in% colnames(rctd_res)]
rctd_res <- rctd_res[,colnames(sfe)]
sfe <- sfe[,colnames(rctd_res)]

# Bioconductor weights are cell types x pixels; purify wants pixels x cell types
W <- t(as.matrix(assay(rctd_res, "weights")))      # pixels (rows) x cell types (cols)

# primary label per pixel (doublet-mode first_type lives in colData)
primary <- as.character(colData(rctd_res)$first_type)
names(primary) <- colnames(rctd_res)

###### Prep snRNA-seq data (same for all samples)
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Remove the neuronal ambiguous population from further analysis. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]
#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

sce

rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq
ref_counts <- as(assay(sce, "counts"), "dgCMatrix")


reference_se <- SummarizedExperiment(
  assays = list(counts = ref_counts),
  colData = colData(sce)$CellType.Final
)

# reference profile matrix: genes x cell types (mean expression per type).
# Use the SAME reference you gave RCTD. If it was a SummarizedExperiment `ref_se`:
#   ref_profiles <- <genes x cell_types matrix of mean expr>

library(Matrix)
library(SummarizedExperiment)

# ref_se: the SummarizedExperiment given to createRctd()
#   assay "counts" = genes x cells; colData holds the cell-type column
ref_counts <- assays(reference_se)[["counts"]]  # genes x cells
cell_types <- as.factor(colData(reference_se)$X)           # adjust column name
names(cell_types) <- colnames(reference_se)

# 1. library-size normalize each cell (each column sums to 1)
nUMI     <- colSums(ref_counts)
ref_norm <- sweep(ref_counts, 2, nUMI, "/")                  # genes x cells

# 2. mean normalized expression within each cell type -> genes x cell types
ct_levels <- levels(cell_types)
ref_profiles <- vapply(ct_levels, function(ct) {
  Matrix::rowMeans(ref_norm[, cell_types == ct, drop = FALSE])
}, numeric(nrow(ref_norm)))
rownames(ref_profiles) <- rownames(ref_counts)
colnames(ref_profiles) <- ct_levels                          # genes x cell types

## align names/genes so the three inputs are mutually consistent
xe_counts <- counts(sfe)
common    <- intersect(rownames(ref_profiles), rownames(xe_counts))
ref_profiles <- ref_profiles[common, , drop = FALSE]

# keep only cell types that appear in all three, in one consistent order
cts <- Reduce(intersect, list(colnames(W), colnames(ref_profiles), unique(primary)))
W            <- W[, cts, drop = FALSE]
ref_profiles <- ref_profiles[, cts, drop = FALSE]

res_split <- SPLIT::purify(
  counts                = xe_counts,
  reference             = ref_profiles,
  deconvolution_weights = W,
  primary_cell_type     = primary,
  DO_purify_singlets    = TRUE
)

