# Convert snRNA-seq SingleCellExperiment to h5ad for scDRS.
#
# Exports raw counts + key colData (Sample, CellType.Final) so that
# scDRS can normalize and score cells.
#
# Usage: sbatch 01_prep_h5ad.sh

library(SingleCellExperiment)
library(zellkonverter)
library(here)

sce_path <- "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds"
out_h5ad <- here("processed-data", "16_scDRS", "snRNA", "h5ad", "snRNA_counts.h5ad")

message("Loading SCE: ", sce_path)
sce <- readRDS(sce_path)
message("Dimensions: ", nrow(sce), " genes x ", ncol(sce), " cells")

# Keep only counts assay for scDRS (it will normalize internally)
sce <- sce[, !is.na(sce$CellType.Final)]
assays(sce) <- list(X = counts(sce))

# Slim colData to what scDRS needs
cd <- colData(sce)[, c("Sample", "CellType.Final"), drop = FALSE]
colData(sce) <- cd

message("Writing h5ad: ", out_h5ad)
writeH5AD(sce, file = out_h5ad)
message("Done.")

sessioninfo::session_info()
