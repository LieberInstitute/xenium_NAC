#Goal: Convert SpatialExperiment object to h5ad format
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5

library(SpatialExperiment)
library(here)
library(Seurat)
library(zellkonverter)
library(ggplot2)

spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))
spe$cell_id <- colnames(spe)

# Nucleus area normalization
spe_nuc <- spe
assays(spe_nuc) <- list(
  nucleus_normcounts = assay(spe_nuc, "nucleus_normcounts")
)
writeH5AD(
    spe_nuc,
    file  = here("processed-data","02_build_spe","h5ad","spe_NormCounts_nucleus_normcounts.h5ad"),
    X_name     = "nucleus_normcounts",
    colData = c("Sample","Barcode","cell_id","sample_id","Xenium_Run_ID","Slide_ID","Donor", "Age", "Sex", "Race")
)
# Cell area normalization
spe_cell <- spe
assays(spe_cell) <- list(
  cell_normcounts = assay(spe_cell, "cell_normcounts")
)
writeH5AD(
    spe_cell,
    file  = here("processed-data","02_build_spe","h5ad","spe_NormCounts_cell_normcounts.h5ad"),
    X_name     = "cell_normcounts",
    colData = c("Sample","Barcode","cell_id","sample_id","Xenium_Run_ID","Slide_ID","Donor", "Age", "Sex", "Race")
)
# Raw counts
spe_counts <- spe
assays(spe_counts) <- list(
  counts = assay(spe_counts, "counts")
)
writeH5AD(
    spe_counts,
    file  = here("processed-data","02_build_spe","h5ad","spe_NormCounts_counts.h5ad"),
    X_name     = "counts",
    colData = c("Sample","Barcode","cell_id","sample_id","Xenium_Run_ID","Slide_ID","Donor", "Age", "Sex", "Race")
)
# Full SPE with all assays and colData
writeH5AD(
    spe,
    file  = here("processed-data","02_build_spe","h5ad","spe_NormCounts_full.h5ad"),
    X_name     = "counts"
)





