#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5

library(SpatialExperiment)
library(here)
library(Seurat)
library(zellkonverter)
library(ggplot2)
library(HDF5Array)

###### Read in the filtered spe object (bin2cell segmented object)
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)
spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

spe

# Raw counts
spe_counts <- spe
assays(spe_counts) <- list(
  counts = assay(spe_counts, "counts")
)
writeH5AD(
    spe_counts,
    file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_counts.h5ad"),
    X_name     = "counts"
)

# logcounts
spe_log <- spe
assays(spe_log) <- list(
  logcounts = assay(spe_log, "logcounts")
)
writeH5AD(
    spe_log,
    file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_logcounts.h5ad"),
    X_name     = "logcounts"
)

# Full SPE with all assays and colData
writeH5AD(
    spe,
    file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_full.h5ad"),
    X_name     = "counts"
)


###### Read in the filtered sfe object (SpaceRanger cell segmentation output with snRNA labels)
library(SpatialFeatureExperiment)
# sfe_path      <- here("processed-data", "HD_Full_Analysis", "sfe_with_labels")
# sfe <- loadHDF5SummarizedExperiment(sfe_path)

# # Raw counts
# sfe_counts <- sfe
# assays(sfe_counts) <- list(
#   counts = assay(sfe_counts, "counts")
# )
# writeH5AD(
#     sfe_counts,
#     file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_sfe_counts.h5ad"),
#     X_name     = "counts"
# )

# # logcounts
# sfe_log <- sfe
# assays(sfe_log) <- list(
#   logcounts = assay(sfe_log, "logcounts")
# )
# writeH5AD(
#     sfe_log,
#     file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_sfe_logcounts.h5ad"),
#     X_name     = "logcounts"
# )

###### Read in the filtered sfe object (SpaceRanger cell segmentation output with spatial domain labels and snRNA labels)
sfe_domain_path      <- here("processed-data","HD_Full_Analysis",  "sfe_spatial_annotated")
sfe_domain <- loadHDF5SummarizedExperiment(sfe_domain_path)

# Raw counts
sfe_counts <- sfe_domain
assays(sfe_counts) <- list(
  counts = assay(sfe_counts, "counts")
)
writeH5AD(
    sfe_counts,
    file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_sfe_counts.h5ad"),
    X_name     = "counts"
)

# logcounts
sfe_log <- sfe_domain
assays(sfe_log) <- list(
  logcounts = assay(sfe_log, "logcounts")
)
writeH5AD(
    sfe_log,
    file  = here("processed-data","HD_Full_Analysis","h5ad","VHD_sfe_logcounts.h5ad"),
    X_name     = "logcounts"
)
