#Goal: Use spatransfer to perform label transfer between snRNA-seq objects and the xenium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

#Load libraries 
library(SingleCellExperiment)
library(SpatialExperiment)
library(nmfLabelTransfer)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(here)

#Load and prep visium from Ravichandran, Bach et al 2025
vis_spe <- loadHDF5SummarizedExperiment("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/03-filter_normalize_spe/spe_filtered_hdf5")
rownames(vis_spe) <- rowData(vis_spe)$gene_name
rownames(vis_spe) <- make.names(rownames(vis_spe),unique = TRUE)

# Add the final spatial domains
clusters_file <- "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data//07_spatial_domains/01_precast/nnSVG_precast/final_clusters/precast_clusters.csv"
vis_spe[["spatial_domains"]] = colData(vis_spe) |>
  as_tibble() |>
  left_join(read.csv(clusters_file), by = 'key') |>
  pull(cluster) |>
  as.factor()

vis_spe

#Convert the matrices
message(paste0("Converting the counts matrix -",Sys.time()))
counts(vis_spe) <- as(counts(vis_spe),"dgCMatrix")

message(paste0("Converting the logcounts matrix -",Sys.time()))
logcounts(vis_spe) <- as(logcounts(vis_spe),"dgCMatrix")


#Load and prep the spe object containing the xenium NAc data
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))
logcounts(spe) <- assay(spe,"nucleus_normcounts")
rowData(spe)$gene_name <- rownames(spe)

spe

#Target is the spe object, source is the visium  object
transfer_res <- transfer_labels(targets = spe,
                                source = vis_spe,
                                assay = "logcounts",
                                annotationsName = "spatial_domains",
                                technicalVarName = "Sample",
                                seed = 1800,
                                save_nmf = TRUE,
                                nmf_path  = here("processed-data","06_label_transfer","spatransfer","visium_to_Xenium.Rds"),
                                k = NULL, #Set null to run cross-validation to identify optimal number of factors
                                tol = 1e-5,
                                alpha = 0) #pure ridge regression

# save the results
saveRDS(transfer_res, file = here("processed-data", "06_label_transfer", "snRNA_to_Xenium_target_predicitons.Rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
