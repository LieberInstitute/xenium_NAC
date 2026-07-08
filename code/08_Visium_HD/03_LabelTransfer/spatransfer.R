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
library(here)


message("Prepping snRNA-seq data - ", Sys.time())
#Load and prep snRNA-seq object from Ravichandran, Bach et al 2025
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq
rowData(sce)$gene_name <- rownames(sce)

#Remove Neuron_Ambig and 0 count genes
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

sce

###### Prep spatial data
message("Prepping spatial data - ", Sys.time())
spe_filtered_path <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)
spe <- loadHDF5SummarizedExperiment(spe_filtered_path)


message("Running spatransfer - ", Sys.time())
#Run spatransfer
transfer_res <- transfer_labels(targets = spe,
                                source = sce,
                                assay = "logcounts",
                                annotationsName = "CellType.Final",
                                technicalVarName = "Sample",
                                seed = 1713,
                                save_nmf = TRUE,
                                nmf_path  = here("processed-data","06_label_transfer","spatransfer","snRNA_to_HD_k105.Rds"),
                                k = 105,
                                tol = 1e-5,
                                alpha = 0) #pure ridge regression

# save the results
saveRDS(transfer_res, file = here("processed-data", "HD_Full_Analysis","LabelTransfer","spatransfer","snRNA_to_HD_target_predicitons.Rds"))

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()

