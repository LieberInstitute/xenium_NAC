#Goal: Use RCTD to transfer cell types from snRNA to VisiumHD
library(SingleCellExperiment)
library(SpatialExperiment)
library(zellkonverter)
library(HDF5Array)
library(spacexr)
library(Matrix)
library(here)

###### Prep snRNA-seq data
#Load snRNA-seq data. 
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")
sce

#From: https://github.com/LieberInstitute/spatialdACC/blob/main/code/04_build_sce/04_raw_sce.R
# -> To avoid getting into sticky situations, let's 'uniquify' these names:
#This allows us to search for genes by their common gene name. If a common gene name is duplicated, then the gene will be
#in the format of genesymbol_ensemblgeneID
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq

# Create SummarizedExperiment.
reference_se <- SummarizedExperiment(
    assays = list(counts = assay(sce,"counts")),
    colData = colData(sce)$CellType.Final
)

reference_se

###### Prep spatial data
#Read in the filtered spe object
spe_filtered_path <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_path)


spatial_spe <- SpatialExperiment(
    assays = list(counts = assay(spe, "counts")),
    spatialCoords = cbind(
        x = spe$pxl_col_in_fullres,
        y = spe$pxl_row_in_fullres
    )
)

spatial_spe

############# Preprocess and run RCTD
message("Creating Rctd object - ",Sys.time())
rctd_data <- createRctd(spatial_spe, reference_se,cell_type_col = "X")

message("Running Rctd - ",Sys.time())
results <- runRctd(rctd_data, rctd_mode = "doublet")

saveRDS(results,here("processed-data","HD_Full_Analysis","LabelTransfer","snRNA_VisiumHD_RCTD.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
