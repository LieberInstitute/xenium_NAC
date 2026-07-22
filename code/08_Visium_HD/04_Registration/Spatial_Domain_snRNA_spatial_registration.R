# Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(dplyr)
library(here)

sfe_filtered_dir <- here("processed-data","HD_Full_Analysis","sfe_spatial_annotated")
sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

# Get cell metadata with barcodes
cell_df <- data.frame(barcode = colnames(sfe),
                      sample = sfe$sample_id,
                      Domain = sfe$Spatial_Domain)

# Sample 33% per sample per cluster
set.seed(1309)
keep <- cell_df %>%
  group_by(sample,Domain) %>%
  slice_sample(prop = 0.33) %>%
  pull(barcode)

# Subset the object
sfe_sub <- sfe[, keep]
sfe_sub

## Perform the spatial registration
HD_modeling_results <- registration_wrapper(
  sce = sfe_sub,
  var_registration = "Spatial_Domain",
  var_sample_id = "sample_id",
  gene_ensembl = "ID",
  gene_name = "Symbol"
)

saveRDS(object = HD_modeling_results,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","HD_modeling_results_Spatial_Domain.Rds"))


#Load sce object containing normalized coutns
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Remove the neuronal ambiguous population 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

sce

## Perform the spatial registration
sce_modeling_results <- registration_wrapper(
  sce = sce,
  var_registration = "CellType.Final",
  var_sample_id = "Sample",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

saveRDS(object = sce_modeling_results,
        file   = here("processed-data","HD_Full_Analysis","Spatial_Registration","sce_modeling_results.Rds"))

cor_res <- layer_stat_cor(
  stats = HD_modeling_results$enrichment, #query
  modeling_results = sce_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   =here("processed-data","HD_Full_Analysis",
                     "Spatial_Registration","HD_snRNA_results.Rds"))


pdf(here("plots","HD_Full_Analysis","Spatial_Registration",
         "snRNA_CellType.Final_Spatial_Domain.pdf"),
    height = 8, 
    width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
