#Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(dplyr)
library(here)

HD_modeling_results <- readRDS(here("processed-data","HD_Full_Analysis",
                                    "Spatial_Registration","HD_modeling_results_Spatial_Domain.Rds"))

#Load sce object containing normalized coutns
spe_dir <-  here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds")
spe <- readRDS(spe_dir)

## Perform the spatial registration
spe_modeling_results <- registration_wrapper(
  sce = spe,
  var_registration = "CellTypes",
  var_sample_id = "Sample",
  gene_ensembl = "ID",
  gene_name = "Symbol"
)

saveRDS(object = spe_modeling_results,
        file   = here("processed-data","HD_Full_Analysis","Spatial_Registration","xenium_modeling_results.Rds"))

cor_res <- layer_stat_cor(
  stats = HD_modeling_results$enrichment, #query
  modeling_results = spe_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","HD_xenium_results.Rds"))


pdf(here("plots","HD_Full_Analysis","Spatial_Registration",
         "Xenium_CellType_Spatial_Domain.pdf"),
    height = 8, 
    width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
