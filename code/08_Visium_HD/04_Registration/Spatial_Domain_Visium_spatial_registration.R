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


vis_modeling_results <- readRDS(here("processed-data",
				     "HD_Full_Analysis",
                                     "Spatial_Registration",
				     "visium_modeling_results_spatial_domains.Rds"))



#
cor_res <- layer_stat_cor(
  stats = HD_modeling_results$enrichment, #query
  modeling_results = vis_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","HD_Visium_results.Rds"))


pdf(here("plots","HD_Full_Analysis","Spatial_Registration","visium_spatial_domains_HD.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
