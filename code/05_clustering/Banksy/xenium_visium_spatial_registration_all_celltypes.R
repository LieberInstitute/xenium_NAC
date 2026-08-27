#Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(here)

xen_modeling_results <- readRDS(here("processed-data","spatial_registration","xen_modeling_results_CellTypes_v2.Rds"))


vis_modeling_results <- readRDS(here("processed-data",
				     "HD_Full_Analysis",
                                     "Spatial_Registration",
				     "visium_modeling_results_spatial_domains.Rds"))


cor_res <- layer_stat_cor(
  stats = xen_modeling_results$enrichment, #query
  modeling_results = vis_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data",
                      "spatial_registration","xenium_visium_registration_results.Rds"))


pdf(here("plots","spatial_registration","xenium_visium_registration_celltypes_v2.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()


sessionInfo()
