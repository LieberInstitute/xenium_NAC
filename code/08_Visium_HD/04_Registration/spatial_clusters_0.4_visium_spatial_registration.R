#Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(dplyr)
library(here)

#Read in the Visium data
# #Read in spe that contains the raw counts matrix. 
spe <- loadHDF5SummarizedExperiment("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/05_harmony_BayesSpace/03-filter_normalize_spe/spe_filtered_hdf5")
rownames(spe) <- rowData(spe)$gene_name
rownames(spe) <- make.names(rownames(spe),unique = TRUE)

# Add the final spatial domains
clusters_file <- "/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data//07_spatial_domains/01_precast/nnSVG_precast/final_clusters/precast_clusters.csv"
spe[["spatial_domains"]] = colData(spe) |>
  as_tibble() |>
  left_join(read.csv(clusters_file), by = 'key') |>
  pull(cluster) |>
  as.factor()

#Remove any spot with no spatial domains
spe <- spe[ ,!is.na(spe[["spatial_domains"]])]


spe


logcounts(spe)[1:5,1:5]


## Perform the spatial registration
vis_modeling_results <- registration_wrapper(
  sce = spe,
  var_registration = "spatial_domains",
  var_sample_id = "sample_id",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)


saveRDS(object = vis_modeling_results,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","visium_modeling_results_spatial_domains.Rds"))



#Load HD_modeling_results
HD_modeling_results <- readRDS(here("processed-data","HD_Full_Analysis",
                                    "Spatial_Registration","HD_modeling_results_spatial_clusters_0.4.Rds"))

#
cor_res <- layer_stat_cor(
  stats = HD_modeling_results$enrichment, #query
  modeling_results = vis_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","HD_Visium_results.Rds"))


pdf(here("plots","HD_Full_Analysis","Spatial_Registration","visium_spatial_domains_HD_spatial_clusters_0.4.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
