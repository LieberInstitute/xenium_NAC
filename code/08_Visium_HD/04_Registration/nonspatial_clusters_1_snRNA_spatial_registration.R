#Perform spatial registration between HD and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(dplyr)
library(here)

#Read in the filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

spe

#Read in the nonspatial clusters
nonspatial_clusters <- readRDS(here("processed-data","HD_Full_Analysis",
                                    "Banksy","NonSpatial_colData",
                                    "nonspatial_banksy_res_1.Rds"))

stopifnot(identical(rownames(nonspatial_clusters),colnames(spe)))

spe$non_spatial_1 <- nonspatial_clusters$clust_M0_lam0_k50_res1

# Get cell metadata with barcodes
cell_df <- data.frame(barcode = colnames(spe),
                      sample = spe$sample_id,
                      cluster = spe$non_spatial_1)

# Sample 33% per sample per cluster
set.seed(2105)
keep <- cell_df %>%
  group_by(sample,cluster) %>%
  slice_sample(prop = 0.33) %>%
  pull(barcode)

# Subset the object
spe_sub <- spe[, keep]
spe_sub

## Perform the spatial registration
HD_modeling_results <- registration_wrapper(
  sce = spe_sub,
  var_registration = "non_spatial_1",
  var_sample_id = "sample_id",
  gene_ensembl = "gene_id",
  gene_name = "gene_name"
)

saveRDS(object = HD_modeling_results,
        file   = here("processed-data","HD_Full_Analysis",
                      "Spatial_Registration","HD_modeling_results_non_spatial_1.Rds"))


#Load sce modeling results
## Perform the spatial registration
sce_modeling_results <- readRDS(here("processed-data","HD_Full_Analysis","Spatial_Registration","sce_modeling_results.Rds"))

cor_res <- layer_stat_cor(
  stats = HD_modeling_results$enrichment, #query
  modeling_results = sce_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   =here("processed-data","HD_Full_Analysis",
                     "Spatial_Registration","HD_nonspatial_snRNA_results.Rds"))


pdf(here("plots","HD_Full_Analysis","Spatial_Registration","snRNA_CellType.Final_nonspatial_clusters_1.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
