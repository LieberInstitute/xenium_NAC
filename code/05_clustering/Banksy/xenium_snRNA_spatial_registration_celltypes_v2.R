#Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(here)

#Load the xenium object 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))

spe

#Add the Banksy cell types to the domain. 
#Load the annotated banksy cell types
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

stopifnot(identical(Banksy_celltypes$cell_id,colnames(spe)))

spe$Banksy_celltypes <- Banksy_celltypes$CellType

## Perform the spatial registration
xen_modeling_results <- registration_wrapper(
  sce = spe,
  var_registration = "CellTypes",
  var_sample_id = "Sample",
  gene_ensembl = "ID",
  gene_name = "Symbol"
)

saveRDS(object = xen_modeling_results,
        file   = here("processed-data","spatial_registration","xen_modeling_results_CellTypes_v2.Rds"))

#Load the sce modeling results
sce_modeling_results <- readRDS(here("processed-data","spatial_registration","sce_modeling_results.Rds"))


cor_res <- layer_stat_cor(
  stats = xen_modeling_results$enrichment, #query
  modeling_results = sce_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data",
                      "spatial_registration","xenium_snRNA_registration_results_celltypes_v2.Rds"))


pdf(here("plots","spatial_registration","xenium_snRNA_registration_v2_celltypes.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
