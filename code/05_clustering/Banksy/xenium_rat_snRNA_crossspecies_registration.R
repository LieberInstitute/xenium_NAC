#Perform spatial registration between xenium and snRNA-seq/Visium data
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(spatialLIBD)
library(orthogene)
library(Seurat)
library(here)

rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

#Convert to 121 orthologs
rat_orthologs <- convert_orthologs(gene_df = counts(rat_sce),
                                   gene_input = "rownames",
                                   gene_output = "dict",
                                   input_species = "rat",
                                   output_species = "human",
                                   non121_strategy = "drop_both_species",
                                   method = "gprofiler")
head(rat_orthologs)

#Generate a dataframe from the above output
rat_ortho_df <- data.frame(rat = names(rat_orthologs),
                           human = unname(rat_orthologs))

#Subset rat_sce to contain only the 121 orthologs
rat_sce <- rat_sce[rownames(rat_sce) %in% rat_ortho_df$rat,]

#Put the rat ortholog datafarme in same order as rat_sce
rat_ortho_df <- rat_ortho_df[match(rownames(rat_sce),rat_ortho_df$rat),]

identical(rownames(rat_sce),rat_ortho_df$rat)

stopifnot(identical(rownames(rat_sce),rat_ortho_df$rat))

#convert the rownames of the rat_sce object to the human orthologs
rownames(rat_sce) <- rat_ortho_df$human

#Load the xenium object
#Load the xenium object 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))

#subset for 121 orthos
spe_sub <- spe[rownames(spe) %in% rownames(rat_sce),]
dim(spe_sub)

#subset rat object for genes in the spe
rat_sce <- rat_sce[rownames(rat_sce) %in% rownames(spe_sub),]

#Put rat sce in same order as spe_sub
rat_sce <- rat_sce[rownames(spe_sub),]
dim(rat_sce)
identical(rownames(spe_sub),rownames(rat_sce))
stopifnot(identical(rownames(spe_sub),rownames(rat_sce)))

#Add a gene name and gene_id column
rowData(rat_sce)[,c("Symbol","ID")] <- rowData(spe_sub)[,c("Symbol","ID")]

#Perform registration
#Rat 
rat_modeling_results <- registration_wrapper(
  sce = rat_sce,
  var_registration = "cellType",
  var_sample_id = "sample",
  gene_ensembl = "ID",
  gene_name = "Symbol"
)

saveRDS(object = rat_modeling_results,
        file   = here("processed-data","spatial_registration","rat_sce_modeling_results.Rds"))


## Perform the spatial registration
xen_modeling_results <- registration_wrapper(
  sce = spe_sub,
  var_registration = "CellTypes",
  var_sample_id = "Sample",
  gene_ensembl = "ID",
  gene_name = "Symbol"
)

saveRDS(object = xen_modeling_results,
        file   = here("processed-data","spatial_registration","xen_modeling_results_CellTypes_v2_rat.Rds"))


cor_res <- layer_stat_cor(
  stats = rat_modeling_results$enrichment, #query
  modeling_results = xen_modeling_results,#reference
  model_type = "enrichment"
)

saveRDS(object = cor_res,
        file   = here("processed-data",
                      "spatial_registration","xenium_rat_snRNA_crossspecies_registration_results_celltypes_v2.Rds"))


pdf(here("plots","spatial_registration","xenium_rat_snRNA_crossspecies_registration.pdf"),height = 8, width = 8)
layer_stat_cor_plot(cor_res)
dev.off()

sessionInfo()
