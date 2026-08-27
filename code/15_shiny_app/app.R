library(spatialLIBD)
library(markdown)
library(here)
library(qs2)

#   For interactive testing at JHPCE
# setwd(here('code', '12_apps_and_sharing', 'shiny_HD_app'))

posit_connect_file <- "/r_data/lcollado/Posit_Connect_shiny_apps/nac/NAC_AP/VisiumHD_sfe"
if (file.exists(posit_connect_file)) {
    ## Location for the https://conn1.libd.org/ server
    spe <- HDF5Array::loadHDF5SummarizedExperiment(posit_connect_file)
} else {
    spe <- HDF5Array::loadHDF5SummarizedExperiment(here::here("nac", "MD_thalamus", "processed-data", "final_spe_clustered"))
}

discrete_vars = c(
    'ManualAnnotation', 'ficture_cluster', 'banksy_cluster', 'cell_type',
    'cell_type_colors'
)
continuous_vars = c(
    'bin_count', 'sum_umi', 'sum_gene', 'expr_chrM', 'expr_chrM_ratio'
)
all_vars_pb = c(
    'cell_type', 'cell_type_colors', 'sum_umi', 'sum_gene', 'expr_chrM',
    'expr_chrM_ratio', 'ncells', 'sample_id', 'donor'
)

## spatialLIBD uses golem
options("golem.app.prod" = TRUE)

## You need this to enable shinyapps to install Bioconductor packages
options(repos = BiocManager::repositories())

#   Load objects
spe = qs_read('spe_shiny.qs2')
modeling_results = qs_read('modeling_results.qs2')
sce_pb = qs_read('sce_pb_shiny.qs2')
sig_genes = qs_read('sig_genes_shiny.qs2')

stopifnot(all(all_vars_pb %in% colnames(colData(sce_pb))))
colData(sce_pb) = colData(sce_pb)[, all_vars_pb]

## Deploy the website
run_app(
    spe,
    sce_layer = sce_pb,
    modeling_results = modeling_results,
    sig_genes = sig_genes,
    title = "habenula_atlas_Visium_HD",
    spe_discrete_vars = discrete_vars,
    spe_continuous_vars = continuous_vars,
    default_cluster = "cell_type",
    docs_path = 'www',
    is_stitched = TRUE
)
