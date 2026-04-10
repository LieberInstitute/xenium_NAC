#Goal: Use RCTD to transfer cell types from snRNA to VisiumHD (per-sample)
library(SingleCellExperiment)
library(SpatialExperiment)
library(zellkonverter)
library(HDF5Array)
library(spacexr)
library(Matrix)
library(here)

###### Prep snRNA-seq data (same for all samples)
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Remove the neuronal ambiguous population from further analysis. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]
#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

sce

rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq

reference_se <- SummarizedExperiment(
    assays = list(counts = assay(sce, "counts")),
    colData = colData(sce)$CellType.Final
)

###### Prep spatial data
spe_filtered_path <- here(
    "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)
spe <- loadHDF5SummarizedExperiment(spe_filtered_path)

###### Get sample from SLURM array task ID
samples <- unique(spe$sample_id)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message("Processing sample: ", current_sample, " (", task_id, " of ", length(samples), ")")

###### Subset to current sample
spe_sub <- spe[, spe$sample_id == current_sample]
message("Number of spots for this sample: ", ncol(spe_sub))

spatial_spe <- SpatialExperiment(
    assays = list(counts = assay(spe_sub, "counts")),
    spatialCoords = cbind(
        x = spe_sub$pxl_col_in_fullres,
        y = spe_sub$pxl_row_in_fullres
    )
)

############# Run RCTD
message("Creating Rctd object - ", Sys.time())
rctd_data <- createRctd(spatial_spe, reference_se, cell_type_col = "X",UMI_min = 1,pixel_count_min = 1, UMI_min_sigma = 1) #Turn off any filtering RCTD does

message("Running Rctd - ", Sys.time())
results <- runRctd(rctd_data,rctd_mode = "full",max_cores=15)

outdir <- here("processed-data", "HD_Full_Analysis", "LabelTransfer","RCTD")
saveRDS(results, file.path(outdir, paste0("snRNA_VisiumHD_RCTD_", current_sample, ".Rds")))

message("Done with sample: ", current_sample, " - ", Sys.time())

### Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
