#Goal: Use RCTD to transfer cell types from snRNA to VisiumHD (per-sample)
library(SpatialFeatureExperiment)
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
ref_counts <- as(assay(sce, "counts"), "dgCMatrix")


reference_se <- SummarizedExperiment(
  assays = list(counts = ref_counts),
  colData = colData(sce)$CellType.Final
)

###### Prep spatial data
#Read in the filtered sfe object
#Now save.
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

###### Get sample from SLURM array task ID
samples <- unique(sfe$sample_id)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message("Processing sample: ", current_sample, " (", task_id, " of ", length(samples), ")")

###### Subset to current sample
sfe_sub <- sfe[, sfe$sample_id == current_sample]
message("Number of spots for this sample: ", ncol(sfe_sub))
# to avoid the sparse->dense coercion
counts_mat <- as(assay(sfe_sub, "counts"), "dgCMatrix")

spatial_sfe <- SpatialExperiment(
  assays = list(counts = counts_mat),
  spatialCoords = cbind(
    x = spatialCoords(sfe_sub)[,1],
    y = spatialCoords(sfe_sub)[,2]
  )
)

############# Run RCTD
message("Creating Rctd object - ", Sys.time())
rctd_data <- createRctd(spatial_sfe, reference_se, cell_type_col = "X",UMI_min = 20)

message("Running Rctd - ", Sys.time())
results <- runRctd(rctd_data,rctd_mode = "doublet",max_cores=4)

outdir <- here("processed-data", "HD_Full_Analysis", "LabelTransfer","SpaceRanger","RCTD")
saveRDS(results, file.path(outdir, paste0("doublet_snRNA_VisiumHD_RCTD_", current_sample, ".Rds")))

message("Done with sample: ", current_sample, " - ", Sys.time())

### Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
