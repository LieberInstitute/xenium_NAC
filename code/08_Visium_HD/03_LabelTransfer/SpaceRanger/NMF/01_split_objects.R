library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(scuttle)
library(here)

# ---- CONFIG --------------------------------------------------
sfe_path      <- here("processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2") 
loadings_path <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/LabelTransfer/NMF/snRNA_nmf_results.rds"       # genes x patterns matrix (W), gene rownames
sample_col    <- "sample_id"                   
assay_name    <- "logcounts"                   
workdir       <- here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF","projector_work")
# ---------------------------------------------------------------------------

manifest_path <- file.path(workdir, "sample_manifest.txt")
split_dir     <- file.path(workdir, "per_sample")
proj_dir      <- file.path(workdir, "proj_out")

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]

tag <- function(i) sprintf("%03d", as.integer(i))   # zero-padded, collision-proof file ids

dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(proj_dir,  recursive = TRUE, showWarnings = FALSE)

#Load the object
sfe <- loadHDF5SummarizedExperiment(sfe_path)
#nmf loadings have ensembl gene ids so make the sfe have them too
rownames(sfe) <- rowData(sfe)$gene_id

#calculate log-normalized counts 
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

if (!assay_name %in% assayNames(sfe))
  stop(sprintf("Assay '%s' not found. Available: %s",
               assay_name, paste(assayNames(sfe), collapse = ", ")))
if (!sample_col %in% colnames(colData(sfe)))
  stop(sprintf("Column '%s' not found in colData.", sample_col))
if (anyDuplicated(colnames(sfe)))
  stop("Cell names (colnames) are not unique; recombine relies on unique names.")

samples <- sort(unique(as.character(colData(sfe)[[sample_col]])))
writeLines(samples, manifest_path)
message(sprintf("Found %d samples. Manifest -> %s", length(samples), manifest_path))

for (i in seq_along(samples)) {
  s   <- samples[i]
  sub <- sfe[, as.character(colData(sfe)[[sample_col]]) == s]
  saveHDF5SummarizedExperiment(sub, file.path(split_dir, paste0("sample_", tag(i))))
  message(sprintf("  [%s] %-20s %d cells", tag(i), s, ncol(sub)))
}

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
