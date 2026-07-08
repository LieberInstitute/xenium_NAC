#!/usr/bin/env Rscript
# ============================================================================
# Project single-cell NMF patterns into Visium HD data, one sample per task.
#
# Three modes (normally driven by run_projectr.sh):
#   Rscript project_nmf_array.R presplit            # run once, splits the SPE
#   Rscript project_nmf_array.R project <index>     # one array task per sample
#   Rscript project_nmf_array.R recombine           # run once, merges results
# ============================================================================

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(HDF5Array)
  library(here)
  library(SpatialExperiment)
  library(projectR)
})

# ---- CONFIG (edit these) ---------------------------------------------------
spe_path      <- here("processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2") # full SpatialExperiment (all samples)
loadings_path <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/LabelTransfer/NMF/snRNA_nmf_results.rds"       # genes x patterns matrix (W), gene rownames
sample_col    <- "sample_id"                   # colData column with the sample label
assay_name    <- "logcounts"                   # MUST be globally/consistently normalized
out_spe_path  <- here("processed-data", "HD_Full_Analysis","SPEs", "spe_with_NMFproj") # where recombine writes the result
reduceddim_nm <- "NMFproj"                      # reducedDim slot name for the projection

workdir       <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/LabelTransfer/NMF/projectr_work"               # scratch for intermediate files
# ---------------------------------------------------------------------------

manifest_path <- file.path(workdir, "sample_manifest.txt")
split_dir     <- file.path(workdir, "per_sample")
proj_dir      <- file.path(workdir, "proj_out")

args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
if (is.na(mode) || !mode %in% c("presplit", "project", "recombine"))
  stop("First argument must be one of: presplit | project | recombine")

tag <- function(i) sprintf("%03d", as.integer(i))   # zero-padded, collision-proof file ids

# ===========================================================================
if (mode == "presplit") {

  dir.create(split_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(proj_dir,  recursive = TRUE, showWarnings = FALSE)

  spe <- loadHDF5SummarizedExperiment(spe_path)
  rownames(spe) <- rowData(spe)$gene_id

  if (!assay_name %in% assayNames(spe))
    stop(sprintf("Assay '%s' not found. Available: %s",
                 assay_name, paste(assayNames(spe), collapse = ", ")))
  if (!sample_col %in% colnames(colData(spe)))
    stop(sprintf("Column '%s' not found in colData.", sample_col))
  if (anyDuplicated(colnames(spe)))
    stop("Cell names (colnames) are not unique; recombine relies on unique names.")

  samples <- sort(unique(as.character(colData(spe)[[sample_col]])))
  writeLines(samples, manifest_path)
  message(sprintf("Found %d samples. Manifest -> %s", length(samples), manifest_path))

  for (i in seq_along(samples)) {
    s   <- samples[i]
    sub <- spe[, as.character(colData(spe)[[sample_col]]) == s]
    saveHDF5SummarizedExperiment(sub, file.path(split_dir, paste0("sample_", tag(i))))
    message(sprintf("  [%s] %-20s %d cells", tag(i), s, ncol(sub)))
  }
  message("presplit done. Submit the array as --array=1-", length(samples))

# ===========================================================================
} else if (mode == "project") {

  i <- as.integer(args[2])
  if (is.na(i)) stop("project mode needs an array index, e.g. 'project 3'")

  samples <- readLines(manifest_path)
  s       <- samples[i]
  message(sprintf("[%s] projecting sample: %s", tag(i), s))

  nmf_res <- readRDS(loadings_path)
  W <- nmf_res@w
  if (is.null(rownames(W))) stop("Loadings matrix must have gene rownames.")

  sub  <- readRDS(file.path(split_dir, paste0("sample_", tag(i), ".rds")))
  expr <- as.matrix(assay(sub, assay_name))         


  # patterns x cells. Use full = TRUE if you also want Wald p-values (returns a list).
  proj <- projectR(data     = expr,
                   loadings = W)

  saveRDS(proj, file.path(proj_dir, paste0("proj_", tag(i), ".rds")))
  message(sprintf("[%s] done: %d patterns x %d cells", tag(i), nrow(proj), ncol(proj)))

# ===========================================================================
} else if (mode == "recombine") {

  samples <- readLines(manifest_path)
  files   <- file.path(proj_dir, paste0("proj_", tag(seq_along(samples)), ".rds"))

  missing <- files[!file.exists(files)]
  if (length(missing))
    stop("Missing projection outputs (rerun just those array indices):\n",
         paste(missing, collapse = "\n"))

  proj <- do.call(cbind, lapply(files, readRDS))       # patterns x all cells

  spe <- loadHDF5SummarizedExperiment(spe_path)
  if (!all(colnames(spe) %in% colnames(proj)))
    stop("Some SPE cells have no projection — was the split complete?")
  proj <- proj[, colnames(spe), drop = FALSE]           # realign to SPE column order

  reducedDim(spe, reduceddim_nm) <- t(proj)             # cells x patterns
  saveHDF5SummarizedExperiment(spe, out_spe_path)
  message(sprintf("recombine done: %s  (reducedDim '%s', %d cells x %d patterns)",
                  out_spe_path, reduceddim_nm, ncol(proj), nrow(proj)))
}
