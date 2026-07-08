library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(projectR)
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

###### Get sample from SLURM array task ID
sfe <- loadHDF5SummarizedExperiment(sfe_path)
samples <- unique(sfe$sample_id)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message("Projecting sample: ", current_sample, " (", task_id, " of ", length(samples), ")")

message("Loading NMF results from snRNA-seq - ", Sys.time())
nmf_res <- readRDS(loadings_path)
W <- nmf_res@w
if (is.null(rownames(W))) stop("Loadings matrix must have gene rownames.")

message("Subsetting the sfe object - ", Sys.time())
sub  <- loadHDF5SummarizedExperiment(file.path(split_dir,paste0("sample_00", task_id)))
rownames(sub) <- rowData(sub)$ID

message("Extracting logcounts - ", Sys.time())
expr <- as.matrix(assay(sub, assay_name))         


# patterns x cells
message("Projecting - ", Sys.time())
proj <- projectR(data     = expr,
                 loadings = W)

saveRDS(proj, file.path(proj_dir, paste0("proj_", task_id, ".rds")))
message(sprintf("[%s] done: %d patterns x %d cells", task_id, nrow(proj), ncol(proj)))


#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

