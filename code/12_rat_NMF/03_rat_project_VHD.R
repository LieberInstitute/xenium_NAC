# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library("RcppML", lib.loc = "/users/rphillip/R/4.3.x")
library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(orthogene)
library(HDF5Array)
library(Matrix)
library(here)

# ---- CONFIG --------------------------------------------------
sfe_path       <- here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
loadings_path  <- here("processed-data", "rat_NMF", "NMF_Results_k53.Rds")   
sample_col     <- "sample_id"
assay_name     <- "logcounts"

workdir        <- here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF_rat","projector_work")
# ---------------------------------------------------------------------------

proj_dir <- file.path(workdir, "proj_out")
dir.create(proj_dir, recursive = TRUE, showWarnings = FALSE)

###### Get sample from SLURM array task ID
sfe <- loadHDF5SummarizedExperiment(sfe_path)
samples <- unique(sfe[[sample_col]])
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
current_sample <- samples[task_id]
message("Projecting sample: ", current_sample, " (", task_id, " of ", length(samples), ")")

message("Loading rat NMF results - ", Sys.time())
nmf_res <- readRDS(loadings_path)
W <- nmf_res@w

message("Subsetting the sfe object to this sample - ", Sys.time())
sub <- sfe[, sfe[[sample_col]] == current_sample]

gene_map <- convert_orthologs(
  gene_df         = rownames(sub),
  gene_input      = "listA",
  gene_output     = "dict",              
  input_species   = "human",
  output_species  = "rat",
  non121_strategy = "drop_both_species",
  method          = "gprofiler"
)


message("Subsetting and renaming to rat orthologs - ", Sys.time())
keep    <- rownames(sub) %in% names(gene_map)
sub_rat <- sub[keep, ]
rownames(sub_rat) <- unlist(gene_map[rownames(sub_rat)])

#Subset the object and the loading matrix to shared genes. 
common_genes <- intersect(rownames(W), rownames(sub_rat))


W        <- W[common_genes, , drop = FALSE]
sub_rat  <- sub_rat[common_genes, , drop = FALSE]

message("Extracting logcounts - ", Sys.time())
expr <- as(assay(sub_rat, assay_name), "CsparseMatrix")

# patterns x cells
message("Projecting with RcppML::project - ", Sys.time())
proj <- project(expr, w = W, L1 = 0)

# remove rowSums == 0
proj1 <- proj[rowSums(proj) == 0, , drop = FALSE]
# keep the rowsums==0 in separate object
proj2 <- proj[rowSums(proj) != 0, , drop = FALSE]

# normalize
proj2 <- apply(proj2, 1, function(x){x / sum(x)})
proj1 <- t(proj1)

# combine and force into same order
proj_final <- cbind(proj2, proj1)
proj_final <- proj_final[, match(rownames(proj), colnames(proj_final))]

# Save raw and normalized
saveRDS(proj,       file.path(proj_dir, paste0("proj_", task_id, ".rds")))
saveRDS(proj_final, file.path(proj_dir, paste0("proj_final_", task_id, ".rds")))

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
