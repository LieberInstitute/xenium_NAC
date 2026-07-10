library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(patchwork) 
library(HDF5Array)
library(ggplot2)
library(scuttle)
library(escheR)
library(scales) 
library(here)

# ---- CONFIG --------------------------------------------------
sfe_path      <- here("processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2")
loadings_path <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/HD_Full_Analysis/LabelTransfer/NMF/snRNA_nmf_results.rds"       # genes x patterns matrix (W), gene rownames
sample_col    <- "sample_id"
assay_name    <- "logcounts"
workdir       <- here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF","projector_work")
n_factors     <- 66                            # number of NMF patterns
# ---------------------------------------------------------------------------
manifest_path <- file.path(workdir, "sample_manifest.txt")
split_dir     <- file.path(workdir, "per_sample")
proj_dir      <- file.path(workdir, "proj_out")

###### Get sample from SLURM array task ID
sfe <- loadHDF5SummarizedExperiment(sfe_path)
samples <- readLines(manifest_path)

raw_files  <- file.path(proj_dir, paste0("proj_",       seq_along(samples), ".rds"))
norm_files <- file.path(proj_dir, paste0("proj_final_", seq_along(samples), ".rds"))



#Combine the raw and normalized nmf weights across all samples
proj      <- do.call(cbind, lapply(raw_files,  readRDS))   
proj_norm <- do.call(rbind, lapply(norm_files, readRDS))  

stopifnot(identical(rownames(proj), colnames(proj_norm)))
stopifnot(identical(colnames(proj), rownames(proj_norm)))   

proj      <- proj[,      colnames(sfe), drop = FALSE] 
proj_norm <- proj_norm[colnames(sfe), , drop = FALSE] 
stopifnot(identical(colnames(sfe), colnames(proj)))
stopifnot(identical(colnames(sfe), rownames(proj_norm)))

reducedDim(sfe, "NMF_Proj")      <- t(proj)       
reducedDim(sfe, "NMF_Proj_norm") <- proj_norm      

#Add raw and normalized NMF weights to column data
for (i in seq_len(n_factors)) {
  print(i)
  col      <- paste0("nmf", i)
  col_norm <- paste0("nmf", i, "_norm")
  sfe[[col]]      <- reducedDim(sfe, "NMF_Proj")[, col]
  sfe[[col_norm]] <- reducedDim(sfe, "NMF_Proj_norm")[, col]
}



# ---- Side-by-side raw vs normalized escheR plots ----------------------------
fill_cols <- c("white","lightgrey","red","black")

for (i in seq_len(n_factors)) {
  print(i)
  col      <- paste0("nmf", i)
  col_norm <- paste0("nmf", i, "_norm")
  
  #Create a directory path for the plots to go into
  out_dir <- here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF", col)
  dir.create(path = out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sample in unique(sfe$sample_id)) {
    print(sample)
    sfe_sub <- sfe[, sfe$sample_id == sample]
    
    p_raw <- make_escheR(sfe_sub) |>
      add_fill(col) +
      scale_fill_gradientn(colours = fill_cols)
    
    ggsave(plot = p_raw,
           filename = file.path(out_dir, paste0("raw_",sample, ".png")),
           height = 12, width = 12, dpi = 200) 
    
    p_norm <- make_escheR(sfe_sub) |>
      add_fill(col_norm) +
      scale_fill_gradientn(colours = fill_cols) +
      ggtitle("normalized")
    
    ggsave(plot = p_norm,
           filename = file.path(out_dir, paste0("normalized_",sample, ".png")),
           height = 12, width = 12, dpi = 200) 
  }
}

#Calculate log-normalized counts so the object contains everything needed for continued analysis.
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

#Save object (now contains nmf{i} and nmf{i}_norm columns, plus NMF_Proj / NMF_Proj_norm reducedDims)
saveHDF5SummarizedExperiment(sfe, here("processed-data", "HD_Full_Analysis", "sfe_with_NMF"),replace = TRUE)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
