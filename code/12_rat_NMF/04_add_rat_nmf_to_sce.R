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
sfe_path      <- here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
loadings_path <- here("processed-data", "rat_NMF", "NMF_Results_k53.Rds")
assay_name    <- "logcounts"
sample_col    <- "sample_id"
workdir       <- here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF_rat","projector_work")
n_factors     <- 53                            
# ---------------------------------------------------------------------------
proj_dir <- file.path(workdir, "proj_out")

###### Load the same sfe object and derive the same sample order used in 02
sfe <- loadHDF5SummarizedExperiment(sfe_path)
samples <- unique(sfe[[sample_col]])   

raw_files  <- file.path(proj_dir, paste0("proj_",       seq_along(samples), ".rds"))
norm_files <- file.path(proj_dir, paste0("proj_final_", seq_along(samples), ".rds"))


proj      <- do.call(cbind, lapply(raw_files,  readRDS))   
proj_norm <- do.call(rbind, lapply(norm_files, readRDS))  

stopifnot(identical(rownames(proj), colnames(proj_norm)))   # same patterns, same order
stopifnot(identical(colnames(proj), rownames(proj_norm)))   # same cells, same order

proj      <- proj[,      colnames(sfe), drop = FALSE]   # patterns x cells
proj_norm <- proj_norm[colnames(sfe), , drop = FALSE]   # cells x patterns
stopifnot(identical(colnames(sfe), colnames(proj)))
stopifnot(identical(colnames(sfe), rownames(proj_norm)))

#Add both raw and normalized projections to the sfe object as separate reducedDims
reducedDim(sfe, "NMF_Proj_rat")      <- t(proj)       
reducedDim(sfe, "NMF_Proj_rat_norm") <- proj_norm      

for (i in seq_len(n_factors)) {
  col      <- paste0("nmf", i, "_rat")
  col_norm <- paste0("nmf", i, "_rat_norm")
  sfe[[col]]      <- reducedDim(sfe, "NMF_Proj_rat")[, colnames(proj_norm)[i]]
  sfe[[col_norm]] <- reducedDim(sfe, "NMF_Proj_rat_norm")[, colnames(proj_norm)[i]]
}

# ---- Side-by-side raw vs normalized escheR plots ----------------------------
fill_cols <- c("lightgrey","red","black")

for (i in seq_len(n_factors)) {
  print(i)
  col      <- paste0("nmf", i, "_rat")
  col_norm <- paste0("nmf", i, "_rat_norm")
  
  #Create a directory path for the plots to go into
  out_dir <- here("plots","rat_NMF", col)
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
      scale_fill_gradientn(colours = fill_cols) 
    
    ggsave(plot = p_norm,
           filename = file.path(out_dir, paste0("normalizezd_",sample, ".png")),
           height = 12, width = 12, dpi = 200)  
  }
}

#Calculate log-normalized counts so the object contains everything needed for continued analysis.
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

#Save object (now contains nmf{i}_rat / nmf{i}_rat_norm columns, plus NMF_Proj_rat[_norm] reducedDims)
saveHDF5SummarizedExperiment(sfe, here("processed-data", "HD_Full_Analysis", "sfe_with_rat_NMF"), replace = TRUE)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
