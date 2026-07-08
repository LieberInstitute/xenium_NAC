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
files   <- file.path(proj_dir, paste0("proj_", seq_along(samples), ".rds"))
missing <- files[!file.exists(files)]
if (length(missing))
  stop("Missing projection outputs (rerun just those array indices):\n",
       paste(missing, collapse = "\n"))

#Combine the nmf weights for each cell
proj <- do.call(cbind, lapply(files, readRDS))      # patterns x all cells

#Make sure everything is in the projections and force same order
if (!all(colnames(sfe) %in% colnames(proj)))
  stop("Some SFE cells have no projection — was the split complete?")
proj <- proj[, colnames(sfe), drop = FALSE]           # realign to SPE column order
stopifnot(identical(colnames(sfe),colnames(proj)))

#Add projections to the sfe object
reducedDim(sfe, "NMF_Proj") <- t(proj)             # cells x patterns

#Add raw and scaled NMF weights to column data
for (i in seq_len(n_factors)) {
  col <- paste0("nmf", i)
  w  <- reducedDim(sfe, "NMF_Proj")[, col]
  sfe[[col]] <- w                 
}

#Identify the top NMF weight and assign the cell type
Top_NMF <- max.col(m = reducedDim(sfe,"NMF_Proj"),ties.method = "first")
colData(sfe)$Top_NMF <- paste0("nmf",Top_NMF)

#Load NMF cell types with associated class based on:https://github.com/LieberInstitute/spatial_NAc/blob/145f74fb821473fb36d6845c5f6331ea0024f089/code/19_sLDSC/04-make_plots_human_NAc_NMF.R#L60
NMF_CellTypes <- read.csv(here("processed-data","NMF_CellTypes.csv"))
colData(sfe)$Pred_NMF_CellType <- NMF_CellTypes[match(colData(sfe)$Top_NMF,NMF_CellTypes$NMF_Factor),"Class"]

# ---- Side-by-side raw vs scaled escheR plots --------------------------------
fill_cols <- c("white","lightgrey","red","black")

for (i in seq_len(n_factors)) {
  print(i)
  col    <- paste0("nmf", i)
  
  #Create a directory path for the plots to go into
  out_dir <- here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger","NMF", col)
  dir.create(path = out_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (sample in unique(sfe$sample_id)) {
    print(sample)
    sfe_sub <- sfe[, sfe$sample_id == sample]
    
    p <- make_escheR(sfe_sub) |>
      add_fill(col) +
      scale_fill_gradientn(colours = fill_cols)
    
    ggsave(plot = p,
           filename = file.path(out_dir, paste0(sample, ".png")),
           height = 12, width = 24, dpi = 200)   # width doubled for 2 panels
  }
}

#Calculate log-normalized counts so the object contains everything needed for continued analysis.
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

#Save object (now contains nmf{i} and nmf{i}_scaled columns)
saveHDF5SummarizedExperiment(sfe, here("processed-data", "HD_Full_Analysis", "sfe_with_NMF"),replace = TRUE)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
