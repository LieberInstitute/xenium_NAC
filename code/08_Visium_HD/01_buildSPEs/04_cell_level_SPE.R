#Goal: Build a cell-level SPE using output from bin2cell 
#Code based on script from Nick Eagles: https://github.com/LieberInstitute/lc_visium_hd/blob/a13f88b6b66bc0fd643a8d690694c454905b3eae/code/03_build_spe/03_build_cell_spe.R
library(SpatialExperiment)
library(zellkonverter)
library(BiocParallel)
library(sessioninfo)
library(tidyverse)
library(HDF5Array)
library(escheR)
library(scran)
library(here)

sample_ids <- c("H1-XKYDCP3_A1","H1-XKYDCP3_D1","H1-M3TCP9V_A1","H1-M3TCP9V_D1")
adata_in_paths <- here("processed-data","HD_Full_Analysis","bin2cell",sprintf('%s.h5ad', sample_ids))

spe_bin_dir <- here("processed-data","HD_Full_Analysis","SPEs","spe_raw.Rds")

spe_norm_dir <- here("processed-data","HD_Full_Analysis","SPEs","spe_cell_norm")
spe_raw_dir  <- here("processed-data","HD_Full_Analysis","SPEs","spe_cell_raw")

#   Given one sample ID, the path to a single-sample AnnData from bin2cell, and
#   a potentially multi-sample bin-level SpatialExperiment, return a
#   single-sample SpatialExperiment formed from the AnnData
anndata_to_spe <- function(sample_id, ad_in_path, spe_bin) {
  message(
    Sys.time(),
    sprintf(
      " - Forming basic SPE from the bin2cell-output AnnData for sample %s...",
      sample_id
    )
  )
  
  #   Read AnnData into a SingleCellExperiment via zellkonverter
  sce <- readH5AD(ad_in_path, use_hdf5 = TRUE)
  
  #   readH5AD stores the main matrix as assay "X" by default.
  #   Rename to "counts" so that counts(sce) works downstream.
  if (!"counts" %in% assayNames(sce) && "X" %in% assayNames(sce)) {
    assayNames(sce)[assayNames(sce) == "X"] <- "counts"
  }
  
  #   Extract spatial coordinates from obsm. bin2cell stores coordinates in
  #   obsm['spatial'] (or similar). readH5AD puts obsm entries in reducedDims.
  spatial_keys <- grep(
    "spatial", reducedDimNames(sce), value = TRUE, ignore.case = TRUE
  )
  if (length(spatial_keys) > 0) {
    coords <- reducedDim(sce, spatial_keys[1])
  } else {
    #   Fallback: check for coordinate columns in colData
    coord_cols <- grep(
      "spatial|array_row|array_col|pxl_", colnames(colData(sce)),
      value = TRUE, ignore.case = TRUE
    )
    stopifnot(
      "Could not find spatial coordinates in AnnData" = length(coord_cols) >= 2
    )
    coords <- as.matrix(colData(sce)[, coord_cols[1:2]])
  }
  
  #   Ensure coords are a numeric matrix with proper column names
  coords <- as.matrix(coords[, 1:2])
  colnames(coords) <- c("pxl_col_in_fullres", "pxl_row_in_fullres")
  
  #   Build the SpatialExperiment from the SCE
  spe <- SpatialExperiment(
    assays = list(counts = counts(sce)),
    colData = colData(sce),
    rowData = rowData(sce),
    spatialCoords = coords
  )
  
  #   Set rownames to match the SCE (Ensembl IDs from bin2cell)
  rownames(spe) <- rownames(sce)
  
  #   Fix sample ID and key
  spe$sample_id <- sample_id
  spe$key <- paste(colnames(spe), sample_id, sep = "_")
  colnames(spe) <- spe$key
  
  #   Transfer over image-related data from the bin-level object (so far, it
  #   doesn't seem there is a simpler method for retaining this data when
  #   converting from AnnData). Also transfer rowData. For unclear reasons,
  #   some genes are dropped during creation of the bin-level object, and so
  #   this step may drop genes in the cell-level object where rowData doesn't
  #   exist
  imgData(spe) <- imgData(spe_bin)
  
  stopifnot(all(rownames(spe_bin) %in% rownames(spe)))
  num_dropped <- length(setdiff(rownames(spe), rownames(spe_bin)))
  if (num_dropped > 0) {
    warning(
      sprintf(
        "Dropping %s of %s genes for %s not present in 'spe_bin'...",
        num_dropped, nrow(spe), sample_id
      )
    )
  }
  spe <- spe[intersect(rownames(spe), rownames(spe_bin))]
  rowData(spe) <- rowData(spe_bin[rownames(spe),])
  
  message(
    Sys.time(),
    sprintf(" - Adding spatialLIBD metrics for sample %s...", sample_id)
  )
  
  is_mito <- as.logical(seqnames(spe_bin[rownames(spe),]) == "chrM")

  spe <- addPerCellQCMetrics(spe, subsets = list(mito = is_mito))

  #   spatialLIBD-specific columns
  rowData(spe)$gene_search <- paste0(
    rowData(spe)$gene_name, "; ", rowData(spe)$gene_id
  )
  spe$ManualAnnotation <- "NA"
  return(spe)
}

################################################################################
#   Main
################################################################################

#-------------------------------------------------------------------------------
#   Build and save raw SpatialExperiment
#-------------------------------------------------------------------------------

spe_bin <- readRDS(spe_bin_dir)

#   Individually build single-sample SPEs from the individual AnnDatas, then
#   merge
spe_list <- list()
for (i in seq_len(length(sample_ids))) {
  spe_list[[sample_ids[[i]]]] <- anndata_to_spe(
    sample_ids[[i]], adata_in_paths[[i]],
    spe_bin[, spe_bin$sample_id == sample_ids[[i]]]
  )
}

message(Sys.time(), " - Merging SPEs across samples")
gene_sets <- unname(lapply(spe_list, rownames))
stopifnot(all(sapply(gene_sets[-1], identical, gene_sets[[1]])))
spe <- do.call(cbind, spe_list)

#-------------------------------------------------------------------------------
#   Sanity check: plot spatial coordinates with escheR to verify orientation
#-------------------------------------------------------------------------------

message(Sys.time(), " - Plotting spatial coordinate sanity checks with escheR")

plot_dir <- here("plots", "HD_Full_Analysis", "cell_spe_build")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

for (sid in sample_ids) {
    spe_sub <- spe[, spe$sample_id == sid]

    #   Trim colData for plotting speed
    colData(spe_sub) <- colData(spe_sub)[, c("sample_id", "sum")]

    p <- make_escheR(spe_sub, spot_size = 0.5) |>
        add_fill(var = "sum", point_size = 0.5) +
        scale_fill_gradientn(
            colors = viridisLite::plasma(256),
            name = "Total UMI"
        ) +
        ggtitle(paste(sid, "- Spatial coordinate check")) +
        theme_void() +
        theme(
            plot.title = element_text(hjust = 0.5, size = 14),
            legend.position = "right"
        )

    ggsave(
        file.path(plot_dir, sprintf("%s_coord_check.png", sid)), p,
        width = 10, height = 8, dpi = 200
    )
}

message(
    "  >> Inspect plots in ", plot_dir, "\n",
    "  >> If tissue appears transposed or flipped, swap or negate the\n",
    "  >> coordinate columns in the anndata_to_spe() function."
)

#   Save now to allow assays to become HDF5-backed in hopes of driving memory
#   down during log-normalization
message(Sys.time(), " - Saving raw SPE")
spe <- saveHDF5SummarizedExperiment(
  spe, dir = spe_raw_dir, replace = TRUE, as.sparse = TRUE
)

#-------------------------------------------------------------------------------
#   Build and save log-normalized and filtered SpatialExperiment
#-------------------------------------------------------------------------------

#   Filter SPE: drop cells with 0 counts for all genes, and drop genes with 0
#   counts in every cell
message(Sys.time(), " - Filtering genes and spots")
spe <- spe[rowSums(assays(spe)$counts) > 0, colSums(assays(spe)$counts) > 0]

#   Use library-size normalization (normalization by deconvolution is not
#   computationally feasible with data this large)
message(Sys.time(), ' - Performing log normalization...')
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

#   Save normalized object
message(Sys.time(), " - Saving normalized SPE")
spe <- saveHDF5SummarizedExperiment(
  spe, dir = spe_norm_dir, replace = TRUE, as.sparse = TRUE
)

message("Memory usage:")
gc()
session_info()
