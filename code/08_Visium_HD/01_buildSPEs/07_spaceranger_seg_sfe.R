library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(SpatialExperiment)
library(DropletUtils)
library(here)
library(sf)

## ---- pick the sample (SLURM array) -----------------------------------------
sample_info <- read.csv(
  here("processed-data", "visiumHD_sample_info_NAc.csv"),
  stringsAsFactors = FALSE
)
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
sample_id <- sample_info$sample_id[task_id]
message(format(Sys.time()), " | sample: ", sample_id)

out_dir <- here("processed-data", "HD_Full_Analysis", "spaceranger_sfe")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

seg_dir <- here("processed-data", "01_spaceranger", sample_id,
                "outs", "segmented_outputs")
h5 <- file.path(seg_dir, "filtered_feature_cell_matrix.h5")
gj <- file.path(seg_dir, "cell_segmentations.geojson")
stopifnot(file.exists(h5), file.exists(gj))

## ---- counts (cell x gene) --------------------------------------------------
sce <- read10xCounts(h5, col.names = TRUE)
#   read10xCounts returns an HDF5-backed DelayedArray for .h5 input; pull the
#   counts into memory (sparse) so the saved .rds is self-contained / portable.
counts(sce) <- as(counts(sce), "dgCMatrix")
message(format(Sys.time()), " | counts: ",
        nrow(sce), " genes x ", ncol(sce), " cells")
message("  example cell name: ", colnames(sce)[1])

#   numeric cell id from the matrix colnames: take the first run of digits, so
#   both a bare "12345" and a prefixed "cellid_12345-1" resolve to 12345
cell_num <- suppressWarnings(
  as.numeric(sub("^\\D*([0-9]+).*$", "\\1", colnames(sce)))
)
if (all(is.na(cell_num))) {
  stop("could not parse a numeric cell id from colnames; example: ",
       colnames(sce)[1])
}

## ---- polygons --------------------------------------------------------------
seg <- st_read(gj, quiet = TRUE)

#   The coordinates are microscope-image pixels, but the geojson declares a
#   WGS84 (lon/lat) CRS. Drop it so sf treats the polygons as planar -- otherwise
#   st_centroid() below assumes spherical geometry and returns wrong centroids.
st_crs(seg) <- NA

message(format(Sys.time()), " | ", nrow(seg), " cell polygons")

#   The geojson has exactly two columns: cell_id (bare integers) and geometry.
stopifnot("cell_id" %in% names(seg))
seg_id <- as.numeric(seg$cell_id)

## ---- align cells <-> polygons ---------------------------------------------
idx <- match(cell_num, seg_id)
keep <- !is.na(idx)
if (any(!keep)) message("  dropping ", sum(!keep), " cells with no polygon")
sce <- sce[, keep]
seg2 <- seg[idx[keep], ]

#   per-cell coordinate = polygon centroid (microscope-image pixels)
suppressWarnings(cent <- st_coordinates(st_centroid(st_geometry(seg2))))

## ---- assemble SpatialExperiment -------------------------------------------
spe <- SpatialExperiment(
  assays = list(counts = counts(sce)),
  rowData = rowData(sce),
  colData = colData(sce),
  spatialCoords = as.matrix(cent)
)
rownames(spe) <- rowData(sce)$ID            # Ensembl IDs (match across samples)
colnames(spe) <- colnames(sce)

#   tag sample + globally-unique cell names (for the merge)
spe$sample_id <- sample_id
colnames(spe) <- paste(sample_id, colnames(spe), sep = "_")

## ---- to SFE + attach cell polygons ----------------------------------------
rownames(seg2) <- colnames(spe)
sfe <- toSpatialFeatureExperiment(spe)
colGeometries(sfe) <- list(cellseg = seg2)
message(format(Sys.time()), " | SFE: ", ncol(sfe), " cells; colGeometries: ",
        paste(colGeometryNames(sfe), collapse = ", "))

## ---- save ------------------------------------------------------------------
out_path <- file.path(out_dir, paste0(sample_id, "_sfe.rds"))
saveRDS(sfe, out_path)
message(format(Sys.time()), " | wrote ", out_path)

sessionInfo()
