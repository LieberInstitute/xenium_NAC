# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
 
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(here)
 
here::i_am("code/13_app/00_update_VHD.R")
 
##### Visium-HD
message("Loading SFE - ", Sys.time())
sfe_path <- here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
sfe <- loadHDF5SummarizedExperiment(sfe_path)
 
message("Making object smaller - ", Sys.time())
rowData(sfe) <- rowData(sfe)[, c("ID", "Symbol")]
colData(sfe) <- colData(sfe)[, c(1:17, 130:137)]
reducedDim(sfe, "NMF_Proj") <- NULL
 
########
## Convert SFE -> SPE. This drops the sf geometry columns and any terra
## images, leaving a lighter object that serializes cleanly. The
## SpatialFeatureExperiment package provides a coercion method:
message("Converting SFE -> SPE - ", Sys.time())
spe <- as(sfe, "SpatialExperiment")
 
## Confirm spatialCoords carried over; if not, set them from the SFE centroids.
if (is.null(spatialCoords(spe)) || ncol(spatialCoords(spe)) == 0) {
  spatialCoords(spe) <- as.matrix(spatialCoords(sfe))
}
 
## Drop the raw SFE from memory before saving.
rm(sfe); gc()
 
########
## Add lowres images as SpatialImage objects. SPE stores these via
## SpatialExperiment::addImg(), which wraps a raster or a file path -- no
## external pointers, so they survive save/reload.
message("Adding image data - ", Sys.time())
for (sid in unique(spe$sample_id)) {
  message(sid)
  spatial_dir <- here("processed-data", "01_spaceranger", sid,
                      "outs", "segmented_outputs", "spatial")
  sf_path  <- file.path(spatial_dir, "scalefactors_json.json")
  img_path <- file.path(spatial_dir, "tissue_lowres_image.png")
  stopifnot(file.exists(sf_path), file.exists(img_path))
 
  sf <- jsonlite::fromJSON(sf_path)
  spe <- SpatialExperiment::addImg(
    spe,
    sample_id    = sid,
    image_id     = "lowres",
    imageSource  = img_path,
    scaleFactor  = sf$tissue_lowres_scalef,
    load         = TRUE   # read PNG into a LoadedSpatialImage now
  )
 
  ## sanity check on alignment (values should be full-res pixel ranges)
  coord_range <- apply(spatialCoords(spe)[spe$sample_id == sid, ], 2, range)
  message(sid, " centroid range: ",
          paste(round(as.vector(coord_range)), collapse = ", "))
}
 
########
message("Saving SPE - ", Sys.time())
saveHDF5SummarizedExperiment(
  spe,
  here("processed-data", "13_app", "VisiumHD_spe"),
  replace = TRUE
)
 
## Verify round trip in this session.
message("Verifying reload - ", Sys.time())
spe_check <- loadHDF5SummarizedExperiment(
  here("processed-data", "13_app", "VisiumHD_spe")
)
stopifnot(nrow(imgData(spe_check)) == length(unique(spe_check$sample_id)))
message("Reload OK. ", nrow(imgData(spe_check)), " images attached.")
 
## spatialCoords check: must exist, be 2D numeric, match ncol(spe), and
## have a sensible range per sample (non-zero spread, no NAs).
sc <- spatialCoords(spe_check)
stopifnot(
  !is.null(sc),
  is.numeric(sc),
  ncol(sc) == 2,
  nrow(sc) == ncol(spe_check),
  !anyNA(sc)
)
message("spatialCoords OK. ", nrow(sc), " x ", ncol(sc),
        " (", paste(colnames(sc), collapse = ", "), ")")
 
for (sid in unique(spe_check$sample_id)) {
  rng <- apply(sc[spe_check$sample_id == sid, , drop = FALSE], 2, range)
  spread <- rng[2, ] - rng[1, ]
  stopifnot(all(spread > 0))
  message(sid, " coord range: ",
          paste(round(as.vector(rng)), collapse = ", "))
}
 
#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
