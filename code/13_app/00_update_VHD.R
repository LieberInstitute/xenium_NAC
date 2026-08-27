# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(here)

here::i_am("code/13_app/00_update_VHD.R")

##### Visium-HD 
message("Loading SFE - ",Sys.time())
sfe_path <- here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
sfe <- loadHDF5SummarizedExperiment(sfe_path)

message("Making object smaller - ",Sys.time())
rowData(sfe) <- rowData(sfe)[,c("ID","Symbol")]
colData(sfe) <- colData(sfe)[,c(1:17,130:137)]
reducedDim(sfe,"NMF_Proj") <- NULL

######## 
message("Adding Img Data - ",Sys.time())
for (sid in unique(sfe$sample_id)) {
  print(sid)
  spatial_dir <- here("processed-data", "01_spaceranger", sid,
                      "outs","segmented_outputs","spatial")
  stopifnot(file.exists(file.path(spatial_dir, "scalefactors_json.json")))
  sf <- jsonlite::fromJSON(file.path(spatial_dir, "scalefactors_json.json"))
  sfe <- addImg(sfe,
    sample_id   = sid,
    image_id    = "lowres",
    imageSource = file.path(spatial_dir, "tissue_lowres_image.png"),
    scale_fct = sf$tissue_lowres_scalef)
 ### Check the alignment. Value should be the full-res pixel range of centroids
 terra::ext(getImg(sfe, sample_id = sid, image_id = "lowres"))
 apply(spatialCoords(sfe)[sfe$sample_id == sid, ], 2, range) 
}

saveHDF5SummarizedExperiment(sfe,here("processed-data","13_app","VisiumHD_sfe"),replace = TRUE)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
