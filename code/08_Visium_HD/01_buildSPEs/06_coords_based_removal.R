library(SpatialExperiment) 
library(HDF5Array)
library(ggplot2)
library(escheR)
library(here)

#load filtered spe object
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

#Susbet for sample that includes targetting issue. 
spe_sub <- spe[,spe$sample_id == "H1-XNQ4F2B_A1"]

#Remove cells based on x-y coordinates. 
spe_sub$removal <- ifelse(spe_sub$array_col >= 3000,
                          "remove",
                          "keep")
p <- make_escheR(spe_sub) |>
  add_fill("removal")
ggsave(filename = here("plots","HD_Full_Analysis","cell_QC",
                       "H1-XNQ4F2B_A1_coords_based_removal.png"),
       height = 10, width = 8, dpi = 200)

#Subset the cells
removal_cells <- rownames(colData(spe_sub)[which(spe_sub$removal == "remove"),])
cells_to_keep <- setdiff(colnames(spe), removal_cells)
spe <- spe[, cells_to_keep]

#Now save. 
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- saveHDF5SummarizedExperiment(
  spe, dir = spe_filtered_dir, replace = TRUE, as.sparse = TRUE
)

message(sprintf("Saved filtered SPE to: %s", spe_filtered_dir))

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
