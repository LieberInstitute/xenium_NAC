# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(here)

#Load the object
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

sfe


#Add spatial coordinates to the object quickly to help identify certain cells
sfe$X_coords <- spatialCoords(sfe)[,1]
sfe$Y_coords <- spatialCoords(sfe)[,2]


#Subset for sample that needs coord-based removal 
samp <-  "H1-XNQ4F2B_A1"
sfe_sub <- sfe[,sfe$sample_id == samp]

#Cells with x_coord greater than 130000 need to be removed
sfe_sub$removal <- ifelse(sfe_sub$X_coords >= 130000,
                          "remove",
                          "keep")

table(sfe_sub$removal)

coords_removal <- ggplot(data = colData(sfe_sub),aes(x = X_coords, y = Y_coords,color = removal)) +
  geom_point()

ggsave(plot = coords_removal,
       filename = here("plots","HD_Full_Analysis","cell_QC",
                       "H1-XNQ4F2B_A1_coords_based_removal.png"),
       height = 12, width = 10, dpi = 200)


#Subset the cells
removal_cells <- rownames(colData(sfe_sub)[which(sfe_sub$removal == "remove"),])
cells_to_keep <- setdiff(colnames(sfe), removal_cells)
sfe <- sfe[, cells_to_keep]

sfe

#Now save. 
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)

saveHDF5SummarizedExperiment(sfe, dir = sfe_filtered_dir, replace = TRUE, as.sparse = TRUE)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
