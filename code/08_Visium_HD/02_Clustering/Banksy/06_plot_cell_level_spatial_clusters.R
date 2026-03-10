#Goal: Plot cell level spatial clusters
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(escheR)
library(here)


################################################################################
#   Set up
################################################################################
spe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs",
  "spe_cell_bansky_spatial_clustered"
))

#   Identify cluster columns from BANKSY
clust_cols <- grep("^clust_", colnames(colData(spe)), value = TRUE)

samples <- unique(spe$sample_id)

################################################################################
#   Generate colors
################################################################################
message(Sys.time(), " | Generating Colors")

cluster_colors <- vector(mode = "list", length = 3)
names(cluster_colors) <- clust_cols
### Generate cluster colors for each clust_fol
for(clust_col in clust_cols){
  #Make colors for CellType
  cluster_colors[[clust_col]] <- Polychrome::createPalette(length(unique(spe[[clust_col]])),
                                                           c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
  names(cluster_colors[[clust_col]]) <- unique(spe[[clust_col]])
  saveRDS(object = cluster_colors[[clust_col]],
          file = here("processed-data","HD_Full_Analysis","Cluster_colors",paste(clust_col,"_cell_level_spatial_colors.Rds")))
}

################################################################################
#   Plot clusters on the samples
################################################################################
message(Sys.time(), " | Plotting clusters on samples")
for(sample in samples){
  print(sample)
  spe_sub <- spe[,spe$sample_id == sample]
  for(clust_col in clust_cols){
    print(clust_col)
    p <- make_escheR(spe_sub) |>
      add_fill(clust_col) +
      scale_fill_manual(values = cluster_colors[[clust_col]])
    ggsave(
      here("plots","HD_Full_Analysis","Banksy","Cell_Level","Spatial",paste0(sample,"_",clust_col,".png")), 
      p,
      width = 10, height = 8, dpi = 200
    )
  }
}


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
