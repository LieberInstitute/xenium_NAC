library(SpatialExperiment)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(here)


files <- list.files(path = here("processed-data","HD_Full_Analysis",
                                "Banksy","NonSpatial_colData"),full.names = TRUE)

col_data_list <- lapply(files, readRDS)

names(col_data_list) <- list.files(path = here("processed-data","HD_Full_Analysis",
                                               "Banksy","NonSpatial_colData"))

names(col_data_list) <- gsub(x = names(col_data_list),pattern = ".Rds",replacement = "")

#Is everything identical? 
all(sapply(col_data_list, function(x) identical(rownames(x), rownames(col_data_list[[1]]))))
stopifnot(all(sapply(col_data_list, function(x) identical(rownames(x), rownames(col_data_list[[1]])))))


cluster_df <- data.frame(
  lapply(col_data_list, function(x) x[[24]])  # first column of each
)

colnames(cluster_df) <- names(col_data_list) 

rownames(cluster_df) <- rownames(col_data_list[[1]])


#Load SPE
spe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered_v2"
)

spe <- loadHDF5SummarizedExperiment(spe_filtered_dir)

stopifnot(identical(rownames(cluster_df),colnames(spe)))

#Add clusters to the object
colData(spe) <- cbind(colData(spe), cluster_df)


cluster_cols <- vector(mode = "list",length = ncol(cluster_df))
names(cluster_cols) <- colnames(cluster_df)

#Generate colors
for(cluster in names(cluster_cols)){
  print(cluster)
  #Generate and save colors
  cluster_cols[[cluster]] <- Polychrome::createPalette(
    length(unique(spe[[cluster]])),
    c("#D81B60", "#1E88E5", "#ffcd14", "#ff7814", "#004D40")
  )
  names(cluster_cols[[cluster]]) <- as.factor(unique(spe[[cluster]]))
  saveRDS(cluster_cols[[cluster]],
          here("processed-data", "HD_Full_Analysis", "Cluster_colors",
               paste0(cluster, "_cell_level_non-spatial_colors.Rds")))
}

#Make plots
for(sample in unique(spe$sample_id)){
  print(sample)
  spe_sub <- spe[,spe$sample_id == sample]
  for(cluster in names(cluster_cols)){
    p <- make_escheR(spe_sub) |>
      add_fill(cluster) +
      scale_fill_manual(values = cluster_cols[[cluster]])
    ggsave(plot = p,
           filename = here("plots","HD_Full_Analysis","Banksy",
                           "Cell_Level","Nonspatial",
                           paste0(cluster,"_",sample,".png")),
           width = 20, height = 14, dpi = 200)
  }
}


#Plot each cluster singulalry on each sample
for (sample in unique(spe$sample_id)) {
  print(sample)
  dir.create(path = here("plots","HD_Full_Analysis","Banksy",
                         "Cell_Level","Nonspatial",sample))
  spe_sub <- spe[, spe$sample_id == sample]
  for (cluster in names(cluster_cols)) {
    print(cluster)
    dir.create(path = here("plots","HD_Full_Analysis","Banksy",
                           "Cell_Level","Nonspatial",sample,cluster))
    for (ct in unique(spe_sub[[cluster]])) {
      print(ct)
      # Create a highlight column: this cell type gets its color, all others grey
      spe_sub$highlight <- ifelse(spe_sub[[cluster]] == ct, as.character(ct), "Other")
      
      # Build color vector: the cell type color + grey for Other
      plot_colors <- c(cluster_cols[[cluster]][as.character(ct)], "Other" = "grey90")
      
      p <- make_escheR(spe_sub) |>
        add_fill("highlight") +
        scale_fill_manual(values = plot_colors) +
        ggtitle(paste(sample, cluster, "- Cluster", ct))
      
      ggsave(filename = here("plots","HD_Full_Analysis","Banksy",
                             "Cell_Level","Nonspatial",sample,cluster,
                             paste0(ct,".png")),
             plot = p,
             width = 20, height = 14, dpi = 200)
    }
  }
}


sessionInfo()



