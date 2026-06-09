#Goal: Run BANKSY leiden clustering at a given resolution and plot per sample
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(HDF5Array)
library(optparse)
library(escheR)
library(Banksy)
library(here)

# Parse resolution from job script
option_list <- list(
  make_option("--res", type = "double", help = "Resolution for clusterBanksy")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$res)) {
  stop("Error: --res must be provided")
}

res <- opt$res
message(sprintf("Running clusterBanksy with resolution=%.2f", res))

# Read in sfe
sfe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis",
  "sfe_banksy_0.8_cell_level"
))

sfe

samples <- unique(sfe$sample_id)

# Cluster
message(paste0("Running Banksy clustering - ", Sys.time()))
sfe <- clusterBanksy(sfe, use_agf = TRUE, lambda = 0.8,
                     algo = "leiden", resolution = res, seed = 948)
message(paste0("Finished Banksy clustering - ", Sys.time()))

# Cluster name
clust_cols <- paste0("clust_M1_lam0.8_k50_res",res)
message("Created columns: ", clust_cols)

# Save cluster assignments
cluster_assign <- cbind(colData(sfe)[, clust_cols], rownames(colData(sfe)))
write.csv(cluster_assign,
          here("processed-data", "HD_Full_Analysis",
               paste0("sr_spatial_banksy_clusters_res_", res, ".csv")))

# Generate and save colors
cluster_colors <- Polychrome::createPalette(
  length(unique(sfe[[clust_cols]])),
  c("#D81B60", "#1E88E5", "#ffcd14", "#ff7814", "#004D40")
)

names(cluster_colors) <- unique(sfe[[clust_cols]])
saveRDS(cluster_colors,
        here("processed-data", "HD_Full_Analysis", "Cluster_colors",
             paste0(clust_cols, "_cell_level_spatial_colors_sr.Rds")))

# Plot per sample
for (sample in samples) {
  message(Sys.time(), " | Plotting: ", sample)
  sfe_sub <- sfe[, sfe$sample_id == sample]
  p <- make_escheR(sfe_sub) |>
    add_fill(clust_cols) +
    scale_fill_manual(values = cluster_colors)
  ggsave(
    here("plots", "HD_Full_Analysis", "Banksy_sr",
         "Cell_Level", "Spatial",
         paste0(sample, "_", clust_cols, ".png")),
    p, width = 20, height = 14, dpi = 200
  )
}

# Reproducibility
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
