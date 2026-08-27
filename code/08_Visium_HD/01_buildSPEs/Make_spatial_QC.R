# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(sessioninfo)
library(HDF5Array)
library(tidyverse)
library(scuttle)
library(escheR)
library(here)

##Load the sfe with QC
sfe_qc_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_with_QC"
)
sfe <- loadHDF5SummarizedExperiment(sfe_qc_dir)

plot_dir <- here("plots", "HD_Full_Analysis", "spaceranger_cell_QC_spatial")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

#--- Spatial plot of flagged cells with escheR ---
#   Convert to factor for categorical coloring with add_fill
sfe$qc_discard_factor <- factor(
  sfe$qc_discard, levels = c(FALSE, TRUE), labels = c("Pass", "Flagged")
)

message(Sys.time(), " | Plotting QC-flagged cells spatially")

sample_ids <- unique(sfe$sample_id)
for (sid in sample_ids) {
  sfe_sub <- sfe[, sfe$sample_id == sid]
  
  #   Trim colData for sfeed
  colData(sfe_sub) <- colData(sfe_sub)[, c("sample_id", "qc_discard_factor")]
  
  p_flag <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
    add_fill(var = "qc_discard_factor",point_size = 1) +
    scale_fill_manual(
      values = c("Pass" = "grey50", "Flagged" = "red"),
      name = "QC Status"
    ) +
    ggtitle(paste(sid, "- QC-flagged cells")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  ggsave(
    file.path(plot_dir, sprintf("%s_spatial_qc_flagged.png", sid)), p_flag,
    width = 10, height = 8, dpi = 200
  )
}



sessioninfo()
