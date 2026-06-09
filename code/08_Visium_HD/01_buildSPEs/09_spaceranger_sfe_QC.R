# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(sessioninfo)
library(HDF5Array)
library(tidyverse)
library(scuttle)
library(escheR)
library(here)


sfe_dir <- here("processed-data", "HD_Full_Analysis", "spaceranger_sfe")
out_path <- file.path(sfe_dir, "merged_sfe.rds")
sfe <- readRDS(out_path)

#Now that we have the merged sfe, only keep the non-empty cells
sfe <- sfe[, colSums(counts(sfe)) > 0]

plot_dir <- here("plots", "HD_Full_Analysis", "spaceranger_cell_QC")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- unique(sfe$sample_id)

message(sprintf(
  "Loaded cell-level SFE: %d genes x %d cells across %d samples",
  nrow(sfe), ncol(sfe), length(sample_ids)
))

sfe

################################################################################
#   Uniquify rownames to gene symbols
################################################################################

message(Sys.time(), " | Setting rownames to gene symbols via uniquifyFeatureNames")

rownames(sfe) <- uniquifyFeatureNames(
  rowData(sfe)$ID,
  rowData(sfe)$Symbol
)

sfe

# begin QC 
sub <- list(mt=grep("^MT-", rownames(sfe)))
sfe <- addPerCellQCMetrics(sfe, subsets=sub)

#compute values for spatial visualization
sfe$log10_sum <- log10(sfe$sum + 1)
sfe$log10_detected <- log10(sfe$detected + 1)


################################################################################
#   Per-sample distributional QC plots
################################################################################

message(Sys.time(), " | Generating per-sample distribution plots")

qc_cols <- c("sample_id", "sum", "detected", "subsets_mt_percent")

qc_df <- colData(sfe) |>
  as.data.frame() |>
  select(all_of(qc_cols))

#   Compute per-sample MAD thresholds for the histogram vlines
#   UMI and genes use lower-tail outliers on the log scale;
#   mito uses upper-tail outliers on the original scale
mad_thresholds <- qc_df |>
  group_by(sample_id) |>
  summarise(
    #   Lower threshold for log(sum): median - #*MAD on log scale, then back-transform
    sum_cutoff = exp(
      median(log(sum)) - 2.5 * mad(log(sum))
    ),
    #   Lower threshold for log(detected)
    detected_cutoff = exp(
      median(log(detected)) - 2.5 * mad(log(detected))
    ),
    #   Upper threshold for mito percent (not log-transformed)
    mito_cutoff = median(subsets_mt_percent, na.rm = TRUE) +
      5 * mad(subsets_mt_percent, na.rm = TRUE),
    .groups = "drop" #required to avoid invisible grouping in the table 
  )

message("Per-sample MAD thresholds:")
print(mad_thresholds)

write.csv(
  mad_thresholds,
  file.path(plot_dir, "mad_thresholds.csv")
)

#--- Histogram: UMI counts per cell (log10 scale) ---
p_umi_hist <- ggplot(qc_df, aes(x = sum + 1)) +
  geom_histogram(bins = 100, fill = "steelblue", color = "black", linewidth = 0.1) +
  geom_vline(
    data = mad_thresholds, aes(xintercept = sum_cutoff),
    color = "red", linetype = "dashed", linewidth = 0.7
  ) +
  scale_x_log10() +
  facet_wrap(~ sample_id, scales = "free_y") +
  labs(
    title = "UMI counts per cell (red = 2.5 MAD lower cutoff)",
    x = "Total UMI (log10)", y = "Number of cells"
  ) +
  theme_bw()

ggsave(
  file.path(plot_dir, "hist_sum_umi.png"), p_umi_hist,
  width = 10, height = 8, dpi = 150
)

#--- Histogram: Genes detected per cell ---
p_gene_hist <- ggplot(qc_df, aes(x = detected + 1)) +
  geom_histogram(bins = 100, fill = "darkgreen", color = "black", linewidth = 0.1) +
  geom_vline(
    data = mad_thresholds, aes(xintercept = detected_cutoff),
    color = "red", linetype = "dashed", linewidth = 0.7
  ) +
  scale_x_log10() +
  facet_wrap(~ sample_id, scales = "free_y") +
  labs(
    title = "Genes detected per cell (red = 2.5 MAD lower cutoff)",
    x = "Genes detected (log10)", y = "Number of cells"
  ) +
  theme_bw()

ggsave(
  file.path(plot_dir, "hist_sum_gene.png"), p_gene_hist,
  width = 10, height = 8, dpi = 150
)

#--- Histogram: Mito percent ---
p_mito_hist <- ggplot(qc_df, aes(x = subsets_mt_percent)) +
  geom_histogram(bins = 100, fill = "firebrick", color = "black", linewidth = 0.1) +
  geom_vline(
    data = mad_thresholds, aes(xintercept = mito_cutoff),
    color = "red", linetype = "dashed", linewidth = 0.7
  ) +
  facet_wrap(~ sample_id, scales = "free_y") +
  labs(
    title = "Mitochondrial percent per cell (red = 5 MAD upper cutoff)",
    x = "Mito UMI %", y = "Number of cells"
  ) +
  theme_bw()

ggsave(
  file.path(plot_dir, "hist_mito_percent.png"), p_mito_hist,
  width = 10, height = 8, dpi = 150
)


#--- Scatter: genes vs UMI, colored by mito percent ---
p_scatter <- ggplot(
  qc_df |> slice_sample(n = min(nrow(qc_df), 200000)),
  aes(x = sum, y = detected, color = subsets_mt_percent)
) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_viridis_c(option = "inferno", limits = c(0, 50), oob = scales::squish) +
  facet_wrap(~ sample_id) +
  labs(
    title = "Genes detected vs UMI (colored by mito %)",
    x = "Total UMI (log10)", y = "Genes detected (log10)",
    color = "Mito %"
  ) +
  theme_bw()

ggsave(
  file.path(plot_dir, "scatter_gene_vs_umi.png"), p_scatter,
  width = 10, height = 8, dpi = 150
)

#--- Per-sample summary table ---
summary_cols <- rlang::exprs(
  n_cells = n(),
  median_umi = median(sum),
  median_genes = median(detected),
  median_mito_pct = median(subsets_mt_percent, na.rm = TRUE),
  pct_below_5_umi = mean(sum < 5) * 100,
  pct_mito_above_20 = mean(subsets_mt_percent > 20, na.rm = TRUE) * 100
)

summary_df <- qc_df |>
  group_by(sample_id) |>
  summarise(!!!summary_cols, .groups = "drop") |>
  left_join(mad_thresholds, by = "sample_id") |>
  as.data.frame()

write.csv(
  summary_df,
  file.path(plot_dir, "sample_qc_summary.csv"),
  row.names = FALSE
)

message("Per-sample QC summary:")
print(summary_df)



################################################################################
#   Spatial QC plots with escheR
################################################################################

message(Sys.time(), " | Generating spatial QC plots with escheR")

escher_spot_size <- 0.5
escher_point_size <- 0.5

for (sid in sample_ids) {
  message(sprintf("  Plotting sample: %s", sid))
  
  sfe_sub <- sfe[, sfe$sample_id == sid]
  
  #   Trim colData to only QC-relevant columns to sfeed up escheR plotting
  keep_cols <- c(
    "sample_id", "sum", "detected", "log10_sum", "log10_detected",
    "subsets_mt_sum", "subsets_mt_percent"
  )
  
  keep_cols <- intersect(keep_cols, colnames(colData(sfe_sub)))
  colData(sfe_sub) <- colData(sfe_sub)[, keep_cols]
  
  #--- Spatial: log10(UMI counts) ---
  p1 <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
    add_fill(var = "log10_sum", point_size = escher_point_size) +
    scale_fill_gradientn(
      colors = viridisLite::plasma(256),
      name = "log10(UMI)"
    ) +
    ggtitle(paste(sid, "- log10(UMI)")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  ggsave(
    file.path(plot_dir, sprintf("%s_spatial_umi.png", sid)), p1,
    width = 10, height = 8, dpi = 200
  )
  
  #--- Spatial: Mito percent ---
  p2 <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
    add_fill(var = "subsets_mt_percent", point_size = escher_point_size) +
    scale_fill_gradientn(
      colors = viridisLite::inferno(256),
      limits = c(0, 50),
      oob = scales::squish,
      name = "Mito %"
    ) +
    ggtitle(paste(sid, "- Mito %")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  ggsave(
    file.path(plot_dir, sprintf("%s_spatial_mito.png", sid)), p2,
    width = 10, height = 8, dpi = 200
  )
  
  #--- Spatial: log10(Genes detected) ---
  p3 <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
    add_fill(var = "log10_detected", point_size = escher_point_size) +
    scale_fill_gradientn(
      colors = viridisLite::viridis(256),
      name = "log10(Genes)"
    ) +
    ggtitle(paste(sid, "- log10(Genes detected)")) +
    theme_void() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      legend.position = "right"
    )
  
  ggsave(
    file.path(plot_dir, sprintf("%s_spatial_genes.png", sid)), p3,
    width = 10, height = 8, dpi = 200
  )
  
  #--- Spatial: log10(Bin count) (if available) ---
  if ("log10_bin_count" %in% keep_cols) {
    p4 <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
      add_fill(var = "log10_bin_count", point_size = escher_point_size) +
      scale_fill_gradientn(
        colors = viridisLite::mako(256),
        name = "log10(Bins)"
      ) +
      ggtitle(paste(sid, "- log10(Bins per cell)")) +
      theme_void() +
      theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        legend.position = "right"
      )
    
    ggsave(
      file.path(plot_dir, sprintf("%s_spatial_bin_count.png", sid)), p4,
      width = 10, height = 8, dpi = 200
    )
  }
}

################################################################################
#   Flag outlier cells using MAD-based detection
################################################################################

message(Sys.time(), " | Flagging outlier cells (per sample, MAD-based)")

#   Run outlier detection per sample to account for sample-level differences
sfe$qc_low_umi <- FALSE
sfe$qc_low_gene <- FALSE
sfe$qc_high_mito <- FALSE
sfe$qc_discard <- FALSE

for (sid in sample_ids) {
  idx <- sfe$sample_id == sid
  
  low_umi <- isOutlier(
    sfe$sum[idx], type = "lower", log = TRUE, nmads = 2.5
  )
  low_gene <- isOutlier(
    sfe$detected[idx], type = "lower", log = TRUE, nmads = 2.5
  )
  high_mito <- isOutlier(
    sfe$subsets_mt_percent[idx], type = "higher", nmads = 5
  )
  
  
  sfe$qc_low_umi[idx] <- low_umi
  sfe$qc_low_gene[idx] <- low_gene
  sfe$qc_high_mito[idx] <- high_mito
  sfe$qc_discard[idx] <- low_umi | low_gene | high_mito 
  
  message(sprintf(
    "  %s: %d/%d cells flagged (%.1f%%) — low_umi: %d, low_gene: %d, high_mito: %d",
    sid, sum(low_umi | low_gene | high_mito), sum(idx),
    mean(low_umi | low_gene | high_mito ) * 100,
    sum(low_umi), sum(low_gene), sum(high_mito)
  ))
}


#--- Spatial plot of flagged cells with escheR ---
#   Convert to factor for categorical coloring with add_fill
sfe$qc_discard_factor <- factor(
  sfe$qc_discard, levels = c(FALSE, TRUE), labels = c("Pass", "Flagged")
)

message(Sys.time(), " | Plotting QC-flagged cells spatially")

for (sid in sample_ids) {
  sfe_sub <- sfe[, sfe$sample_id == sid]
  
  #   Trim colData for sfeed
  colData(sfe_sub) <- colData(sfe_sub)[, c("sample_id", "qc_discard_factor")]
  
  p_flag <- make_escheR(sfe_sub, spot_size = escher_spot_size) |>
    add_fill(var = "qc_discard_factor", point_size = escher_point_size) +
    scale_fill_manual(
      values = c("Pass" = "grey80", "Flagged" = "red"),
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

################################################################################
#   Save updated sfe with QC columns
################################################################################

message(Sys.time(), " | Saving sfe with QC annotations")

sfe_qc_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_with_QC"
)
sfe <- saveHDF5SummarizedExperiment(
  sfe, dir = sfe_qc_dir, replace = TRUE, as.sparse = TRUE
)

message(sprintf("Saved to: %s", sfe_qc_dir))
message(sprintf(
  "Total cells flagged: %d / %d (%.1f%%)",
  sum(sfe$qc_discard), ncol(sfe), mean(sfe$qc_discard) * 100
))

################################################################################
#   Remove discarded cells and save filtered sfe
################################################################################

message(Sys.time(), " | Filtering out QC-flagged cells")

n_before <- ncol(sfe)
sfe_filtered <- sfe[, !sfe$qc_discard]
n_after <- ncol(sfe_filtered)

message(sprintf(
  "  Removed %d cells: %d -> %d (%.1f%% retained)",
  n_before - n_after, n_before, n_after,
  (n_after / n_before) * 100
))

#   Per-sample breakdown
for (sid in sample_ids) {
  idx_before <- sum(sfe$sample_id == sid)
  idx_after <- sum(sfe_filtered$sample_id == sid)
  message(sprintf(
    "  %s: %d -> %d cells (%.1f%% retained)",
    sid, idx_before, idx_after, (idx_after / idx_before) * 100
  ))
}

#   Drop the QC flag columns that are no longer needed
sfe_filtered$qc_discard_factor <- NULL

sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered"
)
sfe_filtered <- saveHDF5SummarizedExperiment(
  sfe_filtered, dir = sfe_filtered_dir, replace = TRUE, as.sparse = TRUE
)

message(sprintf("Saved filtered sfe to: %s", sfe_filtered_dir))

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
