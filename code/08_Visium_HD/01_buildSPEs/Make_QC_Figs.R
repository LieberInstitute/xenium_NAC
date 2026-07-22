# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(sessioninfo)
library(HDF5Array)
library(tidyverse)
library(scuttle)
library(escheR)
library(here)


sfe_dir <- here("processed-data","HD_Full_Analysis","spaceranger_sfe")
out_path <- file.path(sfe_dir, "merged_sfe.rds")
sfe <- readRDS(out_path)

#Now that we have the merged sfe, only keep the non-empty cells
sfe <- sfe[, colSums(counts(sfe)) > 0]

plot_dir <- here("plots", "HD_Full_Analysis", "spaceranger_cell_QC_supp_figure")
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
  file.path(plot_dir, "scatter_gene_umi.png"), p_scatter,
  width = 10, height = 8, dpi = 150
)


################################################################################
#   Per-sample, per-metric QC PDFs (for figure assembly)
#
#   One PDF per sample per metric (4 metrics x N samples), sized consistently
#   so they can be dropped straight into a multi-panel figure (patchwork,
#   cowplot, Illustrator, etc.).
################################################################################

message(Sys.time(), " | Generating per-sample, per-metric QC PDFs")

pdf_dir <- file.path(plot_dir, "per_sample_pdfs")
dir.create(pdf_dir, recursive = TRUE, showWarnings = FALSE)

for (sid in sample_ids) {
  
  df_i  <- qc_df |> filter(sample_id == sid)
  thr_i <- mad_thresholds |> filter(sample_id == sid)
  
  #--- UMI counts per cell ---
  p <- ggplot(df_i, aes(x = sum + 1)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black", linewidth = 0.1) +
    geom_vline(xintercept = thr_i$sum_cutoff, color = "red", linetype = "dashed", linewidth = 0.7) +
    scale_x_log10() +
    labs(
      title = sid, subtitle = "UMI counts per cell (red = 2.5 MAD lower cutoff)",
      x = "Total UMI (log10)", y = "Number of cells"
    ) +
    theme_bw()
  ggsave(file.path(pdf_dir, sprintf("hist_sum_umi_%s.pdf", sid)), p, width = 5, height = 4)
  
  #--- Genes detected per cell ---
  p <- ggplot(df_i, aes(x = detected + 1)) +
    geom_histogram(bins = 100, fill = "darkgreen", color = "black", linewidth = 0.1) +
    geom_vline(xintercept = thr_i$detected_cutoff, color = "red", linetype = "dashed", linewidth = 0.7) +
    scale_x_log10() +
    labs(
      title = sid, subtitle = "Genes detected per cell (red = 2.5 MAD lower cutoff)",
      x = "Genes detected (log10)", y = "Number of cells"
    ) +
    theme_bw()
  ggsave(file.path(pdf_dir, sprintf("hist_sum_gene_%s.pdf", sid)), p, width = 5, height = 4)
  
  #--- Mito percent per cell ---
  p <- ggplot(df_i, aes(x = subsets_mt_percent)) +
    geom_histogram(bins = 100, fill = "firebrick", color = "black", linewidth = 0.1) +
    geom_vline(xintercept = thr_i$mito_cutoff, color = "red", linetype = "dashed", linewidth = 0.7) +
    labs(
      title = sid, subtitle = "Mitochondrial percent per cell (red = 5 MAD upper cutoff)",
      x = "Mito UMI %", y = "Number of cells"
    ) +
    theme_bw()
  ggsave(file.path(pdf_dir, sprintf("hist_mito_percent_%s.pdf", sid)), p, width = 5, height = 4)
  
  #--- Genes vs UMI scatter, colored by mito percent ---
  df_scatter <- df_i |> slice_sample(n = min(nrow(df_i), 200000))
  p <- ggplot(df_scatter, aes(x = sum, y = detected, color = subsets_mt_percent)) +
    geom_point(size = 0.3, alpha = 0.3) +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_viridis_c(option = "inferno", limits = c(0, 50), oob = scales::squish) +
    labs(
      title = sid, subtitle = "Genes detected vs UMI (colored by mito %)",
      x = "Total UMI (log10)", y = "Genes detected (log10)", color = "Mito %"
    ) +
    theme_bw()
  ggsave(file.path(pdf_dir, sprintf("scatter_gene_umi_%s.pdf", sid)), p, width = 5, height = 4.5)
  
  message(sprintf("  wrote PDFs for sample %s", sid))
}

message(Sys.time(), " | Done. PDFs written to: ", pdf_dir)


################################################################################
#   Reproducibility
################################################################################

session_info()
