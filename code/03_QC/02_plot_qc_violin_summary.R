#!/usr/bin/env Rscript
# Recreate the four QC violin plots and the per-sample QC summary plot from
# saved QC outputs. This script does not recalculate or change QC flags.

library(SpatialExperiment)
library(scater)
library(tidyverse)

# Resolve paths from this script rather than from the R working directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) != 1) {
  stop("Run this file with Rscript so its project path can be determined.")
}
script_path <- normalizePath(sub("^--file=", "", script_arg))
code_dir <- dirname(dirname(script_path))
project_dir <- dirname(code_dir)

spe_path <- file.path(
  project_dir,
  "processed-data", "02_build_spe", "SPEs", "spe_withQC.Rds"
)
spe <- readRDS(spe_path)

required_cols <- c(
  "Sample",
  "subsets_any_neg_percent",
  "exclude_any_neg",
  "sum_gex",
  "sum_gex_4MAD",
  "detected_gex",
  "detected_gex_4MAD",
  "cell_area",
  "cell_area_4MAD",
  "global_outliers"
)
missing_cols <- setdiff(required_cols, colnames(colData(spe)))
if (length(missing_cols) > 0) {
  stop("Missing required colData columns: ", paste(missing_cols, collapse = ", "))
}

plotdir <- file.path(project_dir, "plots", "03_qc")
dir.create(plotdir, recursive = TRUE, showWarnings = FALSE)

# Save plots with a fixed data-panel size. This keeps the plotting region
# identical even when legend widths differ or a legend is hidden.
save_qc_plot <- function(p, filename, panel_width = 5.5,
                         panel_height = 5.5, dpi = 300) {
  gt <- ggplotGrob(p)
  panel_rows <- unique(gt$layout$t[grepl("^panel", gt$layout$name)])
  panel_cols <- unique(gt$layout$l[grepl("^panel", gt$layout$name)])

  gt$widths[panel_cols] <- grid::unit(
    panel_width / length(panel_cols), "in"
  )
  gt$heights[panel_rows] <- grid::unit(
    panel_height / length(panel_rows), "in"
  )

  total_width <- grid::convertWidth(
    sum(gt$widths), "in", valueOnly = TRUE
  )
  total_height <- grid::convertHeight(
    sum(gt$heights), "in", valueOnly = TRUE
  )

  ggsave(
    filename = filename,
    plot = gt,
    width = total_width,
    height = total_height,
    dpi = dpi,
    limitsize = FALSE,
    bg = "white"
  )
}

plot_qc_violin <- function(value_col, flag_col, y_label, filename,
                           log_y = FALSE, hide_legend = FALSE) {
  p <- plotColData(
    spe,
    x = "Sample",
    y = value_col,
    color_by = flag_col
  ) +
    stat_summary(
      fun = median,
      fun.min = median,
      fun.max = median,
      geom = "crossbar",
      width = 0.3
    ) +
    labs(
      x = "Sample",
      y = y_label
    ) +
    scale_x_discrete(expand = expansion(add = 0.8)) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.margin = margin(5.5, 5.5, 5.5, 25)
    )

  if (log_y) {
    p <- p + scale_y_log10()
  }
  if (hide_legend) {
    p <- p + theme(legend.position = "none")
  }

  save_qc_plot(
    p,
    file.path(plotdir, filename)
  )
}

plot_qc_violin(
  value_col = "subsets_any_neg_percent",
  flag_col = "exclude_any_neg",
  y_label = "Negative control and unassigned transcripts (%)",
  filename = "exclude_any_neg_outliers_violin.png"
)

plot_qc_violin(
  value_col = "sum_gex",
  flag_col = "sum_gex_4MAD",
  y_label = "Gene expression transcripts per cell",
  filename = "sum_gex_4MAD_outliers_violin.png",
  log_y = TRUE
)

plot_qc_violin(
  value_col = "detected_gex",
  flag_col = "detected_gex_4MAD",
  y_label = "Detected gene expression genes per cell",
  filename = "detected_gex_4MAD_outliers_violin.png",
  log_y = TRUE
)

plot_qc_violin(
  value_col = "cell_area",
  flag_col = "cell_area_4MAD",
  y_label = expression("Cell area ("*mu*"m"^2*")"),
  filename = "cell_area_4MAD_outliers_violin.png",
  hide_legend = TRUE
)

# Rebuild the summary directly from the saved cell-level QC columns.
summary_df <- as.data.frame(colData(spe)) %>%
  group_by(Sample) %>%
  summarize(
    CellsRetained = sum(!global_outliers, na.rm = TRUE),
    Outliers = sum(global_outliers, na.rm = TRUE),
    PercentRemoved = Outliers / (CellsRetained + Outliers) * 100,
    exclude_any_neg = sum(exclude_any_neg, na.rm = TRUE),
    sum_gex_4MAD = sum(sum_gex_4MAD, na.rm = TRUE),
    detected_gex_4MAD = sum(detected_gex_4MAD, na.rm = TRUE),
    cell_area_4MAD = sum(cell_area_4MAD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))

metric_levels <- c(
  "Negative controls ≥25% (cells)",
  "Low GEX transcript count (lower 4 MAD; cells)",
  "Low detected GEX genes (lower 4 MAD; cells)",
  "Cell area outliers (two sided 4 MAD; cells)",
  "Cells removed (%)"
)

summary_long <- summary_df %>%
  select(
    Sample,
    exclude_any_neg,
    sum_gex_4MAD,
    detected_gex_4MAD,
    cell_area_4MAD,
    PercentRemoved
  ) %>%
  pivot_longer(
    cols = -Sample,
    names_to = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = recode(
      Metric,
      exclude_any_neg = metric_levels[[1]],
      sum_gex_4MAD = metric_levels[[2]],
      detected_gex_4MAD = metric_levels[[3]],
      cell_area_4MAD = metric_levels[[4]],
      PercentRemoved = metric_levels[[5]]
    ),
    Metric = factor(Metric, levels = metric_levels)
  )

p_summary <- ggplot(summary_long, aes(x = Sample, y = Value)) +
  geom_col() +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 20)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Per sample QC flags and cells removed"
  )

ggsave(
  filename = file.path(plotdir, "outlier_summary_hist.png"),
  plot = p_summary,
  width = 7,
  height = 10,
  dpi = 300,
  bg = "white"
)
