# Combined visualization of drug/compound results:
#   1. drug2cell: top enriched drugs per spatial domain / cell type (dotplot)
#      + highlight panel of drugs of abuse / psychiatric drugs
#   2. Curated target sets: enrichment heatmaps
#   3. LINCS: top mimicking / reversing compounds per domain
#
# Usage: sbatch 06_drug_viz.sh

library(here)
library(data.table)
library(ggplot2)
library(scales)

d2c_dir <- here("processed-data", "17_drug_analysis", "drug2cell")
cur_dir <- here("processed-data", "17_drug_analysis", "curated_sets")
plot_dir <- here("plots", "17_drug_analysis")

## ---- 0. Shared color scale ---------------------------------------------------
# Color/fill values outside +/- LOGFC_LIM are squished onto the endpoints so a
# few extreme drugs don't wash out the rest of the panel. Clipping is
# display-only: ranking and top-N selection still use the unclipped logFC.
LOGFC_LIM <- 5

scale_logfc <- function(aesthetic = c("colour", "fill"),
                        lim = LOGFC_LIM, mid = "grey85", name = "logFC") {
  aesthetic <- match.arg(aesthetic)
  f <- if (aesthetic == "colour") scale_color_gradient2 else scale_fill_gradient2
  f(low = "#2166AC", mid = mid, high = "#B2182B",
    limits = c(-lim, lim),
    oob = scales::squish,
    breaks = seq(-lim, lim, length.out = 5),
    labels = function(x) {
      out <- format(x, trim = TRUE)
      out[!is.na(x) & x <= -lim] <- sprintf("\u2264 -%g", lim)
      out[!is.na(x) & x >=  lim] <- sprintf("\u2265 %g", lim)
      out
    },
    name = name)
}

clip_caption <- function(lim = LOGFC_LIM) {
  sprintf("logFC clipped at \u00b1%g for display", lim)
}

## ---- 1. drug2cell -----------------------------------------------------------
plot_d2c <- function(group_col, n_top = 8) {
  res <- fread(file.path(d2c_dir, sprintf("drug_enrichment_%s.tsv.gz", group_col)))
  res[, drug_name := sub("^[^|]*\\|", "", drug)]

  # Top N drugs per group by absolute logFC (unclipped)
  top <- res[order(-abs(logfoldchanges)), head(.SD, n_top), by = group_col]
  top[, drug_name := factor(drug_name, levels = rev(unique(drug_name)))]
  top[, neg_log10_fdr := pmin(-log10(pvals_adj + 1e-300), 300)]

  p <- ggplot(top, aes(.data[[group_col]], drug_name,
                       size = neg_log10_fdr, color = logfoldchanges)) +
    geom_point() +
    scale_logfc("colour") +
    labs(
      x = NULL, y = NULL, size = expression(-log[10]~FDR),
      title = sprintf("Top drug2cell drugs by %s", group_col),
      caption = clip_caption()
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plot_dir, sprintf("drug2cell_top_%s.pdf", group_col)),
         p, width = 9, height = 1 + 0.16 * length(unique(top$drug_name)),
         limitsize = FALSE)

  # Highlight panel (drugs of abuse / psychiatric drugs)
  hl <- fread(file.path(d2c_dir, sprintf("drug_enrichment_%s_highlight.tsv", group_col)))
  hl[, drug_name := sub("^[^|]*\\|", "", drug)]
  # One representative row (max |logFC|) per drug name x group
  hl <- hl[order(-abs(logfoldchanges)), head(.SD, 1), by = c("drug_name", group_col)]

  p2 <- ggplot(hl, aes(.data[[group_col]], drug_name, fill = logfoldchanges)) +
    geom_tile() +
    facet_grid(highlight_class ~ ., scales = "free_y", space = "free_y") +
    scale_logfc("fill", mid = "white") +
    labs(x = NULL, y = NULL,
         title = sprintf("Drugs of abuse & psychiatric drugs by %s", group_col),
         caption = clip_caption()) +
    theme_bw(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
  ggsave(file.path(plot_dir, sprintf("drug2cell_highlight_%s.pdf", group_col)),
         p2, width = 9, height = 2 + 0.13 * length(unique(hl$drug_name)),
         limitsize = FALSE)
}

plot_d2c("Spatial_Domain")
plot_d2c("labels")

# Highlighted drugs only, dotplot faceted by drug class
plot_d2c_highlight_dot <- function(group_col, n_top = 5) {
  hl <- fread(file.path(d2c_dir, sprintf("drug_enrichment_%s_highlight.tsv", group_col)))
  hl[, drug_name := sub("^[^|]*\\|", "", drug)]
  # One representative row (max |logFC|) per drug name x group
  hl <- hl[order(-abs(logfoldchanges)), head(.SD, 1), by = c("drug_name", group_col)]
  # Top N drugs per group within each class
  top <- hl[order(-abs(logfoldchanges)), head(.SD, n_top), by = c("highlight_class", group_col)]
  top[, drug_name := reorder(drug_name, logfoldchanges)]
  top[, neg_log10_fdr := pmin(-log10(pvals_adj + 1e-300), 300)]

  p <- ggplot(top, aes(.data[[group_col]], drug_name,
                       size = neg_log10_fdr, color = logfoldchanges)) +
    geom_point() +
    facet_grid(highlight_class ~ ., scales = "free_y", space = "free_y") +
    scale_logfc("colour") +
    labs(
      x = NULL, y = NULL, size = expression(-log[10]~FDR),
      title = sprintf("Highlighted drugs by %s (grouped by class)", group_col),
      caption = clip_caption()
    ) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text.y = element_text(angle = 0)
    )
  ggsave(file.path(plot_dir, sprintf("drug2cell_highlight_dot_%s.pdf", group_col)),
         p, width = 10, height = 2 + 0.16 * length(unique(top$drug_name)),
         limitsize = FALSE)
}

plot_d2c_highlight_dot("Spatial_Domain")
plot_d2c_highlight_dot("labels")

## ---- 2. Curated target sets --------------------------------------------------
for (group_col in c("Spatial_Domain", "labels")) {
  cur_file <- file.path(cur_dir, sprintf("curated_set_enrichment_%s.csv", group_col))
  if (!file.exists(cur_file)) {
    message(sprintf("Curated set results not found for %s; skipping.", group_col))
    next
  }
  res <- fread(cur_file)
  res[, sig := fifelse(wilcox_fdr < 0.05, "*", "")]

  # delta_median is a z-score difference, so +/- 3 is far too wide here.
  # Round the observed max up to the next 0.1; replace with a fixed value
  # (e.g. d_lim <- 1) once the range is stable across runs.
  d_lim <- ceiling(max(abs(res$delta_median), na.rm = TRUE) * 10) / 10
  d_lim <- max(d_lim, 0.1)

  p <- ggplot(res, aes(group, drug_set, fill = delta_median)) +
    geom_tile() +
    geom_text(aes(label = sig), size = 3.5) +
    scale_logfc("fill", lim = d_lim, mid = "white",
                name = "\u0394 median\nz-score") +
    labs(x = NULL, y = NULL,
         title = sprintf("Curated drug-target set enrichment by %s", group_col),
         subtitle = "* Wilcoxon (group vs rest) FDR < 0.05",
         caption = clip_caption(d_lim)) +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  ggsave(file.path(plot_dir, sprintf("curated_sets_%s.pdf", group_col)),
         p, width = 8, height = 5)
}

sessioninfo::session_info()
