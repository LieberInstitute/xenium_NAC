# Diagnostic and summary plots for the scDRS risk-gene driver analysis
# (outputs of 09_risk_gene_drivers.py).
#
# Intended to be run INTERACTIVELY, e.g.:
#   source("code/16_scDRS/10_driver_gene_plots.R")
#   dr <- load_drivers("Spatial_Domain")
#   driver_scatter(dr, "AUD")
#   driver_heatmap(dr, "AUD")
#   driver_specificity(dr, "AUD")
#   weight_percentile(dr)
#   ex <- load_expression("Spatial_Domain")
#   st <- load_group_stats("Spatial_Domain")
#   expression_heatmap(ex, "SCZ", stats = st)
#   expression_dotplot(ex, "SCZ", stats = st)
#   save_plot(driver_scatter(dr, "SCZ"), "driver_scatter_SCZ")
#
# ---------------------------------------------------------------------------
# HOW TO READ THESE PLOTS -- three caveats established during QC:
#
# 1. spearman_r is bounded near 1/sqrt(1000) ~= 0.032, because the scDRS score
#    is a weighted mean of ~1000 genes and any single gene contributes ~0.1%.
#    An r of 0.06 is near the arithmetic maximum, NOT weak evidence.
# 2. r is not comparable ACROSS groups: larger groups estimate it more
#    precisely (coupling correlates ~0.5 with log n_cells), and small groups
#    show inflated maxima. Compare gene RANKINGS within a group.
# 3. The trend of r against magma_weight is mechanical (high-weight genes
#    contribute more to the score by construction). The signal is the
#    DEVIATION from that trend, not the trend itself.
# ---------------------------------------------------------------------------

library(data.table)
library(ggplot2)

BASE <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC"
DRIVER_DIR <- file.path(BASE, "processed-data", "16_scDRS", "risk_gene_drivers")
PLOT_DIR <- file.path(BASE, "plots", "16_scDRS", "drivers")
dir.create(PLOT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N <- 20
HAS_REPEL <- requireNamespace("ggrepel", quietly = TRUE)


#' Load a driver table and attach a within-group rank.
#' @param grouping "Spatial_Domain" or "labels"
load_drivers <- function(grouping = c("Spatial_Domain", "labels")) {
  grouping <- match.arg(grouping)
  d <- fread(file.path(DRIVER_DIR, sprintf("risk_gene_drivers_%s.tsv.gz", grouping)))
  d <- d[!is.na(spearman_r)]
  d[, rank := frank(-spearman_r, ties.method = "first"), by = .(trait, group)]
  # Order group levels by size so facets read largest-first.
  lv <- unique(d[, .(group, n_cells)])[order(-n_cells), group]
  d[, group := factor(group, levels = lv)]
  setattr(d, "grouping", grouping)
  d[]
}


#' Scatter of gene-level tracking against GWAS weight, faceted by group.
#' Points above the local trend track risk more than their weight predicts.
driver_scatter <- function(d, trait_key, top_n = TOP_N, ncol = 5) {
  x <- copy(d[trait == trait_key])
  if (!nrow(x)) stop("no rows for trait ", trait_key)
  x[, is_top := rank <= top_n]

  p <- ggplot(x, aes(magma_weight, spearman_r)) +
    geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey60") +
    geom_point(aes(colour = is_top), alpha = 0.4, size = 0.7) +
    scale_colour_manual(values = c(`FALSE` = "grey75", `TRUE` = "#B2182B"),
                        labels = c("other", paste0("top ", top_n)), name = NULL) +
    facet_wrap(~group, scales = "free_y", ncol = ncol) +
    labs(x = "MAGMA gene z-score (GWAS weight)",
         y = "Spearman r (expression vs per-cell scDRS score)",
         title = sprintf("%s: risk-gene drivers by %s", trait_key,
                         attr(d, "grouping")),
         subtitle = "Positive trend is expected by construction; genes above it track risk more than weight predicts. y-scales differ per facet.")

  lab <- x[is_top == TRUE]
  if (HAS_REPEL) {
    p + ggrepel::geom_text_repel(data = lab, aes(label = gene), size = 1.8,
                                 max.overlaps = 40, segment.size = 0.15,
                                 min.segment.length = 0, colour = "black")
  } else {
    p + geom_text(data = lab, aes(label = gene), size = 1.8, vjust = -0.7)
  }
}


#' Heatmap of the union of top-N driver genes across groups, for one trait.
#' Shows whether a gene drives risk everywhere or only in specific contexts.
driver_heatmap <- function(d, trait_key, top_n = TOP_N) {
  x <- d[trait == trait_key]
  keep_genes <- unique(x[rank <= top_n, gene])
  h <- copy(x[gene %in% keep_genes])
  # order genes by how consistently they rank highly
  ord <- h[, .(m = mean(spearman_r)), by = gene][order(m), gene]
  h[, gene := factor(gene, levels = ord)]

  ggplot(h, aes(group, gene, fill = spearman_r)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, name = "Spearman r") +
    labs(x = NULL, y = NULL,
         title = sprintf("%s: top-%d driver genes across %s", trait_key, top_n,
                         attr(d, "grouping")),
         subtitle = sprintf("%d genes = union of each group's top %d",
                            length(keep_genes), top_n)) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 6))
}


#' How group-specific are the top drivers? Counts, for each gene in any
#' group's top-N, how many groups it appears in.
driver_specificity <- function(d, trait_key, top_n = TOP_N) {
  x <- d[trait == trait_key & rank <= top_n]
  n_groups <- uniqueN(d$group)
  cnt <- x[, .(n_groups_in = uniqueN(group)), by = gene][
    , .(n_genes = .N), by = n_groups_in][order(n_groups_in)]

  ggplot(cnt, aes(factor(n_groups_in), n_genes)) +
    geom_col() +
    labs(x = sprintf("number of groups where the gene is a top-%d driver (of %d)",
                     top_n, n_groups),
         y = "number of genes",
         title = sprintf("%s: driver-gene specificity", trait_key),
         subtitle = "Bars at 1 are context-specific drivers; bars at the right are shared across the tissue")
}


#' Diagnostic: for each trait x group, the mean GWAS-weight percentile of the
#' top-N drivers. ~98 means the drivers are simply the highest-weight genes;
#' lower values mean mid-weight genes are tracking risk in that context.
weight_percentile <- function(d, top_n = TOP_N) {
  wp <- d[, {
    pct <- 100 * frank(magma_weight) / .N
    .(weight_pctile = mean(pct[rank <= top_n]), n_cells = first(n_cells))
  }, by = .(trait, group)]

  ggplot(wp, aes(group, trait, fill = weight_pctile)) +
    geom_tile() +
    scale_fill_viridis_c(name = "mean weight\npercentile") +
    labs(x = NULL, y = NULL,
         title = sprintf("Are top-%d drivers just the highest-weight genes?", top_n),
         subtitle = "~98 = weight-driven (adds little over the GWAS ranking); lower = context-specific genes tracking risk") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))
}


#' Load the per-domain risk-gene EXPRESSION table (09_risk_gene_drivers.py).
#' Distinct from load_drivers(): this is expression LEVEL per group, not
#' correlation with the scDRS score.
load_expression <- function(grouping = c("Spatial_Domain", "labels")) {
  grouping <- match.arg(grouping)
  e <- fread(file.path(DRIVER_DIR,
                       sprintf("risk_gene_expression_%s.tsv.gz", grouping)))
  setattr(e, "grouping", grouping)
  e[]
}


#' Load the scDRS group-level association statistics (07/08 outputs).
#' Used to annotate the expression plots with each domain's trait-level
#' association, so gene expression can be read against whether the domain is
#' associated with the trait at all.
load_group_stats <- function(grouping = c("Spatial_Domain", "labels")) {
  grouping <- match.arg(grouping)
  fread(file.path(BASE, "processed-data", "16_scDRS", "group_analysis",
                  sprintf("scdrs_group_stats_%s.tsv", grouping)))
}


#' Build the column annotation strip: one tile per group, coloured by the
#' scDRS association MC z-score for this trait.
#'
#' `group_levels` is passed in so the strip and the main panel share an
#' identical x axis -- otherwise the two would order groups independently and
#' the columns would silently fail to line up.
assoc_strip <- function(stats, trait_key, group_levels) {
  s <- stats[trait == trait_key & group %in% group_levels]
  # Keep groups with no stats row as explicit NA columns rather than dropping
  # them, so the strip stays aligned with the panel below.
  s <- s[data.table(group = group_levels), on = "group"]
  s[, group := factor(group, levels = group_levels)]

  lim <- max(abs(s$assoc_mcz), na.rm = TRUE)
  ggplot(s, aes(group, 1, fill = assoc_mcz)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim),
                         na.value = "grey70", name = "scDRS assoc.\nMC z") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.margin = margin(2, 2, 0, 2))
}


#' Stack an association strip on top of an expression plot.
#'
#' Returns the plot unchanged (with a warning) if patchwork is missing, so the
#' expression figures still work without it.
add_assoc_strip <- function(p, stats, trait_key, group_levels) {
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    warning("patchwork not installed; returning plot without the association strip.")
    return(p)
  }
  strip <- assoc_strip(stats, trait_key, group_levels) +
    theme(legend.position = "left")
  patchwork::wrap_plots(
    strip, p + theme(plot.title = element_blank(), plot.subtitle = element_blank()),
    ncol = 1, heights = c(1, 20)
  ) +
    patchwork::plot_annotation(
      title = sprintf("%s: risk-gene expression by domain", trait_key),
      subtitle = "Top strip: domain-level scDRS association. Below: expression z-scored across domains within each gene."
    )
}


#' Where are a trait's top-weighted risk genes expressed?
#'
#' Genes are chosen by MAGMA weight (a genetics-only criterion, independent of
#' this dataset, so the selection is not circular with what is plotted).
#' Fill is expression z-scored ACROSS groups within a gene, so a red tile means
#' "high for this gene relative to other domains" -- not "abundant". Raw
#' abundance differences between genes are deliberately removed; without this
#' the figure mostly recapitulates cell-type composition.
expression_heatmap <- function(e, trait_key, n_genes = 30, stats = NULL) {
  x <- e[trait == trait_key]
  if (nrow(x) == 0) stop("no rows for trait ", trait_key)

  keep <- unique(x[, .(gene, magma_weight)])[order(-magma_weight)][
    seq_len(min(n_genes, .N)), gene]
  x <- x[gene %in% keep]

  # Order genes by weight so the strongest GWAS evidence sits at the top.
  gene_ord <- unique(x[, .(gene, magma_weight)])[order(magma_weight), gene]
  x[, gene := factor(gene, levels = gene_ord)]
  group_levels <- sort(unique(x$group))
  x[, group := factor(group, levels = group_levels)]

  lim <- max(abs(x$expr_z), na.rm = TRUE)
  p <- ggplot(x, aes(group, gene, fill = expr_z)) +
    geom_tile() +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                         midpoint = 0, limits = c(-lim, lim),
                         name = "expression\nz (across\ngroups)") +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = sprintf("%s: where are the top-%d MAGMA risk genes expressed?",
                         trait_key, length(keep)),
         subtitle = "Red = high for that gene relative to its other domains. Not comparable between genes.") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))

  if (is.null(stats)) p else add_assoc_strip(p, stats, trait_key, group_levels)
}


#' Same gene selection, but showing magnitude and detection together:
#' size = fraction of cells with a nonzero count, colour = z-scored mean.
#'
#' Use this to catch the main failure mode of the heatmap above -- a large
#' z-score in a gene detected in very few cells is noise, and shows up here as
#' a strongly coloured but tiny point.
#'
#' NOTE the units: frac_expressing is a FRACTION on 0-1 (0.2 = 20% of cells).
#' The legend is labelled in percent for readability, so the breaks are
#' specified as fractions and only the labels are scaled.
expression_dotplot <- function(e, trait_key, n_genes = 30, stats = NULL) {
  # Tolerate the pre-rename column name so older output still plots.
  if (!"frac_expressing" %in% names(e) && "pct_expressing" %in% names(e)) {
    setnames(e, "pct_expressing", "frac_expressing")
  }
  x <- e[trait == trait_key]
  if (nrow(x) == 0) stop("no rows for trait ", trait_key)

  keep <- unique(x[, .(gene, magma_weight)])[order(-magma_weight)][
    seq_len(min(n_genes, .N)), gene]
  x <- x[gene %in% keep]
  gene_ord <- unique(x[, .(gene, magma_weight)])[order(magma_weight), gene]
  x[, gene := factor(gene, levels = gene_ord)]
  group_levels <- sort(unique(x$group))
  x[, group := factor(group, levels = group_levels)]

  lim <- max(abs(x$expr_z), na.rm = TRUE)
  p <- ggplot(x, aes(group, gene, size = frac_expressing, colour = expr_z)) +
    geom_point() +
    scale_colour_gradient2(low = "#2166AC", mid = "grey90", high = "#B2182B",
                           midpoint = 0, limits = c(-lim, lim),
                           name = "expression z") +
    # Most genes sit near 3% detection, so a linear area scale renders almost
    # every point invisible; sqrt spreads the low end out.
    scale_size_area(name = "% cells\ndetected", max_size = 7,
                    trans = "sqrt",
                    breaks = c(0.01, 0.05, 0.1, 0.25, 0.5, 1),
                    labels = c("1", "5", "10", "25", "50", "100")) +
    scale_x_discrete(drop = FALSE) +
    labs(x = NULL, y = NULL,
         title = sprintf("%s: expression of top-%d MAGMA risk genes by domain",
                         trait_key, length(keep)),
         subtitle = "Small points are detected in few cells -- treat their z-scores as noise.") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.text.y = element_text(size = 7))

  if (is.null(stats)) p else add_assoc_strip(p, stats, trait_key, group_levels)
}


#' Save a plot to plots/16_scDRS/drivers/ as both PDF and PNG.
save_plot <- function(p, name, width = 20, height = 9) {
  ggsave(file.path(PLOT_DIR, paste0(name, ".pdf")), p, width = width, height = height)
  ggsave(file.path(PLOT_DIR, paste0(name, ".png")), p, width = width, height = height, dpi = 200)
  message("wrote ", file.path(PLOT_DIR, paste0(name, ".{pdf,png}")))
  invisible(p)
}

#' Generate and save every per-trait figure for one grouping.
#'
#' Note the axis convention: all five plot functions take a TRAIT and place the
#' groups (spatial domains) along the x axis. So one figure already covers every
#' domain, and the loop below is over traits, not domains.
#'
#' The expression figures are skipped with a warning if 09 has not been rerun
#' since the expression output was added -- the driver figures do not depend on
#' it and should still be produced.
generate_all <- function(grouping = "Spatial_Domain", n_genes = 30) {
  d <- load_drivers(grouping)
  stats <- load_group_stats(grouping)

  expr_file <- file.path(DRIVER_DIR,
                         sprintf("risk_gene_expression_%s.tsv.gz", grouping))
  has_expr <- file.exists(expr_file)
  if (!has_expr) {
    warning("Missing ", expr_file,
            "\n  Rerun 09_risk_gene_drivers.sh to produce the expression ",
            "figures. Skipping them for now.", call. = FALSE)
  }
  e <- if (has_expr) load_expression(grouping) else NULL

  traits <- sort(unique(d$trait))
  message(sprintf("Generating figures for %d traits (%s).",
                  length(traits), grouping))

  # Not per-trait: one summary across the whole trait x group table.
  save_plot(weight_percentile(d), sprintf("weight_percentile_%s", grouping))

  for (tr in traits) {
    message("  ", tr)
    # One bad trait should not abort the whole batch, so each figure is
    # attempted independently and failures are reported and skipped.
    tryCatch(
      save_plot(driver_scatter(d, tr), sprintf("driver_scatter_%s_%s", grouping, tr)),
      error = function(err) message("    scatter failed: ", conditionMessage(err))
    )
    tryCatch(
      save_plot(driver_heatmap(d, tr), sprintf("driver_heatmap_%s_%s", grouping, tr),
                height = 12),
      error = function(err) message("    heatmap failed: ", conditionMessage(err))
    )
    tryCatch(
      save_plot(driver_specificity(d, tr),
                sprintf("driver_specificity_%s_%s", grouping, tr),
                width = 10, height = 7),
      error = function(err) message("    specificity failed: ", conditionMessage(err))
    )

    if (!has_expr) next
    tryCatch(
      save_plot(expression_heatmap(e, tr, n_genes = n_genes, stats = stats),
                sprintf("expression_heatmap_%s_%s", grouping, tr),
                width = 12, height = 10),
      error = function(err) message("    expr heatmap failed: ", conditionMessage(err))
    )
    tryCatch(
      save_plot(expression_dotplot(e, tr, n_genes = n_genes, stats = stats),
                sprintf("expression_dotplot_%s_%s", grouping, tr),
                width = 12, height = 10),
      error = function(err) message("    expr dotplot failed: ", conditionMessage(err))
    )
  }

  message("Done. Figures in ", PLOT_DIR)
  invisible(TRUE)
}


# Batch entry point. Guarded so that source()-ing this file interactively only
# defines the functions -- previously the trailing calls fired on every source,
# silently rewriting figures as a side effect of loading the helpers.
if (!interactive()) {
  generate_all("Spatial_Domain")
  print(sessioninfo::session_info())
} else {
  message("Loaded. Start with:  dr <- load_drivers(\"Spatial_Domain\")")
}