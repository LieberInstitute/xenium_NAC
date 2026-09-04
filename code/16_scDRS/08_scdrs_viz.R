# Visualize scDRS results: trait x spatial-domain and trait x cell-type
# heatmaps (ComplexHeatmap), plus spatial maps of per-cell disease scores
# for key traits (escheR).
#
# Usage: sbatch 08_scdrs_viz.sh

library(here)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(SpatialFeatureExperiment)
library(HDF5Array)
library(escheR)
library(ggplot2)

ga_dir <- here("processed-data", "16_scDRS", "group_analysis")
tab_dir <- here("processed-data", "16_scDRS", "tables")
plot_dir <- here("plots", "16_scDRS")

manifest <- fread(here("processed-data", "16_scDRS", "trait_manifest.tsv"))

## ---- Trait subset for this figure set ---------------------------------------
# Restrict the figures to a 21-trait subset: drop the three epilepsy phenotypes
# and the two anthropometric negative controls, leaving T2D as the sole
# negative control.
#
# NOTE: trait_manifest.tsv itself is deliberately NOT edited. It is shared with
# steps 02/04/05/06/07 -- including the --array=1-26 ranges in 04 and 06 and the
# completeness guard in 07_scdrs_downstream_merge.py -- so subsetting the file
# would break a full pipeline rerun. The subset is applied here only.
EXCLUDE_TRAITS <- c("EPI_ALL", "EPI_FOCAL", "EPI_GGE", "HEIGHT", "BMI")
manifest <- manifest[!trait_key %in% EXCLUDE_TRAITS]
message("Plotting ", nrow(manifest), " traits; excluded: ",
        paste(EXCLUDE_TRAITS, collapse = ", "))

## ---- Category colors for column annotation ----------------------------------
cat_order <- c("psychiatric", "substance_use", "cognitive",
               "neurological", "negative_control")
cat_colors <- c(
  psychiatric      = "#E41A1C",
  substance_use    = "#377EB8",
  cognitive        = "#4DAF4A",
  neurological     = "#984EA3",
  negative_control = "#999999"
)

## ---- Group-level heatmaps (ComplexHeatmap) ----------------------------------
plot_group_heatmap <- function(groupby) {
  stats <- fread(file.path(ga_dir, sprintf("scdrs_group_stats_%s.tsv", groupby)))

  # Apply the trait subset BEFORE multiple-testing correction, so FDR is
  # computed across only the traits actually displayed (see EXCLUDE_TRAITS).
  stats <- stats[trait %in% manifest$trait_key]

  # FDR across the displayed trait x group table
  stats[, assoc_fdr := p.adjust(assoc_mcp, method = "BH")]
  stats[, hetero_fdr := p.adjust(hetero_mcp, method = "BH")]
  message(groupby, ": FDR computed across ", nrow(stats),
          " trait x group tests (", uniqueN(stats$trait), " traits)")

  # Significance symbols:
  #   Nominal p (assoc_mcp): +, ++, +++
  #   FDR (assoc_fdr):       *, **, ***
  #   FDR takes precedence over nominal when both apply
  stats[, sig := ""]
  stats[assoc_mcp <= 0.05,  sig := "+"]
  stats[assoc_mcp <= 0.01,  sig := "++"]
  stats[assoc_mcp <= 0.001, sig := "+++"]
  stats[assoc_fdr <= 0.1,   sig := "*"]
  stats[assoc_fdr <= 0.05,  sig := "**"]
  stats[assoc_fdr <= 0.01,  sig := "***"]

  stats <- merge(
    stats, manifest[, .(trait_key, trait_name, category)],
    by.x = "trait", by.y = "trait_key"
  )

  # Save FDR table
  fwrite(stats, file.path(tab_dir, sprintf("scdrs_group_stats_%s_fdr.csv", groupby)))

  # ---- Build the z-score matrix (rows = groups, columns = traits) ----
  zmat <- dcast(stats, group ~ trait, value.var = "assoc_mcz")
  zm <- as.matrix(zmat[, -1])
  rownames(zm) <- zmat$group
  zm[is.na(zm)] <- 0

  # Build matching significance matrix
  sigmat <- dcast(stats, group ~ trait, value.var = "sig")
  sm <- as.matrix(sigmat[, -1])
  rownames(sm) <- sigmat$group

  # ---- Order columns by category, then alphabetically within ----
  manifest[, category := factor(category, levels = cat_order)]
  setorder(manifest, category, trait_key)
  trait_ord <- manifest$trait_key[manifest$trait_key %in% colnames(zm)]
  zm <- zm[, trait_ord]
  sm <- sm[, trait_ord]

  # Column split factor (category per trait)
  col_category <- manifest[match(trait_ord, trait_key), category]

  # ---- Column annotation: disease category color bar ----
  ha_col <- HeatmapAnnotation(
    Category = as.character(col_category),
    col = list(Category = cat_colors),
    annotation_name_side = "left"
  )

  # ---- Row clustering by hierarchical clustering of z-scores ----
  row_dend <- hclust(dist(zm))

  # ---- Color scale (diverging blue-white-red) ----
  max_abs <- max(abs(zm), na.rm = TRUE)
  col_fun <- colorRamp2(c(-max_abs, 0, max_abs), c("#2166AC", "white", "#B2182B"))

  # ---- Build heatmap ----
  hm <- Heatmap(
    zm,
    name = "scDRS assoc.\nMC z-score",
    col = col_fun,
    cluster_rows = row_dend,
    cluster_columns = FALSE,
    column_split = col_category,
    column_title_rot = 0,
    top_annotation = ha_col,
    row_names_side = "left",
    column_names_rot = 45,
    column_names_gp = gpar(fontsize = 9),
    row_names_gp = gpar(fontsize = 9),
    rect_gp = gpar(col = "grey80", lwd = 0.5),
    cell_fun = function(j, i, x, y, width, height, fill) {
      txt <- sm[i, j]
      if (!is.na(txt) && nchar(txt) > 0) {
        grid.text(txt, x, y, gp = gpar(fontsize = 7, col = "black"))
      }
    },
    column_title_gp = gpar(fontsize = 10, fontface = "bold"),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 9),
      labels_gp = gpar(fontsize = 8)
    )
  )

  # ---- Save ----
  n_rows <- nrow(zm)
  n_cols <- ncol(zm)
  pdf(
    file.path(plot_dir, sprintf("scdrs_heatmap_%s.pdf", groupby)),
    width = 4 + 0.3 * n_cols,
    height = 3 + 0.3 * n_rows
  )
  # Legend key for significance symbols
  sig_legend <- Legend(
    title = "Significance",
    labels = c("*** FDR <= 0.01", "** FDR <= 0.05", "* FDR <= 0.1",
               "+++ p <= 0.001", "++ p <= 0.01", "+ p <= 0.05"),
    type = "points",
    pch = NA,
    labels_gp = gpar(fontsize = 8),
    title_gp = gpar(fontsize = 9, fontface = "bold")
  )

  draw(hm,
       column_title = sprintf("scDRS trait association by %s", groupby),
       column_title_gp = gpar(fontsize = 13, fontface = "bold"),
       heatmap_legend_side = "right",
       annotation_legend_list = list(sig_legend))
  dev.off()

  invisible(stats)
}

stats_domain <- plot_group_heatmap("Spatial_Domain")
stats_labels <- plot_group_heatmap("labels")

## ---- Spatial maps of per-cell scores (escheR) -------------------------------
key_traits <- manifest$trait_key

# Load SpatialFeatureExperiment object
sfe <- loadHDF5SummarizedExperiment(
  here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
)

# Read normalized scDRS scores and attach to the object
# Scores were computed on a filtered subset (Hypo/Excitatory excluded),
# so subset the SFE to the scored cells.
scores <- fread(file.path(ga_dir, "norm_scores_all_traits.tsv.gz"))
setnames(scores, 1, "cell_id")
sfe <- sfe[, colnames(sfe) %in% scores$cell_id]
scores <- scores[match(colnames(sfe), cell_id)]
stopifnot(identical(scores$cell_id, colnames(sfe)))

samples <- unique(sfe$Sample)

# sfe$Sample holds the full absolute path to each sample's spaceranger .h5
# file, e.g. ".../01_spaceranger/H1-XKYDCP3_A1/outs/.../matrix.h5". Using that
# directly in file.path() recreated the whole absolute tree as nested
# directories (a stray "dcs05/lieber/..." folder under plots/). Extract the
# real sample ID: the directory name immediately under "01_spaceranger".
sample_id <- function(s) basename(dirname(dirname(dirname(as.character(s)))))

# Map each raw Sample value to its clean ID once.
sample_ids <- setNames(sample_id(samples), samples)
sample_dir <- function(s) file.path(plot_dir, sample_ids[[as.character(s)]])

# Create per-sample directories
for (s in samples) {
  dir.create(sample_dir(s), showWarnings = FALSE, recursive = TRUE)
}

# Plot every trait we tested. key_traits comes from the manifest; warn loudly
# if any manifest trait lacks a score column so nothing silently drops out.
plot_traits <- intersect(key_traits, names(scores))
missing_traits <- setdiff(key_traits, names(scores))
if (length(missing_traits) > 0) {
  warning(sprintf(
    "No score column for %d manifest trait(s): %s",
    length(missing_traits), paste(missing_traits, collapse = ", ")
  ))
}
message(sprintf("Plotting spatial maps for %d traits x %d samples.",
                length(plot_traits), length(samples)))

# Sequential white -> red -> black scale. We only care about cells that harbor
# genetic risk (high scores), so scores <= 0 read as white and darkening
# indicates increasing risk. Negatives are squished onto the white floor.
risk_colors <- c("white", "#B2182B", "black")

for (trait in plot_traits) {
  message(sprintf("Plotting spatial maps for %s ...", trait))
  sfe$norm_score <- scores[[trait]]

  # One scale per trait, computed across all samples of that trait. This makes
  # the eight sample maps for a given disease directly comparable to each other
  # and gives every disease the full colour range (maximising within-disease
  # contrast). Colours are therefore NOT comparable BETWEEN diseases -- the
  # legend on each panel reports its own range.
  trait_max <- max(scores[[trait]], na.rm = TRUE)
  message(sprintf("  %s: colour scale [0, %.3f]", trait, trait_max))

  for (s in samples) {
    sfe_sub <- sfe[, sfe$Sample == s]

    p <- make_escheR(sfe_sub) |>
      add_fill("norm_score") +
      scale_fill_gradientn(
        colours = risk_colors,
        limits = c(0, trait_max),
        oob = scales::squish,
        name = "scDRS\nnorm. score"
      ) +
      ggtitle(sprintf("%s — %s", trait, sample_ids[[as.character(s)]])) +
      theme(plot.title = element_text(hjust = 0.5))

    ggsave(
      filename = sprintf("scdrs_spatial_%s.png", trait),
      path = sample_dir(s),
      plot = p,
      width = 16, height = 16, dpi = 300
    )
  }
}

sessioninfo::session_info()
