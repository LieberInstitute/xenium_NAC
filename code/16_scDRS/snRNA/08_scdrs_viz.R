# Visualize snRNA-seq scDRS results: trait x cell-type heatmap (ComplexHeatmap).
#
# No spatial plots — snRNA-seq has no spatial coordinates.
#
# Usage: sbatch 08_scdrs_viz.sh

library(here)
library(data.table)
library(ComplexHeatmap)
library(circlize)

ga_dir <- here("processed-data", "16_scDRS", "snRNA", "group_analysis")
tab_dir <- here("processed-data", "16_scDRS", "snRNA", "tables")
plot_dir <- here("plots", "16_scDRS", "snRNA")

manifest <- fread(here("processed-data", "16_scDRS", "trait_manifest.tsv"))

## ---- Trait subset for this figure set ---------------------------------------
# Mirrors the Visium-HD figures in ../08_scdrs_viz.R: drop the three epilepsy
# phenotypes and the two anthropometric negative controls, leaving T2D as the
# sole negative control. Keeping this list identical across the two scripts is
# what makes the snRNA and Visium-HD heatmaps comparable -- the BH correction
# below is computed over the displayed table, so a different trait subset would
# silently shift every FDR value.
#
# NOTE: trait_manifest.tsv itself is deliberately NOT edited. It is shared with
# steps 02/04/05/06/07 -- including the --array=1-26 ranges and the completeness
# guard in 07_scdrs_downstream_merge.py -- so subsetting the file would break a
# full pipeline rerun. The subset is applied here only.
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

## ---- Heatmap: trait x CellType.Final ----------------------------------------
stats <- fread(file.path(ga_dir, "scdrs_group_stats_CellType.Final.tsv"))

# Apply the trait subset BEFORE multiple-testing correction, so FDR is computed
# across only the traits actually displayed (see EXCLUDE_TRAITS).
stats <- stats[trait %in% manifest$trait_key]

# FDR across the displayed trait x group table
stats[, assoc_fdr := p.adjust(assoc_mcp, method = "BH")]
stats[, hetero_fdr := p.adjust(hetero_mcp, method = "BH")]
message("CellType.Final: FDR computed across ", nrow(stats),
        " trait x group tests (", uniqueN(stats$trait), " traits)")

# Significance symbols
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

fwrite(stats, file.path(tab_dir, "scdrs_group_stats_CellType.Final_fdr.csv"))

# ---- Build z-score and significance matrices ----
zmat <- dcast(stats, group ~ trait, value.var = "assoc_mcz")
zm <- as.matrix(zmat[, -1])
rownames(zm) <- zmat$group
zm[is.na(zm)] <- 0

sigmat <- dcast(stats, group ~ trait, value.var = "sig")
sm <- as.matrix(sigmat[, -1])
rownames(sm) <- sigmat$group

# ---- Order columns by category, then alphabetically ----
manifest[, category := factor(category, levels = cat_order)]
setorder(manifest, category, trait_key)
trait_ord <- manifest$trait_key[manifest$trait_key %in% colnames(zm)]
zm <- zm[, trait_ord]
sm <- sm[, trait_ord]

col_category <- manifest[match(trait_ord, trait_key), category]

# ---- Column annotation ----
ha_col <- HeatmapAnnotation(
  Category = as.character(col_category),
  col = list(Category = cat_colors),
  annotation_name_side = "left"
)

# ---- Row clustering ----
row_dend <- hclust(dist(zm))

# ---- Color scale ----
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

n_rows <- nrow(zm)
n_cols <- ncol(zm)

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

pdf(
  file.path(plot_dir, "scdrs_heatmap_CellType.Final.pdf"),
  width = 4 + 0.3 * n_cols,
  height = 3 + 0.3 * n_rows
)
draw(hm,
     column_title = "scDRS trait association by CellType.Final (snRNA-seq)",
     column_title_gp = gpar(fontsize = 13, fontface = "bold"),
     heatmap_legend_side = "right",
     annotation_legend_list = list(sig_legend))
dev.off()

sessioninfo::session_info()
