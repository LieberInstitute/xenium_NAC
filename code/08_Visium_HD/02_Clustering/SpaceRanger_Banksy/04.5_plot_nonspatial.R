# Goal: For one BANKSY non-spatial clustering result (a single --res), plot the
#       clusters on each tissue section and build a marker-gene ComplexHeatmap.
library(SummarizedExperiment)
library(SpatialExperiment)
library(ComplexHeatmap)
library(HDF5Array)
library(optparse)
library(circlize)
library(scuttle)
library(ggplot2)
library(escheR)
library(scran)
library(here)

# ---------------------------------------------------------------------------
# 1. Parse resolution (same convention as 04_non-spatial_cluster.R)
# ---------------------------------------------------------------------------
option_list <- list(
  make_option("--res", type = "double", help = "Resolution that was clustered")
)
opt <- parse_args(OptionParser(option_list = option_list))
if (is.null(opt$res)) stop("Error: --res must be provided")
res <- opt$res
message(sprintf("Plotting non-spatial BANKSY result for resolution = %.2f", res))

# ---------------------------------------------------------------------------
# 2. Output directories (the "new directory")
# ---------------------------------------------------------------------------
out_base <- here("plots", "HD_Full_Analysis", "Banksy_sr", "Cell_Level",
                 "NonSpatial_Resolution_Sweep")
spatial_dir <- file.path(out_base, "spatial")
heatmap_dir <- file.path(out_base, "heatmaps")
dir.create(spatial_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------------------------------------------------------
# 3. Load full SFE + this resolution's colData, then attach the cluster column
# ---------------------------------------------------------------------------
message(paste0("Loading SFE - ", Sys.time()))

sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)

sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

cd_path <- here("processed-data", "HD_Full_Analysis", "sr_banksy", "NonSpatial_colData",
                paste0("sr_nonspatial_banksy_res_", res, ".Rds"))
if (!file.exists(cd_path)) stop("Cluster colData not found: ", cd_path)
cd_saved <- readRDS(cd_path)

# Robustly identify the clustering column for this resolution.
clust_cols <- grep("^clust", colnames(cd_saved), value = TRUE)
if (length(clust_cols) == 0) stop("No 'clust*' column found in ", cd_path)
res_in_name <- suppressWarnings(as.numeric(sub(".*res", "", clust_cols)))
res_col <- clust_cols[which(abs(res_in_name - res) < 1e-9)]
if (length(res_col) != 1) res_col <- tail(clust_cols, 1)  # fallback: most recent
message("Using cluster column: ", res_col)

# Align cells (same object, so colnames match) and attach.
stopifnot(all(colnames(sfe) %in% rownames(cd_saved)))
sfe[[res_col]] <- as.character(cd_saved[colnames(sfe), res_col])

# Order cluster levels numerically when possible, else alphabetically.
lev_raw <- unique(sfe[[res_col]])
lev_num <- suppressWarnings(as.numeric(lev_raw))
lev <- if (!any(is.na(lev_num))) lev_raw[order(lev_num)] else sort(lev_raw)
sfe[[res_col]] <- factor(sfe[[res_col]], levels = lev)
n_clust <- length(lev)
message(sprintf("Resolution %.2f -> %d clusters", res, n_clust))

samples <- as.character(unique(sfe$sample_id))

# ---------------------------------------------------------------------------
# 4. Cluster color palette (shared across all sections so colors are consistent)
# ---------------------------------------------------------------------------
make_palette <- function(n) {
  if (requireNamespace("Polychrome", quietly = TRUE)) {
    seeds <- c("#1F77B4", "#FF7F0E", "#2CA02C")
    pal <- suppressWarnings(
      Polychrome::createPalette(max(n + length(seeds), 3), seeds)
    )
    pal <- as.vector(pal)[-seq_along(seeds)]
  } else {
    pal <- grDevices::hcl.colors(n, palette = "Dynamic")
  }
  pal <- pal[seq_len(n)]
  setNames(pal, lev)
}
cluster_cols <- make_palette(n_clust)

# ---------------------------------------------------------------------------
# 5. Spatial plots: one panel per tissue section, one PNG for the resolution
# ---------------------------------------------------------------------------
message(paste0("Building spatial plots - ", Sys.time()))

plot_one_sample <- function(s) {
  sub <- sfe[, sfe$sample_id == s]
  make_escheR(sub) |>
    add_fill(var = res_col, point_size = 0.15) +
    scale_fill_manual(values = cluster_cols, drop = FALSE, name = "cluster") +
    guides(fill = guide_legend(override.aes = list(size = 2.5), ncol = 1)) +
    ggtitle(s) +
    theme(plot.title  = element_text(size = 9, hjust = 0.5),
          axis.text   = element_blank(),
          axis.ticks  = element_blank(),
          axis.title  = element_blank())
}

plots <- lapply(samples, plot_one_sample)

ncol_grid <- ceiling(sqrt(length(samples)))
nrow_grid <- ceiling(length(samples) / ncol_grid)

spatial_png <- file.path(spatial_dir, paste0("spatial_clusters_res_", res, ".png"))

if (requireNamespace("patchwork", quietly = TRUE)) {
  combined <- patchwork::wrap_plots(plots, ncol = ncol_grid, guides = "collect") +
    patchwork::plot_annotation(
      title = sprintf("Non-spatial BANKSY clusters  |  resolution %.1f  |  %d clusters",
                      res, n_clust))
  ggsave(spatial_png, combined,
         width = ncol_grid * 3.4 + 1.6, height = nrow_grid * 3.2 + 0.6,
         dpi = 200, limitsize = FALSE)
} else {
  # Fallback without patchwork: one PNG per section.
  for (i in seq_along(samples)) {
    f <- file.path(spatial_dir, paste0("spatial_clusters_res_", res, "_", samples[i], ".png"))
    ggsave(f, plots[[i]], width = 6, height = 5.5, dpi = 200, limitsize = FALSE)
  }
}
message("Wrote spatial plot(s) to: ", spatial_dir)

# ---------------------------------------------------------------------------
# 6. Marker-gene ComplexHeatmap (adapted from the annotated-domain heatmap)
#    Rows = numbered clusters at this resolution; columns = marker genes.
# ---------------------------------------------------------------------------
message(paste0("Building marker heatmap - ", Sys.time()))

if (!"logcounts" %in% assayNames(sfe)) {
  stop("Assay 'logcounts' not found in SFE; available: ",
       paste(assayNames(sfe), collapse = ", "))
}

splitit <- function(x) split(seq(along = x), x)

# Marker panel (unchanged from your code).
markers_all <- c("SNAP25","GAD1","GAD2","SLC32A1",            # GABA neurons
                 "PPP1R1B","FOXP2","BCL11B",                  # Broad neurons
                 "NPY","CORT","SST","CHODL",                  # Inhibitory subset
                 "DRD1","RELN","TAC1","PDYN",                 # D1_MSN
                 "DRD2","ADORA2A","PENK","GPR6",              # D2_MSN
                 "RXFP1","TSHZ1", "OPRM1",                            # D1/D2 islands
                 "SEMA5B","TRHDE","GABRQ","VWC2L","CPNE4",    # D1_Island_A
                 "SEMA3E","PROK2","NPY1R","RXFP2","KCNH5","VIP", # D1_Island_B
                 "SLC18A3","SLC5A7","CHAT","PVALB","GFRA2","KIT", # CHAT
                 "SLC17A7","TBR1",                            # Excitatory
                 "CLDN5","NR2F2",                             # Endo
                 "C3","P2RY13","MS4A6A",                      # Microglia
                 "GJA1","AQP4","GFAP",                        # Astrocyte
                 "CFAP157","CAPS",                            # Ependymal
                 "OPALIN","MOBP","MOG",                       # Oligo / WM
                 "PDGFRA","VCAN")                             # OPC

# Marker labels reordered to MATCH markers_all (your original had WM/OPC/Ependymal
# at the tail, but the genes run Ependymal -> WM -> OPC).
marker_labels <- c(rep("GABA",4),
                   rep("MSN",3),
                   rep("SST",4),
                   rep("DRD1_MSN",4),
                   rep("DRD2_MSN",4),
                   rep("D1_Islands",14),
                   rep("CHAT_PVALB",6),
                   rep("Excitatory",2),
                   rep("Endo",2),
                   rep("Microglia",3),
                   rep("Astrocyte",3),
                   rep("Ependymal",2),
                   rep("WM",3),
                   rep("OPC",2))
stopifnot(length(markers_all) == length(marker_labels))

# Keep only genes present in the (targeted Xenium) panel, keeping labels aligned.
present <- markers_all %in% rownames(sfe)
if (any(!present)) {
  message("Markers not in panel (dropped): ",
          paste(markers_all[!present], collapse = ", "))
}
markers_use <- markers_all[present]
labels_use  <- factor(marker_labels[present],
                      levels = unique(marker_labels[present]))

# Mean logcounts per cluster, then z-score each gene across clusters
# (clusters as rows, genes as columns) -- same transform as your code.
clust <- sfe[[res_col]]
cell_idx <- splitit(clust)
dat <- as.matrix(assay(sfe, "logcounts")[markers_use, , drop = FALSE])

hm_mat <- scale(
  t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[, i, drop = FALSE])))),
  center = TRUE, scale = TRUE)
hm_mat <- hm_mat[levels(clust), , drop = FALSE]   # numeric cluster order

col_fun <- circlize::colorRamp2(c(min(hm_mat), 0, max(hm_mat)),
                                c("blue", "white", "red"))

# Column annotation: marker category.
col_ha <- ComplexHeatmap::columnAnnotation(
  marker = labels_use,
  show_annotation_name = FALSE,
  show_legend = TRUE)

# Row annotation: number of cells per cluster + cluster color key.
cl_sizes <- as.integer(table(clust)[levels(clust)])
row_ha <- ComplexHeatmap::rowAnnotation(
  cluster = levels(clust),
  n_cells = ComplexHeatmap::anno_barplot(cl_sizes,
                                         gp = gpar(fill = "grey40"),
                                         width = unit(2, "cm")),
  col = list(cluster = cluster_cols),
  show_annotation_name = c(cluster = FALSE, n_cells = TRUE),
  annotation_name_gp = gpar(fontsize = 8),
  show_legend = FALSE)

hm <- ComplexHeatmap::Heatmap(
  matrix = hm_mat,
  name = "centered,\nscaled",
  column_title = sprintf("Marker expression across non-spatial clusters (res %.1f)", res),
  column_title_gp = gpar(fontface = "bold"),
  cluster_rows = TRUE,            # group similar numbered clusters (aids annotation)
  cluster_columns = TRUE,        # keep marker order
  show_row_dend = TRUE,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 7),
  bottom_annotation = col_ha,
  right_annotation = row_ha,
  column_split = labels_use,
  column_gap = unit(1, "mm"),
  row_title = NULL,
  rect_gp = gpar(col = "gray50", lwd = 0.5),
  col = col_fun)

heatmap_pdf <- file.path(heatmap_dir, paste0("marker_heatmap_res_", res, ".pdf"))
pdf(file = heatmap_pdf,
    width = max(16, 0.22 * length(markers_use) + 6),
    height = max(8, 0.28 * n_clust + 4))
ComplexHeatmap::draw(hm)
dev.off()
message("Wrote heatmap to: ", heatmap_pdf)

# ---------------------------------------------------------------------------
# 7. Reproducibility
# ---------------------------------------------------------------------------
message("Done - ", Sys.time())
options(width = 120)
sessioninfo::session_info()
