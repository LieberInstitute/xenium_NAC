# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#   Goal: Perform QC based on bin-level object built in 01 script
#. Performed at 8um bin level
#   Outputs:
#     - QC plots (spatial + distributional) per sample
#     - A table of flagged bins (for optional downstream filtering)
#     - Updated SPE with QC columns added

library(SpatialExperiment)
library(sessioninfo)
library(spatialLIBD)
library(patchwork)
library(tidyverse)
library(escheR)
library(scuttle)
library(scater)
library(here)

################################################################################
#   Setup
################################################################################

spe <- readRDS(
    here("processed-data", "HD_Full_Analysis", "SPEs", "spe_norm.Rds")
)

plot_dir <- here("plots", "HD_Full_Analysis", "bin_QC")
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

sample_ids <- unique(spe$sample_id)

################################################################################
#   Uniquify rownames to gene symbols
################################################################################

message(Sys.time(), " | Setting rownames to gene symbols via uniquifyFeatureNames")

rownames(spe) <- scuttle::uniquifyFeatureNames(
    rowData(spe)$gene_id,
    rowData(spe)$gene_name
)

################################################################################
#   Compute QC metrics
################################################################################

message(Sys.time(), " | Computing bin-level QC metrics with addPerCellQCMetrics")

#   Identify mitochondrial genes
is_mito <- grepl("^MT-", rowData(spe)$gene_name, ignore.case = TRUE)

#   If gene_name isn't available, fall back to seqnames if present
if (!any(is_mito) && "seqnames" %in% names(rowData(spe))) {
    is_mito <- seqnames(spe) == "chrM"
}

message(sprintf("  Found %d mitochondrial genes out of %d total", sum(is_mito), nrow(spe)))

#   Compute all QC metrics in one call
#   Adds to colData: sum, detected, subsets_mito_sum, subsets_mito_detected,
#   subsets_mito_percent, total
#   NOTE: Ribosomal protein genes (RPL*/RPS*) are absent from this dataset
#   (likely excluded by the probe panel or filtered due to zero counts), so
#   ribosomal QC is omitted.
spe <- addPerCellQCMetrics(spe, subsets = list(
    mito = is_mito
))

#   Log-transformed versions for spatial visualization
spe$log10_sum <- log10(spe$sum + 1)
spe$log10_detected <- log10(spe$detected + 1)

################################################################################
#   Per-sample distributional QC plots
################################################################################

message(Sys.time(), " | Generating per-sample distribution plots")

qc_df <- colData(spe) |>
    as.data.frame() |>
    select(sample_id, sum, detected, subsets_mito_percent)

#   Compute per-sample 3x MAD thresholds for the histogram vlines
#   UMI and genes use lower-tail outliers on the log scale;
#   mito uses upper-tail outliers on the original scale
mad_thresholds <- qc_df |>
    group_by(sample_id) |>
    summarise(
        #   Lower threshold for log(sum): median - 3*MAD on log scale, then back-transform
        sum_cutoff = exp(
            median(log(sum)) - 3 * mad(log(sum))
        ),
        #   Lower threshold for log(detected)
        detected_cutoff = exp(
            median(log(detected)) - 3 * mad(log(detected))
        ),
        #   Upper threshold for mito percent (not log-transformed)
        mito_cutoff = median(subsets_mito_percent, na.rm = TRUE) +
            3 * mad(subsets_mito_percent, na.rm = TRUE),
        .groups = "drop"
    )

message("Per-sample 3x MAD thresholds:")
print(mad_thresholds)

#--- Histogram: UMI counts per bin (log10 scale) ---
p_umi_hist <- ggplot(qc_df, aes(x = sum + 1)) +
    geom_histogram(bins = 100, fill = "steelblue", color = "black", linewidth = 0.1) +
    geom_vline(
        data = mad_thresholds, aes(xintercept = sum_cutoff),
        color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    scale_x_log10() +
    facet_wrap(~ sample_id, scales = "free_y") +
    labs(
        title = "UMI counts per bin (red = 3 MAD lower cutoff)",
        x = "Total UMI (log10)", y = "Number of bins"
    ) +
    theme_bw()

ggsave(
    file.path(plot_dir, "hist_sum_umi.png"), p_umi_hist,
    width = 10, height = 8, dpi = 150
)

#--- Histogram: Genes detected per bin ---
p_gene_hist <- ggplot(qc_df, aes(x = detected + 1)) +
    geom_histogram(bins = 100, fill = "darkgreen", color = "black", linewidth = 0.1) +
    geom_vline(
        data = mad_thresholds, aes(xintercept = detected_cutoff),
        color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    scale_x_log10() +
    facet_wrap(~ sample_id, scales = "free_y") +
    labs(
        title = "Genes detected per bin (red = 3 MAD lower cutoff)",
        x = "Genes detected (log10)", y = "Number of bins"
    ) +
    theme_bw()

ggsave(
    file.path(plot_dir, "hist_sum_gene.png"), p_gene_hist,
    width = 10, height = 8, dpi = 150
)

#--- Histogram: Mito percent ---
p_mito_hist <- ggplot(qc_df, aes(x = subsets_mito_percent)) +
    geom_histogram(bins = 100, fill = "firebrick", color = "black", linewidth = 0.1) +
    geom_vline(
        data = mad_thresholds, aes(xintercept = mito_cutoff),
        color = "red", linetype = "dashed", linewidth = 0.7
    ) +
    facet_wrap(~ sample_id, scales = "free_y") +
    labs(
        title = "Mitochondrial percent per bin (red = 3 MAD upper cutoff)",
        x = "Mito UMI %", y = "Number of bins"
    ) +
    theme_bw()

ggsave(
    file.path(plot_dir, "hist_mito_percent.png"), p_mito_hist,
    width = 10, height = 8, dpi = 150
)

#--- Scatter: genes vs UMI, colored by mito percent ---
p_scatter <- ggplot(
    qc_df |> slice_sample(n = min(nrow(qc_df), 200000)),
    aes(x = sum, y = detected, color = subsets_mito_percent)
) +
    geom_point(size = 0.2, alpha = 0.3) +
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
summary_df <- qc_df |>
    group_by(sample_id) |>
    summarise(
        n_bins = n(),
        median_umi = median(sum),
        median_genes = median(detected),
        median_mito_pct = median(subsets_mito_percent, na.rm = TRUE),
        pct_below_5_umi = mean(sum < 5) * 100,
        pct_mito_above_20 = mean(subsets_mito_percent > 20, na.rm = TRUE) * 100,
        .groups = "drop"
    )

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
#   These are the most important plots: they reveal tissue-level artifacts
#   (stripes, bubbles, folds, permeabilization failures) that distributional
#   plots cannot catch.
#
#   NOTE on spot_size / point_size: Visium HD 8um bins are much smaller and
#   denser than standard Visium spots. We use very small point sizes (0.1-0.3)
#   to avoid over-plotting. You may need to adjust depending on your tissue
#   coverage and display resolution.

message(Sys.time(), " | Generating spatial QC plots with escheR")

#   escheR needs variables in colData, which we've already added above.
#   We trim colData before piping into make_escheR to speed up plotting
#   (per escheR docs recommendation).

escher_spot_size <- 0.2
escher_point_size <- 0.2

for (sid in sample_ids) {
    message(sprintf("  Plotting sample: %s", sid))

    spe_sub <- spe[, spe$sample_id == sid]

    #   Trim colData to only QC-relevant columns to speed up escheR plotting
    keep_cols <- c(
        "sample_id", "sum", "detected", "log10_sum", "log10_detected",
        "subsets_mito_sum", "subsets_mito_percent"
    )
    keep_cols <- intersect(keep_cols, colnames(colData(spe_sub)))
    colData(spe_sub) <- colData(spe_sub)[, keep_cols]

    #--- Spatial: log10(UMI counts) ---
    p1 <- make_escheR(spe_sub, spot_size = escher_spot_size) |>
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
    p2 <- make_escheR(spe_sub, spot_size = escher_spot_size) |>
        add_fill(var = "subsets_mito_percent", point_size = escher_point_size) +
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
    p3 <- make_escheR(spe_sub, spot_size = escher_spot_size) |>
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
}

################################################################################
#   Flag outlier bins using MAD-based detection
################################################################################

message(Sys.time(), " | Flagging outlier bins (per sample, MAD-based)")

#   Run outlier detection per sample to account for sample-level differences
spe$qc_low_umi <- FALSE
spe$qc_low_gene <- FALSE
spe$qc_high_mito <- FALSE
spe$qc_discard <- FALSE

for (sid in sample_ids) {
    idx <- spe$sample_id == sid

    low_umi <- isOutlier(
        spe$sum[idx], type = "lower", log = TRUE, nmads = 3
    )
    low_gene <- isOutlier(
        spe$detected[idx], type = "lower", log = TRUE, nmads = 3
    )
    high_mito <- isOutlier(
        spe$subsets_mito_percent[idx], type = "higher", nmads = 3
    )

    spe$qc_low_umi[idx] <- low_umi
    spe$qc_low_gene[idx] <- low_gene
    spe$qc_high_mito[idx] <- high_mito
    spe$qc_discard[idx] <- low_umi | low_gene | high_mito

    message(sprintf(
        "  %s: %d/%d bins flagged (%.1f%%) — low_umi: %d, low_gene: %d, high_mito: %d",
        sid, sum(low_umi | low_gene | high_mito), sum(idx),
        mean(low_umi | low_gene | high_mito) * 100,
        sum(low_umi), sum(low_gene), sum(high_mito)
    ))
}

#--- Spatial plot of flagged bins with escheR ---
#   Convert to factor for categorical coloring with add_fill
spe$qc_discard_factor <- factor(
    spe$qc_discard, levels = c(FALSE, TRUE), labels = c("Pass", "Flagged")
)

message(Sys.time(), " | Plotting QC-flagged bins spatially")

for (sid in sample_ids) {
    spe_sub <- spe[, spe$sample_id == sid]

    #   Trim colData for speed
    colData(spe_sub) <- colData(spe_sub)[, c("sample_id", "qc_discard_factor")]

    p_flag <- make_escheR(spe_sub, spot_size = escher_spot_size) |>
        add_fill(var = "qc_discard_factor", point_size = escher_point_size) +
        scale_fill_manual(
            values = c("Pass" = "grey80", "Flagged" = "red"),
            name = "QC Status"
        ) +
        ggtitle(paste(sid, "- QC-flagged bins")) +
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
#   Save updated SPE with QC columns
################################################################################

message(Sys.time(), " | Saving SPE with QC annotations")

spe_qc_path <- here(
    "processed-data", "HD_Full_Analysis", "SPEs", "spe_norm_with_binQC.Rds"
)
saveRDS(spe, spe_qc_path)

message(sprintf("Saved to: %s", spe_qc_path))
message(sprintf(
    "Total bins flagged: %d / %d (%.1f%%)",
    sum(spe$qc_discard), ncol(spe), mean(spe$qc_discard) * 100
))

################################################################################
#   Remove discarded bins and save filtered SPE
################################################################################

message(Sys.time(), " | Filtering out QC-flagged bins")

n_before <- ncol(spe)
spe_filtered <- spe[, !spe$qc_discard]
n_after <- ncol(spe_filtered)

message(sprintf(
    "  Removed %d bins: %d -> %d (%.1f%% retained)",
    n_before - n_after, n_before, n_after,
    (n_after / n_before) * 100
))

#   Per-sample breakdown
for (sid in sample_ids) {
    idx_before <- sum(spe$sample_id == sid)
    idx_after <- sum(spe_filtered$sample_id == sid)
    message(sprintf(
        "  %s: %d -> %d bins (%.1f%% retained)",
        sid, idx_before, idx_after, (idx_after / idx_before) * 100
    ))
}

#   Drop the QC flag columns that are no longer needed
spe_filtered$qc_discard_factor <- NULL

spe_filtered_path <- here(
    "processed-data", "HD_Full_Analysis", "SPEs", "spe_norm_binQC_filtered.Rds"
)
saveRDS(spe_filtered, spe_filtered_path)

message(sprintf("Saved filtered SPE to: %s", spe_filtered_path))


#Save the barcodes that passed filtering
barcode_dir <- here("processed-data", "HD_Full_Analysis", "SPEs")
for (sid in sample_ids) {
    passing_barcodes <- colnames(spe_filtered[, spe_filtered$sample_id == sid])
    write.csv(
        data.frame(barcode = passing_barcodes),
        file.path(barcode_dir, paste0("passing_barcodes_", sid, ".csv")),
        row.names = FALSE
    )
    message(sprintf("  Wrote %d passing barcodes for %s", length(passing_barcodes), sid))
}

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
