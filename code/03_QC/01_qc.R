#!/usr/bin/env Rscript
#Goal: Calculate QC metrics and explore. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
#Modified from:
 # https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/f3b87d7a0356e716468b304196f5cd362fd91515/code/analysis/02_xenium_qc/01_threshold_outliers.R
 # https://github.com/LieberInstitute/spatialAmygdala/blob/32108cdc145bb83aa73822c862dc28299a29b47c/code/Xenium/03_quality_control/01_perCellQC.R
 # https://pachterlab.github.io/voyager/articles/vig5_xenium.html#quality-control
 # http://127.0.0.1:29231/library/SpotSweeper/doc/getting_started.html#identifying-local-outliers-using-spotsweeper 

library(SpatialExperiment)
library(scattermore)
library(tidyverse)
#library(Voyager)
library(scuttle)
library(scater)
library(escheR)
library(scran)
library(dplyr)
library(here)
library(SpotSweeper)
library(sessioninfo)
library(ggplot2)

  
#Read in the RDS file from 01_build_spe. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))
spe$Donor[spe$Donor == "Br6426"] <- "Br6436"
#Make the sample column a factor
spe$Sample <- factor(x = spe$Sample, levels = unique(spe$Sample))

# Remove empty cells
empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]           # 20395 cells
spe <- spe[, colSums(counts(spe)) > 0]                          # 4895438 cells

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel

#Per: https://genomics.uci.edu/wp-content/uploads/sites/30/PDF_UCI_GRT_Hub_Xenium_Workshop_20240118-120fec8a1ce5134f.pdf
# "Negative Control Codewords are
# a random subset of codewords,
# with identical properties to gene
# codewords" 
#and
# "By definition, a call made to negative control codeword is an error"
#Unassigned probes help measure noise and off-target activity because they do not match any gene codewords
#Negative control probes are sequences that should not bind anything and therefore measure false positive rates and specificity

#Add QC metrics 
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX))

# # 0) Remove cells that are 100% negatives (do this BEFORE subsetting rows to GEX)
# spe <- spe[, spe$subsets_any_neg_percent < 100]            # 

# 1) Subset to GEX rows (panel genes)
spe <- spe[is_GEX, ]

# 2) Remove empty cells (based on GEX-only counts)
cts_gex <- counts(spe)
spe$sum_gex      <- as.numeric(Matrix::colSums(cts_gex))
spe$detected_gex <- as.numeric(Matrix::colSums(cts_gex > 0))
spe <- spe[, spe$sum_gex > 0]                           # 58 cells removed (all had 0 counts in GEX)

# # 3) Fixed-threshold flags (GEX-only)
# spe$exclude_low_lib_gex  <- spe$sum_gex < 10.          # 137383 cells
# spe$exclude_detected_gex <- spe$detected_gex < 4       # 49872 cells

# 3) Noise fraction from FULL features (already computed pre-subset)
spe$exclude_any_neg <- spe$subsets_any_neg_percent >= 25              # 324 cells

# 4) Adaptive total counts and detected genes outliers (GEX-only metrics), per sample
# spe$sum_gex_3MAD       <- isOutlier(spe$sum_gex,      log=TRUE, batch=as.factor(spe$Sample), nmads=3, type="lower")       # 38830 cells
# spe$detected_gex_3MAD  <- isOutlier(spe$detected_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=3, type="lower")       # 107546 cells
spe$sum_gex_4MAD       <- isOutlier(spe$sum_gex,      log=TRUE, batch=as.factor(spe$Sample), nmads=4, type="lower")       # 7753 cells
spe$detected_gex_4MAD  <- isOutlier(spe$detected_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=4, type="lower")       # 43750 cells
# spe$sum_gex_5MAD       <- isOutlier(spe$sum_gex,      log=TRUE, batch=as.factor(spe$Sample), nmads=5, type="lower")       # 1531 cells
# spe$detected_gex_5MAD  <- isOutlier(spe$detected_gex, log=TRUE, batch=as.factor(spe$Sample), nmads=5, type="lower")       # 20001 cells

# 5) Cell area outliers (after removing GEX-empties)
# spe$cell_area_6MAD <- isOutlier(spe$cell_area, log=FALSE, batch=as.factor(spe$Sample), nmads=6, type="both")              # 6 MAD: 138 cells
# spe$cell_area_5MAD <- isOutlier(spe$cell_area, log=FALSE, batch=as.factor(spe$Sample), nmads=5, type="both")              # 5 MAD: 549 cells
spe$cell_area_4MAD <- isOutlier(spe$cell_area, log=FALSE, batch=as.factor(spe$Sample), nmads=4, type="both")              # 4 MAD: 2157 cells

table(spe$Sample, spe$sum_gex_4MAD)
table(spe$Sample, spe$detected_gex_4MAD)
table(spe$Sample, spe$cell_area_4MAD)

# global outliers flag
spe$global_outliers <- as.logical(spe$exclude_any_neg) |
                       as.logical(spe$detected_gex_4MAD) | 
                       as.logical(spe$sum_gex_4MAD) |
                       as.logical(spe$cell_area_4MAD)

# save outlier tables
processed_outdir <- here("processed-data","Xenium_QC")

cd <- as.data.frame(colData(spe))
df <- cd[ cd$global_outliers , , drop = FALSE ]
write.csv(df, file = file.path(processed_outdir, "global_outliers.csv"), row.names = TRUE)

as.data.frame(colData(spe)) %>%
  group_by(Sample) %>%
  summarize(
    CellsRetained                        = sum(!global_outliers),
    Outliers                             = sum(global_outliers),
    PercentRemoved                       = Outliers/(CellsRetained + Outliers)*100,
    exclude_any_neg                      = sum(exclude_any_neg),
    detected_gex_4MAD                    = sum(detected_gex_4MAD),
    sum_gex_4MAD                         = sum(sum_gex_4MAD),
    cell_area_4MAD                       = sum(cell_area_4MAD)
  ) %>%
  ungroup() %>%
  write.csv(
    file      = here("processed-data","Xenium_QC","global_outliers_summary.csv"),
    row.names = FALSE,
    quote     = FALSE
  )


######################### PLOTTING #########################
# Plotting outliers by sample and pooled
outdir <- here::here("plots", "03_qc")
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# By sample: outliers only, median lines from FULL spe (per sample)
plot_flag_hist_by_sample <- function(flag_col, value_col, label, bins = 50) {
  df_out <- as.data.frame(colData(spe)) %>%
    filter(!is.na(.data[[flag_col]]), .data[[flag_col]]) %>%
    select(Sample, !!rlang::sym(value_col)) %>%
    rename(value = !!rlang::sym(value_col))

  if (nrow(df_out) == 0) return(invisible(NULL))
  n_total <- nrow(df_out)

  med_tbl <- as.data.frame(colData(spe)) %>%
    select(Sample, !!rlang::sym(value_col)) %>%
    rename(value = !!rlang::sym(value_col)) %>%
    group_by(Sample) %>%
    summarise(med_value = median(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot(df_out, aes(x = value)) +
    geom_histogram(bins = bins, alpha = 0.9) +
    geom_vline(data = med_tbl, aes(xintercept = med_value), linetype = "dashed") +
    geom_text(
      data = med_tbl,
      aes(x = med_value, y = Inf, label = paste0("med=", signif(med_value, 3))),
      vjust = 1.1, size = 3
    ) +
    facet_wrap(~ Sample, scales = "free_y") +
    labs(
      title = paste0(label, " outliers (n = ", n_total, ")"),
      x = value_col, y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text = element_text(face = "bold"),
      plot.title = element_text(hjust = 0.5)
    )

  ggsave(
    filename = file.path(outdir, paste0(label, "_hist_by_sample.png")),
    plot = p, width = 16, height = 12, dpi = 300, bg = "white"
  )

  invisible(p)
}

# POOLED: outliers only, overall median from FULL spe
plot_flag_hist_pooled <- function(flag_col, value_col, label, bins = 50,
                                  width = 10, height = 7, dpi = 300) {
  df_out <- as.data.frame(colData(spe)) %>%
    filter(!is.na(.data[[flag_col]]), .data[[flag_col]]) %>%
    pull(!!rlang::sym(value_col)) %>%
    data.frame(value = .)

  if (nrow(df_out) == 0) return(invisible(NULL))
  n_total     <- nrow(df_out)
  med_overall <- median(colData(spe)[[value_col]], na.rm = TRUE)  # FULL spe

  p <- ggplot(df_out, aes(x = value)) +
    geom_histogram(bins = bins) +
    geom_vline(xintercept = med_overall, linetype = "dashed") +
    annotate("text", x = med_overall, y = Inf,
             label = paste0("median = ", signif(med_overall, 3)),
             vjust = 1.2, size = 3) +
    labs(
      title = paste0(label, " outliers (n = ", n_total, ")"),
      x = value_col, y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      plot.title       = element_text(hjust = 0.5, face = "bold")
    )

  ggsave(
    filename = file.path(outdir, paste0(label, "_hist.png")),
    plot = p, width = width, height = height, dpi = dpi, bg = "white"
  )

  invisible(p)
}

# Make plots for each flag
# sum_gex flags
# plot_flag_hist_by_sample("sum_gex_3MAD", "sum_gex", "sum_gex_3MAD")
# plot_flag_hist_pooled   ("sum_gex_3MAD", "sum_gex", "sum_gex_3MAD")
plot_flag_hist_by_sample("sum_gex_4MAD", "sum_gex", "sum_gex_4MAD")
plot_flag_hist_pooled   ("sum_gex_4MAD", "sum_gex", "sum_gex_4MAD")
# plot_flag_hist_by_sample("sum_gex_5MAD", "sum_gex", "sum_gex_5MAD")
# plot_flag_hist_pooled   ("sum_gex_5MAD", "sum_gex", "sum_gex_5MAD")

# detected_gex flags
# plot_flag_hist_by_sample("detected_gex_3MAD", "detected_gex", "detected_gex_3MAD")
# plot_flag_hist_pooled   ("detected_gex_3MAD", "detected_gex", "detected_gex_3MAD")
plot_flag_hist_by_sample("detected_gex_4MAD", "detected_gex", "detected_gex_4MAD")
plot_flag_hist_pooled   ("detected_gex_4MAD", "detected_gex", "detected_gex_4MAD")
# plot_flag_hist_by_sample("detected_gex_5MAD", "detected_gex", "detected_gex_5MAD")
# plot_flag_hist_pooled   ("detected_gex_5MAD", "detected_gex", "detected_gex_5MAD")

# cell area flags
# plot_flag_hist_by_sample("cell_area_6MAD", "cell_area", "cell_area_6MAD")
# plot_flag_hist_pooled   ("cell_area_6MAD", "cell_area", "cell_area_6MAD")
# plot_flag_hist_by_sample("cell_area_5MAD", "cell_area", "cell_area_5MAD")
# plot_flag_hist_pooled   ("cell_area_5MAD", "cell_area", "cell_area_5MAD")
plot_flag_hist_by_sample("cell_area_4MAD", "cell_area", "cell_area_4MAD")
plot_flag_hist_pooled   ("cell_area_4MAD", "cell_area", "cell_area_4MAD")


# Plot metrics with per-sample medians (from full spe)
plot_metric_by_sample <- function(metric_col, label, bins = 50, logx = FALSE,
                                  width = 16, height = 12, dpi = 300) {
  df <- as.data.frame(colData(spe)) %>%
    select(Sample, !!rlang::sym(metric_col)) %>%
    rename(value = !!rlang::sym(metric_col)) %>%
    filter(!is.na(value))

  if (nrow(df) == 0) return(invisible(NULL))

  # Per-sample medians from full spe
  med_tbl <- df %>%
    group_by(Sample) %>%
    summarise(med_value = median(value, na.rm = TRUE), .groups = "drop")

  p <- ggplot(df, aes(x = value)) +
    geom_histogram(bins = bins, alpha = 0.9) +
    geom_vline(data = med_tbl, aes(xintercept = med_value), linetype = "dashed") +
    geom_text(data = med_tbl,
              aes(x = med_value, y = Inf,
                  label = paste0("med=", signif(med_value, 3))),
              vjust = 1.1, size = 3) +
    facet_wrap(~ Sample, scales = "free_y") +
    labs(
      title = paste0(label, " distribution by Sample"),
      x = label, y = "Count"
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      strip.text       = element_text(face = "bold"),
      plot.title       = element_text(hjust = 0.5)
    )

  if (logx) p <- p + scale_x_log10()

  ggsave(file.path(outdir, paste0(label, "_hist_by_sample.png")),
         plot = p, width = width, height = height, dpi = dpi, bg = "white")
  invisible(p)
}

plot_metric_by_sample("sum_gex",      "sum_gex",      bins = 150, logx = FALSE)
plot_metric_by_sample("detected_gex", "detected_gex", bins = 150, logx = FALSE)
plot_metric_by_sample("cell_area",    "cell_area",    bins = 150, logx = FALSE)

#Visualize outliers using violin plots
# Save each QC plot with an identical data-panel size. Legends and axis text can
# change the total PNG size, but not the width or height of the plotting panel.
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

p <- plotColData(spe, x = "Sample", y = "subsets_any_neg_percent", color_by = "exclude_any_neg") +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 5.5, 25)
  ) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  labs(
    x = "Sample",
    y = "Negative control and unassigned transcripts (%)"
  )
save_qc_plot(
  p,
  here("plots","03_qc","exclude_any_neg_outliers_violin.png")
)


p <- plotColData(spe, x = "Sample", y = "detected_gex", color_by = "detected_gex_4MAD") +
  scale_y_log10() + 
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 5.5, 25)
  ) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  labs(
    x = "Sample",
    y = "Detected gene expression genes per cell"
  )
save_qc_plot(
  p,
  here("plots","03_qc","detected_gex_4MAD_outliers_violin.png")
)

p <- plotColData(spe, x = "Sample", y = "sum_gex", color_by = "sum_gex_4MAD") +
  scale_y_log10() +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.margin = margin(5.5, 5.5, 5.5, 25)
  ) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  labs(
    x = "Sample",
    y = "Gene expression transcripts per cell"
  )
save_qc_plot(
  p,
  here("plots","03_qc","sum_gex_4MAD_outliers_violin.png")
)

p <- plotColData(spe,x = "Sample", y = "cell_area", color_by = "cell_area_4MAD") +
  scale_x_discrete(expand = expansion(add = 0.8)) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.margin = margin(5.5, 5.5, 5.5, 25)
  )  +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) +
  labs(
    x = "Sample",
    y = expression("Cell area ("*mu*"m"^2*")")
  )
save_qc_plot(
  p,
  here("plots","03_qc","cell_area_4MAD_outliers_violin.png")
)

p <- plotColData(spe,x = "Sample", y = "nucleus_area") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(filename = here("plots","03_qc","nucleus_area_by_Sample_violin.png"),plot = p)

# Plot spatial distribution of outliers per sample using SpotSweeper
plotdir <- here("plots","03_qc")
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)

p1 <- SpotSweeper::plotQCpdf(
      spe,
      sample_id  = "Sample",
      metric     = "sum_gex", 
      outliers   = "sum_gex_4MAD",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white"), 
      fname      = paste0(plotdir,"/sum_gex_4MAD_outliers_spatial_by_Sample.pdf")
      )

p2 <- SpotSweeper::plotQCpdf(
      spe, 
      sample_id  = "Sample",
      metric     = "detected_gex", 
      outliers   = "detected_gex_4MAD",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white"),
      fname      = paste0(plotdir,"/detected_gex_4MAD_outliers_spatial_by_Sample.pdf")
      )

p3 <- SpotSweeper::plotQCpdf(
      spe, 
      sample_id  = "Sample",
      metric     = "cell_area",
      outliers   = "cell_area_4MAD",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black"),
      fname      = paste0(plotdir,"/cell_area_4MAD_outliers_spatial_by_Sample.pdf")
      )

# Outlier summary plots
df <- read_csv(here("processed-data","Xenium_QC","global_outliers_summary.csv")) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
df_long <- df %>%
  select(
    Sample,
    exclude_any_neg,
    detected_gex_4MAD,
    sum_gex_4MAD,
    cell_area_4MAD,
    PercentRemoved
  ) %>%
  pivot_longer(
    cols      = -Sample,
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "exclude_any_neg"    ~ "Negative controls ≥25% (cells)",
      Metric == "sum_gex_4MAD"       ~ "Low GEX transcript count (lower 4 MAD; cells)",
      Metric == "detected_gex_4MAD"  ~ "Low detected GEX genes (lower 4 MAD; cells)",
      Metric == "cell_area_4MAD"     ~ "Cell area outliers (two sided 4 MAD; cells)",
      Metric == "PercentRemoved"     ~ "Cells removed (%)",
      TRUE                            ~ Metric
    ),
    Metric = factor(
      Metric,
      levels = c(
        "Negative controls ≥25% (cells)",
        "Low GEX transcript count (lower 4 MAD; cells)",
        "Low detected GEX genes (lower 4 MAD; cells)",
        "Cell area outliers (two sided 4 MAD; cells)",
        "Cells removed (%)"
      )
    )
  )
p <- ggplot(df_long, aes(x = Sample, y = Value)) +
  geom_col() +
  facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5),
    strip.text  = element_text(face = "bold"),
    plot.margin = margin(5, 5, 5, 20)
  ) +
  labs(
    x = NULL,
    y = NULL,
    title = "Per-sample QC flags and % removed"
  )
ggsave(
  here("plots","03_qc","outlier_summary_hist.png"),
  plot   = p,
  width  = 6,
  height = 10,
  dpi    = 300
)

#How do nuclei/cell area correlate with counts + detected? 
#Nuc area by UMIs
png(here("plots","03_qc","Nucleus_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "sum")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","Nucleus_area_by_detected.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "detected")
dev.off()

#Nuc area by total_counts
png(here("plots","03_qc","Nucleus_area_by_totalcounts.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "total_counts")
dev.off()

#Cell area by UMIs
png(here("plots","03_qc","cell_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "sum")
dev.off()

#Cell area by nGenes
png(here("plots","03_qc","cell_area_by_detected.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "detected")
dev.off()

#Cell area by total_counts
png(here("plots","03_qc","Cell_area_by_totalcounts.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "total_counts")
dev.off()


#Split these plots by sample using facet_wrap
coldata_df <- as.data.frame(colData(spe))

#Nuc area by sum
png(here("plots","03_qc","Nucleus_area_by_sum_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = nucleus_area, y = sum)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","Nucleus_area_by_detected_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = nucleus_area, y = detected)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Cell area by sum
png(here("plots","03_qc","cell_area_by_sum_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = cell_area, y = sum)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","cell_area_by_detected_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = cell_area, y = detected)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()


######################### PLOTTING #########################


# # Determine adaptive cutoffs using Gaussian mixture modeling (GMM)
# library(mclust)

# decide_cutoffs <- function(x, nmads_low = 3, nmads_high = 4) {
#   lx <- log1p(x)
#   # Fit 1 vs 2 components
#   m1 <- Mclust(lx, G = 2, verbose = FALSE)
#   m2 <- Mclust(lx, G = 3, verbose = FALSE)
#   bic1 <- m1$bic; bic2 <- m2$bic
#   bimodal <- (bic2 - bic1) >= 10  # prefer 2 components if ΔBIC ≥ 10

#   if (bimodal) {
#     # use the higher-mean component as "cell" cluster
#     k <- which.max(m2$parameters$mean)
#     post <- m2$z[, k]
#     lx_c <- lx[post >= 0.5]  # confident cells
#   } else {
#     lx_c <- lx
#   }

#   med <- median(lx_c, na.rm = TRUE)
#   mad <- median(abs(lx_c - med), na.rm = TRUE)  # unscaled MAD
#   low  <- exp(med - nmads_low  * mad) - 1
#   high <- exp(med + nmads_high * mad) - 1

#   list(bimodal = bimodal, bic1 = bic1, bic2 = bic2,
#        low = low, high = high,
#        center = exp(med) - 1, mad_log = mad)
# }

# # Run per sample for sum_gex and detected_gex
# qc_bounds <- lapply(split(as.data.frame(colData(spe)), spe$Sample), function(df) {
#   s1 <- decide_cutoffs(df$sum_gex)
#   s2 <- decide_cutoffs(df$detected_gex)
#   tibble(
#     Sample        = unique(df$Sample),
#     sum_bimodal   = s1$bimodal, sum_low  = s1$low, sum_high  = s1$high,
#     detected_bimodal = s2$bimodal, detected_low = s2$low, detected_high = s2$high
#   )
# }) %>% bind_rows()

# qc_bounds

######################### Save cleaned SPE #########################
#Save spe before removing low quality cells
message(paste0("Saving spe object with QC - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_withQC.Rds"))

#Remove all of the low quality cells
table(spe$global_outliers)
spe <- spe[,!spe$global_outliers]
dim(spe)

#Save cleaned SPE
message(paste0("Saving cleaned SPE object - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_clean.Rds"))

# Check number of cells with low counts and detected genes after filtering
table(spe$sum_gex < 10)
table(spe$detected_gex < 4)
table(spe$subsets_any_neg_percent >= 25)

spe$sum_gex_tag <- factor(
  ifelse(spe$sum_gex %in% 1:9, as.character(spe$sum_gex), ">=10"),
  levels = c(as.character(1:9), ">=10")
)
table(spe$Sample, spe$sum_gex_tag)

spe$detected_gex_tag <- factor(
  ifelse(spe$detected_gex %in% 1:3, as.character(spe$detected_gex), ">=4"),
  levels = c(as.character(1:3), ">=4")
)
table(spe$Sample, spe$detected_gex_tag)

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
