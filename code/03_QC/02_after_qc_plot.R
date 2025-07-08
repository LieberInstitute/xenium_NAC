#!/usr/bin/env Rscript
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
# library(scattermore)
library(tidyverse)
# library(scuttle)
library(scater)
library(escheR)
# library(scran)
library(dplyr)
library(here)
library(SpotSweeper)

spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))
spe$Sample <- factor(x = spe$Sample, levels = unique(spe$Sample))

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX))

outlier_dir <- here("processed-data","Xenium_QC","Sample_outliers")
global_outliers_files  <- list.files(outlier_dir, pattern = "_global_outliers\\.csv$", full.names = TRUE)
local_outliers_files  <- list.files(outlier_dir, pattern = "_local_outliers\\.csv$", full.names = TRUE)

global_outliers <- unique(unlist(lapply(global_outliers_files, function(f) {
  global_df <- read.csv(f, row.names = 1)
  rownames(global_df)
})))
local_outliers <- unique(unlist(lapply(local_outliers_files, function(f) {
  local_df <- read.csv(f, row.names = 1)
  rownames(local_df)
})))
spe$global_outliers <- colnames(spe) %in% global_outliers
spe$local_outliers <- colnames(spe) %in% local_outliers
spe$outliers <- as.logical(spe$global_outliers) |
                as.logical(spe$local_outliers)

spe$exclude_low_lib                       <- FALSE
spe$exclude_any_neg                       <- FALSE
spe$cell_area_outliers                    <- FALSE
spe$subsets_any_neg_percent_outliers      <- FALSE

for (f in global_outliers_files) {
  df    <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  cells <- rownames(df)
  idx <- match(cells, colnames(spe))
  spe$exclude_low_lib   [idx] <- df$exclude_low_lib
  spe$exclude_any_neg   [idx] <- df$exclude_any_neg
  spe$cell_area_outliers[idx] <- df$cell_area_outliers
}

for (f in local_outliers_files) {
  df    <- read.csv(f, row.names = 1, stringsAsFactors = FALSE)
  cells <- rownames(df)
  idx <- match(cells, colnames(spe))
  spe$subsets_any_neg_percent_outliers[idx] <- df$subsets_any_neg_percent_outliers
}

#Save discard metrics as a csv
# table(spe$Sample,spe$outliers) %>% 
#   as.data.frame.matrix() %>%
#   mutate(Sample = rownames(.),
#          CellsRetained = `FALSE`,
#          Outliers = `TRUE`,
#          PercentRemoved = Outliers/(Outliers+CellsRetained)*100) %>%
#   select(Sample,CellsRetained,Outliers,PercentRemoved) %>%
#   `rownames<-`(NULL) %>%
#   write.csv(file = here("processed-data","Xenium_QC","All_Outliers.csv"),
#             quote = FALSE)

as.data.frame(colData(spe)) %>%
  group_by(Sample) %>%
  summarize(
    CellsRetained                        = sum(!outliers),
    Outliers                             = sum(outliers),
    PercentRemoved                       = Outliers/(CellsRetained + Outliers)*100,
    exclude_low_lib                      = sum(exclude_low_lib),
    exclude_any_neg                      = sum(exclude_any_neg),
    cell_area_outliers                   = sum(cell_area_outliers),
    subsets_any_neg_percent_outliers     = sum(subsets_any_neg_percent_outliers)
  ) %>%
  ungroup() %>%
  write.csv(
    file      = here("processed-data","Xenium_QC","All_Outliers.csv"),
    row.names = FALSE,
    quote     = FALSE
  )

df <- read_csv(here("processed-data","Xenium_QC","All_Outliers.csv")) %>%
  mutate(Sample = factor(Sample, levels = unique(Sample)))
df_long <- df %>%
  select(
    Sample,
    exclude_low_lib,
    exclude_any_neg,
    cell_area_outliers,
    subsets_any_neg_percent_outliers,
    PercentRemoved
  ) %>%
  pivot_longer(
    cols      = -Sample,
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  mutate(
    Metric = case_when(
      Metric == "exclude_low_lib"                       ~ "Low-lib (=0) outliers (count)",
      Metric == "exclude_any_neg"                       ~ "High-neg% (≥25%) outliers (count)",
      Metric == "cell_area_outliers"                    ~ "Cell-area outliers (±6 MAD) (count)",
      Metric == "subsets_any_neg_percent_outliers"      ~ "Local neg% outliers (SpotSweeper) (count)",
      Metric == "PercentRemoved"                        ~ "Percent of cells removed (%)",
      TRUE                                              ~ Metric
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
  here("plots","03_qc","outlier_counts_and_percent_removed_by_sample.png"),
  plot   = p,
  width  = 6,
  height = 10,
  dpi    = 300
)


# Plotting the outliers
p1 <- plotColData(spe, x = "Sample", y = "sum", color_by = "outliers") +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p1, filename = here("plots","03_qc","outliers_total_counts_violin.png"))

# p2 <- plotColData(spe, x = "Sample", y = "detected", color_by = "outliers") +
#   scale_y_log10() +
#   theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
#   stat_summary(fun = median, 
#                fun.min = median, 
#                fun.max = median,
#                geom = "crossbar", 
#                width = 0.3)
# ggsave(p2, filename = here("plots","03_qc","outliers_detected_genes_violin.png"))

p3 <- plotColData(spe, x = "Sample", y = "subsets_any_neg_percent", color_by = "outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p3, filename = here("plots","03_qc","outliers_any_neg_percent_violin.png"))

p4 <- plotColData(spe, x = "Sample", y = "cell_area", color_by = "outliers") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p4, filename = here("plots","03_qc","outliers_cell_area_violin.png"))    


## Check with initial clustering
spe_outlier <- spe[,spe$local_outliers]
# Banksy clustering
banksy_cluster <- read.csv(here("processed-data","05_Clustering","Banksy_clusters.csv"), row.names = 1)
colnames(banksy_cluster) <- c("cluster", "cell_id")
row.names(banksy_cluster) <- banksy_cluster$cell_id

common_keys <- intersect(colnames(spe_outlier), banksy_cluster$cell_id)
spe_common <- spe_outlier[, common_keys]
common_idx <- match(common_keys, banksy_cluster$cell_id)
spe_common$banksy <- banksy_cluster$cluster[common_idx]

outlier_cluster_banksy <- data.frame(cluster = spe_common$banksy)
p6 <- ggplot(outlier_cluster_banksy, aes(x = factor(cluster))) +
  geom_bar(fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Banksy cluster",
    x     = "Cluster ID",
    y     = "Number of outlier cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(
  filename = here("plots","03_qc","banksy_outliers_histogram.png"),
  plot     = p6,
  width    = 6,
  height   = 4,
  dpi      = 300
)

# STAligner clustering
spe_outlier$Sample <- as.character(spe_outlier$Sample)
spe_outlier$Sample[spe_outlier$Sample == "Br6660_Nac10_4080"] <- "Br6660_NAc10_4080"
spe_outlier$Sample[spe_outlier$Sample == "Br6660_Nac11_5580"] <- "Br6660_NAc11_5580"
spe_outlier$Sample <- factor(spe_outlier$Sample, levels = unique(spe_outlier$Sample))

staligner_r0_mclust20 <- read.csv(here("processed-data","05_Clustering","STAligner_initial_r0_mclust20.csv"), row.names = 1)
staligner_r15_mclust20 <- read.csv(here("processed-data","05_Clustering","STAligner_initial_r15_mclust20.csv"), row.names = 1)

spe_keys <- paste(spe_outlier$Sample, spe_outlier$cell_id, sep = "_")
staligner_keys  <- paste(staligner_r0_mclust20$slice_name, staligner_r0_mclust20$cell_id, sep = "_")
common_keys <- intersect(spe_keys, staligner_keys)
common_idx_spe <- which(spe_keys %in% common_keys)
common_idx_st  <- match(common_keys, staligner_keys)
spe_common <- spe_outlier[, common_idx_spe]
spe_common$mclust20_r0 <- staligner_r0_mclust20$mclust_20[common_idx_st]
spe_common$mclust20_r15 <- staligner_r15_mclust20$mclust_20[common_idx_st]

outlier_cluster_r0 <- data.frame(cluster = spe_common$mclust20_r0)
p7 <- ggplot(outlier_cluster_r0, aes(x = factor(cluster))) +
  geom_bar(fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "STAligner cluster (r = 0)",
    x     = "Cluster ID",
    y     = "Number of outlier cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(
  filename = here("plots","03_qc","staligner_r0_mclust20_outliers_histogram.png"),
  plot     = p7,
  width    = 6,
  height   = 4,
  dpi      = 300
)

outlier_cluster_r15 <- data.frame(cluster = spe_common$mclust20_r15)
p8 <- ggplot(outlier_cluster_r15, aes(x = factor(cluster))) +
  geom_bar(fill = "steelblue", color = "white") +
  theme_minimal(base_size = 14) +
  labs(
    title = "STAligner cluster (r = 15)",
    x     = "Cluster ID",
    y     = "Number of outlier cells"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
ggsave(
  filename = here("plots","03_qc","staligner_r15_mclust20_outliers_histogram.png"),
  plot     = p8,
  width    = 6,
  height   = 4,
  dpi      = 300
)

#Remove all of the low quality cells
spe <- spe[,!spe$outliers]

#Save cleaned SPE
message(paste0("Saving cleaned SPE object - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_clean.Rds"))


