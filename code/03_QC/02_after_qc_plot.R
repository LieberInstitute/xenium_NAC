#!/usr/bin/env Rscript
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
library(SpatialExperiment)
# library(scattermore)
# library(tidyverse)
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

outlier_dir <- here("processed-data","Xenium_QC","SpotSweeper_outliers")
global_outliers_files  <- list.files(outlier_dir, pattern = "_global_outliers\\.csv$", full.names = TRUE)
local_outliers_files  <- list.files(outlier_dir, pattern = "_local_outliers\\.csv$", full.names = TRUE)

global_outliers <- unique(unlist(lapply(global_outliers_files, function(f) {
  df <- read.csv(f, row.names = 1)
  rownames(df)
})))
local_outliers <- unique(unlist(lapply(local_outliers_files, function(f) {
  df <- read.csv(f, row.names = 1)
  rownames(df)
})))
spe$global_outliers <- colnames(spe) %in% global_outliers
spe$local_outliers <- colnames(spe) %in% local_outliers
spe$outliers <- as.logical(spe$global_outliers) |
                as.logical(spe$local_outliers)

#Save discard metrics as a csv
table(spe$Sample,spe$outliers) %>% 
  as.data.frame.matrix() %>%
  mutate(Sample = rownames(.),
         CellsRetained = `FALSE`,
         Outliers = `TRUE`,
         PercentRemoved = Outliers/(Outliers+CellsRetained)*100) %>%
  select(Sample,CellsRetained,Outliers,PercentRemoved) %>%
  `rownames<-`(NULL) %>%
  write.csv(file = here("processed-data","Xenium_QC","All_Outliers.csv"),
            quote = FALSE)

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

p2 <- plotColData(spe, x = "Sample", y = "detected", color_by = "outliers") +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p2, filename = here("plots","03_qc","outliers_detected_genes_violin.png"))

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

# plotdir <- here("plots","03_qc", "SpotSweeper_outliers")
# p5 <- SpotSweeper::plotQCmetrics(
#       spe, 
#       metric     = "sum", 
#       outliers   = "outliers",
#       sample_id  = "Sample",
#       sample     = "Br6660_NAc4_2080",
#       point_size = 0.5,
#       stroke     = 0.5) +
#       ggtitle(paste0("All outliers: ", "Br6660_NAc4_2080"))
# ggsave(file.path(plotdir, paste0("Br6660_NAc4_2080", "_all_outliers.png")),
#        p5, width=5, height=5, dpi=300)


#Remove all of the low quality cells
spe <- spe[,!spe$outliers]

#Save cleaned SPE
message(paste0("Saving cleaned SPE object - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_clean.Rds"))