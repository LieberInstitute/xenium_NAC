#!/usr/bin/env Rscript
#Goal: Calculate QC metrics and explore. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
#Modified from:
 # https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/f3b87d7a0356e716468b304196f5cd362fd91515/code/analysis/02_xenium_qc/01_threshold_outliers.R
 # https://github.com/LieberInstitute/spatialAmygdala/blob/32108cdc145bb83aa73822c862dc28299a29b47c/code/Xenium/03_quality_control/01_perCellQC.R
 # https://pachterlab.github.io/voyager/articles/vig5_xenium.html#quality-control
 # http://127.0.0.1:29231/library/SpotSweeper/doc/getting_started.html#identifying-local-outliers-using-spotsweeper 

args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("No sample name provided. Please provide a sample name as an argument.")
} else {
  samp <- args[1]
}

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

  
#Read in the RDS file from 01_build_spe. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))
spe_sub <- spe[, spe$Sample == samp]
# spe

#Make the sample column a factor
spe_sub$Sample <- factor(x = spe_sub$Sample, levels = unique(spe_sub$Sample))
# levels(spe_sub$Sample)

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe_sub), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe_sub), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe_sub), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe_sub)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel


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
spe_sub <- scuttle::addPerCellQCMetrics(spe_sub, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX))

spe_sub$exclude_low_lib  <- colSums(counts(spe_sub)) == 0
spe_sub$exclude_all_neg  <- spe_sub$subsets_any_neg_percent == 100

spe_sub$global_outliers <- as.logical(spe_sub$exclude_low_lib) |
                           as.logical(spe_sub$exclude_all_neg)

# save outlier tables
outdir <- here("processed-data","Xenium_QC","SpotSweeper_outliers")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cd <- as.data.frame(colData(spe_sub))
df <- cd[ cd$global_outliers , , drop = FALSE ]
write.csv(df, file = file.path(outdir, paste0(samp,"_global_outliers.csv")), row.names = TRUE)


# Remove cells with extremely low library size
n0 <- ncol(spe_sub)
spe_sub <- spe_sub[, colSums(counts(spe_sub)) > 0]
n1 <- ncol(spe_sub)
message(sprintf("Dropped %d cells with extremely low library size (from %d down to %d)", 
                n0 - n1, n0, n1))


#There are some cells within this object in which 100% of the reads are negative controls - Remove them
spe_sub <- spe_sub[,spe_sub$subsets_any_neg_percent < 100]
n2 <- ncol(spe_sub)
message(sprintf("Dropped %d cells with 100%% negative control reads (from %d down to %d)", 
                n1 - n2, n1, n2))

colnames(colData(spe_sub))


# subset to this sample
message("Processing sample: ", samp, " (n cells = ", ncol(spe_sub), ")")

# compute local outliers on library size within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "sum",
  direction = "lower",
  log       = TRUE,
  cutoff    = 3
)
# compute local outliers on any_neg within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "subsets_any_neg_percent",
  direction = "higher",
  log       = FALSE,
  cutoff    = 3
)
# compute local outliers on cell area within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "cell_area",
  direction = "both",
  log       = FALSE,
  cutoff    = 3
)
# compute local outliers on detected genes within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "detected",
  direction = "lower",
  log       = TRUE,
  cutoff    = 3
)

spe_sub$local_outliers <- as.logical(spe_sub$sum_outliers) |
                          as.logical(spe_sub$subsets_any_neg_percent_outliers) |
                          as.logical(spe_sub$cell_area_outliers) |
                          as.logical(spe_sub$detected_outliers)

cd <- as.data.frame(colData(spe_sub))
df <- cd[ cd$local_outliers , , drop = FALSE ]
write.csv(df, file = file.path(outdir, paste0(samp,"_local_outliers.csv")), row.names = TRUE)

# plot & save
plotdir <- here("plots","03_qc", "SpotSweeper_outliers")
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)

p1 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "sum_log", 
      outliers   = "sum_outliers",
      point_size = 0.5,
      stroke     = 0.5) +
      ggtitle(paste0("Library size outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_library_size_outliers.png")),
       p1, width=5, height=5, dpi=300)

p2 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "subsets_any_neg_percent",
      outliers   = "subsets_any_neg_percent_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white")) +
      ggtitle(paste0("Any neg outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_any_neg_percent_outliers.png")),
       p2, width=5, height=5, dpi=300)

p3 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "cell_area",
      outliers   = "cell_area_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black")) +
      ggtitle(paste0("Cell area outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_cell_area_outliers.png")),
       p3, width=5, height=5, dpi=300)

p4 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "detected",
      outliers   = "detected_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black")) +
      ggtitle(paste0("Detected genes outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_detected_genes_outliers.png")),
       p4, width=5, height=5, dpi=300)

p5 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "sum_log", 
      outliers   = "local_outliers",
      point_size = 0.5,
      stroke     = 0.5) +
      ggtitle(paste0("Local outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_outliers.png")),
       p5, width=5, height=5, dpi=300)

# reproducibility
message("Session info for ", samp,":")
print(session_info())


