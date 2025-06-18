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
# spe

#Make the sample column a factor
spe$Sample <- factor(x = spe$Sample, levels = unique(spe$Sample))
# levels(spe$Sample)

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

# Remove empty cells
empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]
spe <- spe[, colSums(counts(spe)) > 0]
# dim(spe)

#There are some cells within this object in which 100% of the reads are negative controls - Remove them
spe <- spe[,spe$subsets_any_neg_percent < 100]
# dim(spe)

colnames(colData(spe))


# subset to this sample
spe_sub <- spe[, spe$Sample == samp]
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

spe_sub$local_outliers <- as.logical(spe_sub$sum_outliers) |
    as.logical(spe_sub$subsets_any_neg_percent_outliers)

# save outlier tables
outdir <- here("processed-data","03_qc","outliers")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cd <- as.data.frame(colData(spe_sub))
for(col in c("sum_outliers","subsets_any_neg_percent_outliers","local_outliers")) {
  df <- cd[ which(cd[[col]]), , drop=FALSE ]
  write.csv(
    df,
    file = file.path(outdir, paste0(samp,"_", col, ".csv")),
    row.names = TRUE
  )
}

# plot & save
plotdir <- here("plots","03_qc")
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)

p1 <- plotQCmetrics(spe_sub, "sum_log", outliers = "sum_outliers") +
      ggtitle(paste0("Library‐Size Outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_library_size_outliers.png")),
       p1, width=5, height=5, dpi=300)

p2 <- plotQCmetrics(spe_sub, "subsets_any_neg_percent",
                    outliers="subsets_any_neg_percent_outliers") +
      ggtitle(paste0("Any‐Neg Outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_any_neg_percent_outliers.png")),
       p2, width=5, height=5, dpi=300)

p3 <- plotQCmetrics(spe_sub, "sum_log", outliers="local_outliers") +
      ggtitle(paste0("Local Outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_outliers.png")),
       p3, width=5, height=5, dpi=300)

# 9) reproducibility
message("Session info for ", samp,":")
print(session_info())




