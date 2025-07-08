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
library(sessioninfo)

  
#Read in the RDS file from 01_build_spe. 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))
# spe_highneg <- spe[, !is.na(spe$subsets_any_neg_percent) & spe$subsets_any_neg_percent >= 25 ] # 365
# table(spe_highneg$sum)
# #  1   2   3   4   5   6   7   8 
# # 41  73  85 147   4   7   5   3 
# table(spe_highneg$subsets_any_neg_percent)

#               25 28.5714285714286 33.3333333333333               40 
#              147                5               92                3 
#               50               60              100 
#               76                1               41 

# spe_lowcount <- spe[, spe$sum <= 2] # 47855
# table(spe_lowcount$sum)
# #    0     1     2 
# # 20395 13375 14085 

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
# There are NAs in the 'subsets_any_neg_percent' because of 0 counts.
spe_sub$exclude_any_neg  <- spe_sub$subsets_any_neg_percent >= 25
spe_sub$exclude_any_neg[is.na(spe_sub$exclude_any_neg)] <- FALSE

spe_sub$global_outliers <- as.logical(spe_sub$exclude_low_lib) |
                           as.logical(spe_sub$exclude_any_neg)

# save outlier tables
outdir <- here("processed-data","Xenium_QC","SpotSweeper_outliers")
dir.create(outdir, recursive=TRUE, showWarnings=FALSE)

cd <- as.data.frame(colData(spe_sub))
df <- cd[ cd$global_outliers , , drop = FALSE ]
write.csv(df, file = file.path(outdir, paste0(samp,"_global_outliers.csv")), row.names = TRUE)

# plot & save
plotdir <- here("plots","03_qc", "SpotSweeper_outliers")
dir.create(plotdir, recursive=TRUE, showWarnings=FALSE)

p6 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "sum", 
      outliers   = "exclude_low_lib",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white")) +
      ggtitle(paste0("Global low lib outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_global_library_size_outliers.png")),
       p6, width=5, height=5, dpi=300)

p7 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "subsets_any_neg_percent", 
      outliers   = "exclude_any_neg",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white")) +
      ggtitle(paste0("Global any neg outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_global_any_neg_percent_outliers.png")),
       p7, width=5, height=5, dpi=300)

# Remove cells with extremely low library size
# There are some cells within this object in which a high percentage of the reads are negative controls - Remove them
spe_sub <- spe_sub[, !spe_sub$global_outliers]


# subset to this sample
message("Processing sample: ", samp, " (n cells = ", ncol(spe_sub), ")")

# compute local outliers on library size within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "sum",
  direction = "lower",
  log       = TRUE,
  cutoff    = 3,
  n_neighbors = 1000
)
# compute local outliers on any_neg within this sample
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "subsets_any_neg_percent",
  direction = "higher",
  log       = FALSE,
  cutoff    = 3,
  n_neighbors = 1000
)
# compute local outliers on cell area within this sample
# spe_sub <- SpotSweeper::localOutliers(
#   spe_sub,
#   metric    = "cell_area",
#   direction = "both",
#   log       = FALSE,
#   cutoff    = 3.5,
#   n_neighbors = 1000
# )

spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "cell_area",
  direction = "lower",
  log       = FALSE,
  cutoff    = 3.5,
  n_neighbors = 1000
)
spe_sub$cell_area_lower_outliers <- spe_sub$cell_area_outliers

spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric    = "cell_area",
  direction = "higher",
  log       = FALSE,
  cutoff    = 3.5,
  n_neighbors = 1000
)
spe_sub$cell_area_higher_outliers <- spe_sub$cell_area_outliers

# compute local outliers on detected genes within this sample
# spe_sub <- SpotSweeper::localOutliers(
#   spe_sub,
#   metric    = "detected",
#   direction = "lower",
#   log       = TRUE,
#   cutoff    = 3
# )

spe_sub$local_outliers <- as.logical(spe_sub$sum_outliers) |
                          as.logical(spe_sub$subsets_any_neg_percent_outliers) |
                          as.logical(spe_sub$cell_area_higher_outliers)        |
                          as.logical(spe_sub$cell_area_lower_outliers)        #             |             as.logical(spe_sub$detected_outliers)

cd <- as.data.frame(colData(spe_sub))
df <- cd[ cd$local_outliers , , drop = FALSE ]
write.csv(df, file = file.path(outdir, paste0(samp,"_local_outliers.csv")), row.names = TRUE)

p1 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "sum_log", 
      outliers   = "sum_outliers",
      point_size = 0.5,
      stroke     = 0.5) +
      ggtitle(paste0("Library size outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_library_size_outliers.png")),
       p1, width=5, height=5, dpi=300)

p2 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "subsets_any_neg_percent",
      outliers   = "subsets_any_neg_percent_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white")) +
      ggtitle(paste0("Any neg outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_any_neg_percent_outliers.png")),
       p2, width=5, height=5, dpi=300)

# p3 <- SpotSweeper::plotQCmetrics(
#       spe_sub, 
#       metric     = "cell_area",
#       outliers   = "cell_area_outliers",
#       point_size = 0.5,
#       stroke     = 0.5,
#       colors     = c("white","black")) +
#       ggtitle(paste0("Cell area outliers: ", samp))
# ggsave(file.path(plotdir, paste0(samp, "_local_cell_area_outliers.png")),
#        p3, width=5, height=5, dpi=300)

p3 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "cell_area",
      outliers   = "cell_area_higher_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black")) +
      ggtitle(paste0("Cell area outliers (higher): ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_cell_area_higher_outliers.png")),
       p3, width=5, height=5, dpi=300)

p3 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "cell_area",
      outliers   = "cell_area_lower_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black")) +
      ggtitle(paste0("Cell area outliers (lower): ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_cell_area_lower_outliers.png")),
       p3, width=5, height=5, dpi=300)

# p4 <- SpotSweeper::plotQCmetrics(
#       spe_sub, 
#       metric     = "detected",
#       outliers   = "detected_outliers",
#       point_size = 0.5,
#       stroke     = 0.5,
#       colors     = c("white","black")) +
#       ggtitle(paste0("Detected genes outliers: ", samp))
# ggsave(file.path(plotdir, paste0(samp, "_detected_genes_outliers.png")),
#        p4, width=5, height=5, dpi=300)

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









spe_sub$has_any_neg <- spe_sub$subsets_any_neg_percent != 0
table(spe_sub$has_any_neg, useNA="ifany")



p2 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "subsets_any_neg_percent",
      outliers   = "has_any_neg",
      point_size = 0.1,
      stroke     = 0.1,
      colors     = c("black","white")) +
      ggtitle(paste0("Any neg outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_has_any_neg.png")),
       p2, width=5, height=5, dpi=300)








# 1) subset to cells with any_neg > 0
hasNeg <- spe_sub$subsets_any_neg_percent > 0
spe_nz <- spe_sub[, hasNeg]

# # 2) jitter their percent slightly so MAD ≠ 0
# devs   <- abs(spe_nz$subsets_any_neg_percent - median(spe_nz$subsets_any_neg_percent))
# min_nz <- min(devs[devs>0], na.rm=TRUE)
# eps    <- min_nz/10
# set.seed(42)
# jitter <- sample(seq(0, eps, length.out = ncol(spe_nz)))
# spe_nz$negPct_buf <- spe_nz$subsets_any_neg_percent + jitter

# 3) run SpotSweeper locally
spe_nz <- localOutliers(
  spe_nz,
  metric      = "subsets_any_neg_percent",
  direction   = "higher",
  log         = FALSE,
  n_neighbors = 100,
  cutoff      = 4
)

# 4) pull back the flag into your full SPE
spe_sub$subsets_any_neg_percent_outliers <- FALSE
pos_in_sub <- match(colnames(spe_nz), colnames(spe_sub))
spe_sub$subsets_any_neg_percent_outliers[pos_in_sub] <- spe_nz$subsets_any_neg_percent_outliers



p2 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "subsets_any_neg_percent",
      outliers   = "subsets_any_neg_percent_outliers",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("black","white")) +
      ggtitle(paste0("Any neg outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_local_any_neg_percent_outliers.png")),
       p2, width=5, height=5, dpi=300)




# pull out colData into a data.frame
df <- as.data.frame(colData(spe_sub))
df_out <- df[df$subsets_any_neg_percent_outliers, ]
p <- ggplot(df_out, aes(x = subsets_any_neg_percent)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of % Negative Probes for Outlier Cells",
    x     = "subsets_any_neg_percent",
    y     = "Number of Cells"
  )
ggsave(file.path(plotdir, paste0(samp, "_histogram_local_any_neg_percent_outliers.png")),
       p, width=5, height=5, dpi=300)


table(spe_sub$subsets_any_neg_percent_outliers, useNA="ifany")
spe_remained <- spe_sub[, !spe_sub$subsets_any_neg_percent_outliers]
summary(spe_remained$subsets_any_neg_percent)



####### cell area
df <- as.data.frame(colData(spe_sub))
p <- ggplot(df, aes(x = cell_area)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  theme_minimal() +
  labs(
    title = "Distribution of Cell Area",
    x     = "Cell Area",
    y     = "Number of Cells"
  )
ggsave(file.path(plotdir, paste0(samp, "_histogram_cell_area.png")),
       p, width=5, height=5, dpi=300)



# 1) global thresholds
spe_sub$cell_area_outliers <- isOutlier(
  spe_sub$cell_area,
  log = FALSE,
  # batch = as.factor(spe$Sample),
  nmads = 5.5,
  type = "both"
)

table(spe_sub$cell_area_outliers)

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



library(mgcv)

# 1) extract
df       <- as.data.frame(colData(spe_sub))
df$area  <- spe_sub$cell_area

# 2) fit GAM
fit      <- gam(area ~ s(spatialCoords(spe_sub)[,1], spatialCoords(spe_sub)[,2]), data = df)

# 3) residuals
spe_sub$area_resid <- residuals(fit)

# 4) local outlier on residual
spe_sub <- SpotSweeper::localOutliers(
  spe_sub,
  metric      = "area_resid",
  direction   = "both",
  log         = FALSE,
  n_neighbors = 50,
  cutoff      = 3
)

table(spe_sub$area_resid_outliers)


p3 <- SpotSweeper::plotQCmetrics(
      spe_sub, 
      metric     = "cell_area",
      outliers   = "global_areaMAD3",
      point_size = 0.5,
      stroke     = 0.5,
      colors     = c("white","black")) +
      ggtitle(paste0("Cell area outliers: ", samp))
ggsave(file.path(plotdir, paste0(samp, "_global_areaMAD3_outliers.png")),
       p3, width=5, height=5, dpi=300)