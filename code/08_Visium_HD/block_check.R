#Goal: Run Banksy spatial clustering with parameters equating to spatial clustering for domain finding
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(DelayedArray)
library(sessioninfo)
library(HDF5Array)
library(ggplot2)
library(escheR)
library(Banksy)
library(scran)
library(here)

here::i_am(path = "code/08_Visium_HD/block_check.R")

#Read in the filtered sfe object
message(Sys.time(), " | Loading the object")
sfe_filtered_dir <- here(
  "processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2"
)
sfe <- loadHDF5SummarizedExperiment(sfe_filtered_dir)

####
#Calcualte log-normalized counts
message(Sys.time(), " | Calcualte log-normalized counts")
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)

####HVGs 
message(Sys.time(), " | Selecting top HVGs without block")
n_hvgs <- 2000
gene_var <- modelGeneVar(sfe)
top_hvgs <- getTopHVGs(gene_var, n = n_hvgs)

message(Sys.time(), " | Selecting top HVGs with block")
gene_var2 <- modelGeneVar(sfe,block = sfe$sample_id)
top_hvgs2 <- getTopHVGs(gene_var2, n = n_hvgs)

message(Sys.time(), " | HVG intersection")
length(intersect(top_hvgs,top_hvgs2))

message(Sys.time(), " | HVG difference")
length(setdiff(top_hvgs,top_hvgs2))
setdiff(top_hvgs,top_hvgs2)

message(Sys.time(), " | done")
session_info()
