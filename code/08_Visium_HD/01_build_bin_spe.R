#Goal: Build initial Visium-HD SPE object. 
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
# Code modified from
#.        https://github.com/LieberInstitute/human_VTA_spatial/blob/devel/code/02_Build_Spe/01_Build_Spe.R
# and 
#          https://github.com/LieberInstitute/lc_visium_hd/blob/a13f88b6b66bc0fd643a8d690694c454905b3eae/code/03_build_spe/01_build_bin_spe.sh
library(spatialLIBD)
library(sessioninfo)
library(tidyverse)
library(HDF5Array)
library(here)
library(scran)


#Create a dataframe/csv containing information about the samoples
sample_info <- read.csv(here("processed-data","visiumHD_sample_info_NAc.csv"))
sample_info


## Output directories
spe_raw_path  <- here("processed-data","HD_Full_Analysis","SPEs","spe_raw.Rds")
spe_norm_path <- here("processed-data","HD_Full_Analysis","SPEs","spe_norm.Rds")


## Read sample info
sample_ids <- sample_info$sample_id
sr_out_dirs <- here(sample_info$spaceranger_dir,"outs","binned_outputs","square_008um")
reference_gtf <- "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"

################################################################################
#   Compute 8um bin-level SpatialExperiment object
################################################################################

# Hack around 'read10xVisium's requirement for the 'outs' directory to be the immediate parent to 'spatial' directory and other outputs (create a symlink
# named 'outs'). Note we already handled the other required workaround: create a tissue_positions.csv file (not parquet format) (see sh script)
temp_sr_dirs <- file.path(tempdir(), sample_ids, 'outs')
for (this_dir in temp_sr_dirs) {
  dir.create(dirname(this_dir))
}

file.symlink(sr_out_dirs, temp_sr_dirs) |>
  all() |>
  stopifnot()

# Note providing the reference GTF here is mandatory (without doing so, read10xVisiumWrapper searches for a web summary that doesn't exist to try to infer the GTF)
message(Sys.time(), ' | Building SpatialExperiment...')
spe <- read10xVisiumWrapper(
  samples = temp_sr_dirs,
  sample_id = sample_ids,
  type = "sparse",
  data = "raw",
  images = "lowres",
  load = TRUE,
  reference_gtf = reference_gtf
)

message(Sys.time(), " | Saving raw SPE")
saveRDS(spe, spe_raw_path)

################################################################################
#   Filtering SPE
################################################################################
# Filter raw SPE: take only bins in tissue, drop bins with 0 counts for all genes, and drop genes with 0 counts in every bin
message(Sys.time(), ' | Filtering bins, genes, and to tissue...')
spe <- spe[
  rowSums(assays(spe)$counts) > 0,
  (colSums(assays(spe)$counts) > 0) & spe$in_tissue
]

################################################################################
#   Library-Size Normalization
################################################################################
# Use library-size normalization (normalization by deconvolution is not computationally feasible with data this large)

message(Sys.time(), ' | Performing log normalization...')
spe <- computeLibraryFactors(spe)
spe <- logNormCounts(spe)

################################################################################
#   Save object after normalization
################################################################################

message(Sys.time(), " | Saving normalized SPE")
saveRDS(spe, spe_norm_path)

## Size in Gb
lobstr::obj_size(spe)

#-------------------------- Reproducibility information -----------------------#

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
session_info()
