#Goal: Add RCTD weights and predictions to the SPE
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
library(SpatialExperiment)
library(HDF5Array)
library(escheR)
library(here)

## Load the full SPE
message("Loading SPE object - ", Sys.time())
spe <- loadHDF5SummarizedExperiment(here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_norm_QC_filtered"
))

#Load colData from 02 script  
res <- readRDS( here(
  "processed-data", "HD_Full_Analysis", "SPEs", "spe_cell_RCTD_colData.Rds"
))

res

#replace colData
stopifnot(identical(rownames(res),colnames(spe)))

colData(spe) <- res

## Get max weight and assigned type per cell
weight_cols <- grep("^rctd_weight_", names(colData(spe)), value = TRUE)
weight_mat <- as.matrix(colData(spe)[, weight_cols])
spe$rctd_max_type <- weight_cols[max.col(weight_mat)]
spe$rctd_max_type <- gsub("^rctd_weight_", "", spe$rctd_max_type)
table(spe$rctd_max_type)


#Revise coldata names
names(colData(spe))[match(weight_cols, names(colData(spe)))] <- gsub("^rctd_weight_", "", weight_cols)

#51 out of 1925144 cells seem to be filtered out during RCTD. Given that this is .0026% of the entire dataset, 
#I am going to remove them. 
spe <- spe[,!is.na(spe$rctd_max_type)]
spe

for(sample in unique(spe$sample_id)){
  print(sample)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","RCTD",sample)
  if(dir.exists(paths = dir_path)){
    unlink(dir_path,recursive = TRUE)
    dir.create(dir_path)
  }else{
    dir.create(dir_path)
  }
}

#Plot the weights for each cell type ojne the tissue section 
for(sample in unique(spe$sample_id)){
  print(sample)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","RCTD",sample)
  spe_sub <- spe[,spe$sample_id == sample]
  for(ct in unique(spe$rctd_max_type)){
    print(ct)
    p <- make_escheR(spe_sub) |>
      add_fill(ct) +
      scale_fill_gradientn(colors = c("lightgrey","red","black"))
    ggsave(filename = here(dir_path,paste0(ct,".png")),
           plot = p,
           width = 20, height = 14, dpi = 200)
  }
}

#Plot the max type
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE) 

for(sample in unique(spe$sample_id)){
  print(sample)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","RCTD",sample)
  spe_sub <- spe[,spe$sample_id == sample]
  p <- make_escheR(spe_sub) |>
    add_fill("rctd_max_type") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here(dir_path,paste0(sample,"all_celltypes",".png")),
         plot = p,
         width = 20, height = 14, dpi = 200)
}


#Plot one predicted cell type at a time
for(ct in unique(spe$rctd_max_type)){
  print(ct)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","RCTD",ct)
  if(dir.exists(paths = dir_path)){
    unlink(dir_path,recursive = TRUE)
    dir.create(dir_path)
  }else{
    dir.create(dir_path)
  }
}


for(ct in unique(spe$rctd_max_type)){
  print(ct)
  dir_path <- here("plots","HD_Full_Analysis","LabelTransfer","RCTD",ct)
  for(sample in unique(spe$sample_id)){
    print(sample)
    spe_sub <- spe[,spe$sample_id == sample]
    spe_sub$CellType <- ifelse(spe_sub$rctd_max_type == ct,
                               ct,
                               "Other")
    ct_colors <- c("lightgrey",cluster_cols[ct])
    names(ct_colors) <- c("Other",ct)
    p <- make_escheR(spe_sub) |>
      add_fill("CellType") +
      scale_fill_manual(values = ct_colors)
    ggsave(filename = here(dir_path,paste0(sample,".png")),
           plot = p,
           width = 20, height = 14, dpi = 200)
    
  }
}


sessionInfo()
