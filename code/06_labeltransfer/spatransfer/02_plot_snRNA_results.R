#Goal: Investigate spatransfer results
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

#Load libraries 
library(SingleCellExperiment)
library(SpatialExperiment)
library(nmfLabelTransfer)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)

#Load and prep the spe object containing the xenium NAc data
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))
logcounts(spe) <- assay(spe,"nucleus_normcounts")
rowData(spe)$gene_name <- rownames(spe)

#load in label transfer results
transfer_res <- readRDS(here("processed-data", "06_label_transfer", "snRNA_to_Xenium_target_predicitons.Rds"))

stopifnot(identical(rownames(colData(transfer_res$targets)),rownames(colData(spe))))

#Add predictions to the object
spe$nmf_preds <- as.character(transfer_res$targets$nmf_preds)

#load cluster colors from snRNA-seq paper. 
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)

#Remove neuron ambig
cluster_cols <- cluster_cols[-14]
cluster_cols

#Plot all cell types on the tissues
for(i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[,spe$Sample == i]
  p <- make_escheR(sub_spe) |>
    add_fill("nmf_preds") +
    scale_fill_manual(values = cluster_cols)
  ggsave(filename = here("plots","06_label_transfer","spaTransfer","snRNA",paste0(i,".png")),
         height = 20, width = 20)
  
}


#Plot each cell type on the tissue individually
for(sample in unique(spe$Sample)){
  print(sample)
  sub_spe <- spe[,spe$Sample == sample]
  dir.create(here("plots","06_label_transfer","spaTransfer",
                  "snRNA","PerSample_PerCluster",sample))
  for(celltype in unique(spe$nmf_preds)){
    print(celltype)
    sub_spe$pred_celltype_only <- ifelse(sub_spe$nmf_preds == celltype,
                                         celltype,
                                         "Other")
    celltype_colors <- c(cluster_cols[celltype],"grey")
    names(celltype_colors)[2] <- "Other"
    p <- make_escheR(sub_spe) |>
      add_fill("pred_celltype_only") +
      scale_fill_manual(values = celltype_colors)
    ggsave(filename = here("plots","06_label_transfer","spaTransfer",
                           "snRNA","PerSample_PerCluster",
                           sample,
                           paste0(celltype,"_",sample,".png")),
           height = 20, width = 20)
    
  }
}

## Reproducibility information
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessionInfo()
