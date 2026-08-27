# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialFeatureExperiment)
library(SpatialExperiment)
library(HDF5Array)
library(here)

here::i_am("code/13_app/00_write_objects.R")

#####Xenium 
xen_spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_celltype_v2.Rds"))

xen_spe

#Before doing any size reduction, save the rowData to a txt file to be the supp file for the xenium panel
genes <- as.data.frame(rowData(xen_spe))
genes <- genes[,c("ID","Symbol")]
rownames(genes) <- NULL
colnames(genes) <- c("ENSEMBL_ID","Gene_Name")
write.csv(x = genes,file = here("processed-data","Table_S1.csv"))

#remove cellular normalized counts
assay(xen_spe, "cell_normcounts") <- NULL
logcounts(xen_spe) <- assay(xen_spe,"nucleus_normcounts")
assay(xen_spe, "nucleus_normcounts") <- NULL

#Keep ID and Sybol 
rowData(xen_spe) <- rowData(xen_spe)[,c("ID","Symbol")]
colData(xen_spe) <- colData(xen_spe)[,c(1:13,35:61)]

saveRDS(xen_spe,file = here("processed-data","13_app","Xenium_spe.Rds"))


##### Visium-HD 
sfe_path <- here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
sfe <- loadHDF5SummarizedExperiment(sfe_path)

rowData(sfe) <- rowData(sfe)[,c("ID","Symbol")]
colData(sfe) <- colData(sfe)[,c(1:17,130:137)]
reducedDim(sfe,"NMF_Proj") <- NULL


saveHDF5SummarizedExperiment(sfe,here("processed-data","13_app","VisiumHD_sfe"))

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

