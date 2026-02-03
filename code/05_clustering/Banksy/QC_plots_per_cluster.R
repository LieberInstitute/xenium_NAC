# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SpatialExperiment)
library(here)

#Load the xenium object 
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe


#Add the Banksy cell types to the domain. 
#Load the annotated banksy cell types/u tgufl
Banksy_celltypes <- readRDS(here("processed-data","06_label_transfer","Objects","Banksy_CellTypes_transfer.Rds"))

stopifnot(identical(Banksy_celltypes$cell_id,colnames(spe)))

spe$Banksy_celltypes <- Banksy_celltypes$CellType

#Plot some QC metrics of the clusters
#First factorize the clusters
spe$Banksy_celltypes <- factor(x = spe$Banksy_celltypes,
                               levels = c("D1_Island_A","D1_Island_B","DRD1_MSN","DRD2_MSN","CHAT",
                                          "Inh_PVALB","Inh_SST","Excitatrory","Astro_A","Astro_B",
                                          "Ependymal","Fibroblast_A","Fibroblast_B","WM_A","WM_B",
                                          "WM_C","WM_D","WM_E","OPC","Microglia_A",
                                          "Microglia_B","MSN_Oligo")) 

#Load cluster colors 
CellType_cols <- readRDS(here("processed-data","05_Clustering","CellType_cols.Rds"))

library(scater)
#Number of genes
detected_vln <- plotColData(spe,x = "Banksy_celltypes",
            y = "detected",colour_by = "Banksy_celltypes") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_color_manual(values = CellType_cols) +
  scale_y_log10() +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) 
ggsave(filename = here("plots","05_clustering","Banksy","Detected_Banksy_Annotated.png"),plot = detected_vln)

#library size
library_vln <- plotColData(spe,x = "Banksy_celltypes",
                            y = "sum",colour_by = "Banksy_celltypes") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  scale_color_manual(values = CellType_cols) +
  scale_y_log10() +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3) 
ggsave(filename = here("plots","05_clustering","Banksy","Sum_Banksy_Annotated.png"),plot = library_vln)

sessionInfo()
