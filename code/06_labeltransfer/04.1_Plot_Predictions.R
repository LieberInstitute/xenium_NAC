#Add predicted id information to 6436 spe object
library(SpatialExperiment)
library(escheR)
library(Seurat)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

#Subset to Br6436
spe_6436 <- spe[,spe$Donor == "Br6436"]

#Load suerat object for Br6436 containing the 
seurat_6436 <- readRDS(here("processed-data", "06_label_transfer","Objects", "seurat_6436_predictions_object.Rds"))

#Are the objects in the same order? 
stopifnot(identical(colnames(spe_6436),colnames(seurat_6436)))
stopifnot(identical(rownames(colData(spe_6436)),rownames(seurat_6436@meta.data)))

#Add the predicted id to the column data
spe_6436$CellType <- seurat_6436$predicted.id

#Load cell type colors
CellType_cols <- readRDS(here("processed-data","05_Clustering","CellType_cols.Rds"))

##########
#plot the clusters on the tissue
for(i in unique(spe_6436$Sample)){
  print(i)
  sub_spe <- spe_6436[,spe_6436$Sample == i]
  p <- make_escheR(sub_spe) |>
    add_fill("CellType") +
    scale_fill_manual(values = CellType_cols)
  ggsave(filename = here("plots","06_label_transfer","Predicted_ID_6436",paste0(i,".png")),
         height = 20, width = 20)
}

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
