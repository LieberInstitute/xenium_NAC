#Goal: Perform label transfer with SingleR on a per sample basis. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.4.x

library(SingleCellExperiment)
library(SpatialExperiment)
library(BiocParallel)
library(sessioninfo)
library(SingleR)
library(ggplot2)
library(escheR)
library(scater)
library(scran)
library(here)

## get sample i
sample_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe

#load single cell rds object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

sce

#Subset the spe for a single sample. 
samples <- unique(spe$Sample)
sample_run <- samples[[sample_i]]
message("Running Sample: ", sample_run, " (", sample_i, "/", length(samples), ")")

spe <- spe[, spe$Sample == sample_run]
message("Number of cells in this sample:", ncol(spe))

#create directory to stick plots in
dir <- file.path(here("plots","06_label_transfer",sample_run))
dir.create(dir)

#Subset to genes only in the xenium dataset
sce_sub <- sce[rowData(sce)$gene_id %in% rowData(spe)$ID,]

dim(sce_sub)

#Remove the Neuron_Ambiguous cluster of cells. 
sce_sub <- sce_sub[,sce_sub$CellType.Final != "Neuron_Ambig"]

#set gene symbols as feature names
rownames(sce_sub) <- rowData(sce_sub)$gene_name


#Run label transfer 
res <- SingleR(test = spe,
               ref = sce_sub,
               labels = sce_sub$CellType.Final,
               assay.type.test = "nucleus_normcounts",
               assay.type.ref = "logcounts",
               aggr.ref = TRUE) #Pseudobulks reference data. 


#Save the results
save(res,file = here("processed-data","06_label_transfer",paste0(sample_run,"_SingleR_results.rda")))

#Add results to the object
stopifnot(identical(rownames(res),colnames(spe)))

#Add cell type info
spe$SingleR_labels <- res$labels
spe$SingleR_labels_pruned <- res$pruned.labels

#Remove the NA level
labels_pruned <- unique(res$pruned.labels)
labels_pruned <- labels_pruned[!is.na(labels_pruned)]

#Generate spot plots
for(l in labels_pruned){
  print(l)
  
  #Make 1 vs other column
  spe$CellType_of_Interest <- ifelse(spe$SingleR_labels_pruned == l,
                                     l,
                                     "Other")

  #Generate escheR plot
  p <- make_escheR(spe,spot_size = 0.5) |>
    add_fill(var = "CellType_of_Interest") +
    ggtitle(l) +
    theme(plot.title = element_text(hjust = 0.5))
  
  
  #Save plot 
  ggsave(plot = p,
         filename = file.path(paste0(dir,"/",l,"_",sample_run,".png")),
         height = 20,
         width = 20)
}



###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()


