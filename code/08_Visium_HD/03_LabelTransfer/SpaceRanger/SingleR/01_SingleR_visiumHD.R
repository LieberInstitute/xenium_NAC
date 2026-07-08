#Goal: Perform label transfer with SingleR on a per sample basis. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5

library(SingleCellExperiment)
library(SpatialExperiment)
library(BiocParallel)
library(sessioninfo)
library(HDF5Array)
library(SingleR)
library(ggplot2)
library(escheR)
library(scater)
library(scran)
library(here)

## get sample i
sample_i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

message("Loading sfe object - ", Sys.time())
sfe_dir <- here("processed-data", "HD_Full_Analysis", "sfe_cell_norm_QC_filtered_v2") 

sfe <- loadHDF5SummarizedExperiment(sfe_dir)
sfe

#Calculate log-normalized counts
sfe <- computeLibraryFactors(sfe)
sfe <- logNormCounts(sfe)


#load single cell rds object
message("Loading sce object - ", Sys.time())
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

sce

#Remove the Neuron_Ambiguous cluster of cells. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

#Subset the sfe for a single sample. 
samples <- unique(sfe$sample_id)
sample_run <- samples[[sample_i]]
message("Running Sample: ", sample_run, " (", sample_i, "/", length(samples), ")")

#Subset for the sample
sfe <- sfe[, sfe$sample_id == sample_run]
message("Number of cells in this sample:", ncol(sfe))

#create directory to stick plots in
dir <- file.path(here("plots","HD_Full_Analysis","LabelTransfer","SpaceRanger","SingleR",sample_run))
dir.create(dir)

#set gene symbols as feature names
rownames(sce) <- rowData(sce)$gene_name


#Run label transfer 
res <- SingleR(test = sfe,
               ref = sce,
               labels = sce$CellType.Final,
               assay.type.test = "logcounts",
               assay.type.ref = "logcounts",
               aggr.ref = TRUE) 


#Save the results
saveRDS(res,file = here("processed-data","HD_Full_Analysis","LabelTransfer","SpaceRanger",
                        "SingleR",paste0(sample_run,"_SingleR_results.Rds")))

#Add results to the object
stopifnot(identical(rownames(res),colnames(sfe)))

#Add cell type info
sfe$SingleR_labels <- res$labels
sfe$SingleR_labels_pruned <- res$pruned.labels

#Remove the NA level
labels_pruned <- unique(res$pruned.labels)
labels_pruned <- labels_pruned[!is.na(labels_pruned)]

#Generate spot plots
for(l in labels_pruned){
  print(l)
  
  #Make 1 vs other column
  sfe$CellType_of_Interest <- ifelse(sfe$SingleR_labels_pruned == l,
                                     l,
                                     "Other")
  colors <- c("red","lightgrey")
  names(colors) <- c(l,"Other")
  
  #Generate escheR plot
  p <- make_escheR(sfe,spot_size = 0.5) |>
    add_fill(var = "CellType_of_Interest") +
    ggtitle(l) +
    scale_fill_manual(values = colors) +
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
