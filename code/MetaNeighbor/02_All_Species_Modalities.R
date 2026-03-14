# Use metaneighbor to compare human, NHP, rat snRNA-seq
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(MetaNeighbor)
library(orthogene)
library(Seurat)
library(here)

#############################################
###########     Human   #####################
#############################################

message("Loading human objects |", Sys.time())

##### Xenium
# Load human spe object (Xenium)
spe <- readRDS(here("processed-data", "02_build_spe", "SPEs", "spe_celltype_v2.Rds"))


###### Single nucleus RNA-sequencing object
# Load human sce object (snRNA-seq)
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

# Remove ambiguous neurons and genes with zero counts
sce <- sce[, sce$CellType.Final != "Neuron_Ambig"]
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

# Uniquify gene names to avoid duplicates
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(
  rowData(sce)$gene_id, rowData(sce)$gene_name
)
rownames(sce) <- rowData(sce)$Symbol.uniq

##### Find intersection of genes and subset both objects #####
shared_genes <- intersect(rownames(sce), rownames(spe))


###########################################
###########     RAT   #####################
###########################################

message("Loading Rat object |", Sys.time())

#Plan: Identify 1-to-1 orthologs from rat to human --> Subset the rat object --> Change rownames

rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

#subset for just neurons
rat_sce <- rat_sce[,rat_sce$ident %in% c("Drd1.1","Drd1.2","Drd2.1","Drd2.2","Drd3","GABA","Chst9",
                                         "Sema5a","Pvalb","Sst","ChAT")]

#Subset to just neurons. 


#Convert to 121 orthologs
rat_orthologs <- convert_orthologs(gene_df = counts(rat_sce),
                                   gene_input = "rownames",
                                   gene_output = "dict",
                                   input_species = "rat",
                                   output_species = "human",
                                   non121_strategy = "drop_both_species",
                                   method = "gprofiler")
head(rat_orthologs)

#Generate a dataframe from the above output
rat_ortho_df <- data.frame(rat = names(rat_orthologs),
                           human = unname(rat_orthologs))

#Subset rat_sce to contain only the 121 orthologs
rat_sce <- rat_sce[rownames(rat_sce) %in% rat_ortho_df$rat,]

#Put the rat ortholog datafarme in same order as rat_sce
rat_ortho_df <- rat_ortho_df[match(rownames(rat_sce),rat_ortho_df$rat),]

identical(rownames(rat_sce),rat_ortho_df$rat)

stopifnot(identical(rownames(rat_sce),rat_ortho_df$rat))

#convert the rownames of the rat_sce object to the human orthologs
rownames(rat_sce) <- rat_ortho_df$human

###########################################
###########     NHP   #####################
###########################################

message("Loading NHP object |", Sys.time())

#Plan: Identify 1-to-1 orthologs from NHP to human --> Subset the  NHP object --> Change rownames

nhp_obj <- readRDS(file = here("processed-data","NHP_Data","Results_MSNs_processed_final.rds"))
DefaultAssay(nhp_obj) <- "RNA"
nhp_sce <- as.SingleCellExperiment(x = nhp_obj)
nhp_sce

#Convert to 121 orthologs
nhp_orthologs <- convert_orthologs(gene_df = counts(nhp_sce),
                                   gene_input = "rownames",
                                   gene_output = "dict",
                                   input_species = "macaque",
                                   output_species = "human",
                                   non121_strategy = "drop_both_species",
                                   method = "gprofiler")
head(nhp_orthologs)

#Generate a dataframe from the above output
nhp_ortho_df <- data.frame(NHP = names(nhp_orthologs),
                           human = unname(nhp_orthologs))

nhp_sce <- nhp_sce[rownames(nhp_sce) %in% nhp_ortho_df$NHP,]

nhp_ortho_df <- nhp_ortho_df[match(rownames(nhp_sce),nhp_ortho_df$NHP),]

identical(rownames(nhp_sce),nhp_ortho_df$NHP)

stopifnot(identical(rownames(nhp_sce),nhp_ortho_df$NHP))

#convert the rownames of the rat_sce object to the human orthologs
rownames(nhp_sce) <- nhp_ortho_df$human

#####Find intersection of the genes and subset all objects

message("Set up combo object |", Sys.time())

#Find the intersection of the mouse, human, and rat
rat_nhp_orthos <- intersect(rownames(rat_sce),rownames(nhp_sce))
all_orthos <- intersect(rat_nhp_orthos,shared_genes)
length(all_orthos) 


#now subset
#rat
rat_sub <- rat_sce[all_orthos,]
rat_sub

#nhp
nhp_sub <- nhp_sce[all_orthos,]
nhp_sub

#Subset the human objects
sce <- sce[all_orthos, ]
spe <- spe[all_orthos, ]

# Keep counts and rename nucleus_normcounts -> logcounts in Xenium
assays(spe) <- assays(spe)[c("counts", "nucleus_normcounts")]
assayNames(spe)[2] <- "logcounts"

###############Generate a single object containing both rat and human data. 
#Add column data with dataset information.
sce$dataset <- "Human_snRNA"
spe$dataset <- "Human_Xenium"
nhp_sub$dataset <- "NHP_snRNA"
rat_sub$dataset <- "Rat_snRNA"

#Remove rowRanges to avoid any issues regarding differences in seqnames between species/objects.
#Before though, save the rownaems to reset after.
rowRanges(sce) <- NULL
rowRanges(spe) <- NULL
rowRanges(rat_sub) <- NULL
rowRanges(nhp_sub) <- NULL

#Now subset the column data down to sample and cell type
#Rat 
colData(rat_sub) <- colData(rat_sub)[,c("sample","cellType","dataset")]
colnames(colData(rat_sub)) <- c("Sample","CellType","dataset")
rat_sub

#NHP
colData(nhp_sub) <- colData(nhp_sub)[,c("monkey","MSN_type","dataset")]
colnames(colData(nhp_sub)) <- c("Sample","CellType","dataset")
nhp_sub

#Human
colData(sce) <- colData(sce)[,c("Sample","CellType.Final","dataset")]
colnames(colData(sce)) <- c("Sample","CellType","dataset")
sce


# Xenium (spatial)
colData(spe) <- colData(spe)[, c("Sample", "CellTypes", "dataset")]
colnames(colData(spe)) <- c("Sample", "CellType", "dataset") #SpatialExperiment forces a sample_id column

#Human xenium 
spe_clean <- SingleCellExperiment(
  assays = assays(spe),
  colData = colData(spe)[,c("Sample","CellType","dataset")]
)
rownames(spe_clean) <- rownames(spe)

#   Combine 
stopifnot(identical(rownames(sce),rownames(nhp_sub)))
stopifnot(identical(rownames(rat_sub),rownames(nhp_sub)))
stopifnot(identical(rownames(sce),rownames(rat_sub)))
stopifnot(identical(rownames(spe_clean),rownames(nhp_sub)))
stopifnot(identical(rownames(spe_clean),rownames(rat_sub)))
stopifnot(identical(rownames(spe_clean),rownames(sce)))

#Combine all objects
combo <- cbind(rat_sub,nhp_sub,sce,spe_clean)
rownames(combo) <- all_orthos

combo

message("Combo object finished |", Sys.time())
message("Calculating HVGs |", Sys.time())
#   Pick highly variable genes from the shared ortholog set for the
#   similarity network 
hvgs <- variableGenes(dat = combo, 
                      exp_labels = combo$dataset)

#Keep top 4000 HVGs
hvgs_keep <- hvgs[1:4000]

message("Running MetaNeighbor |", Sys.time())
#   Run unsupervised MetaNeighbor
aurocs <- MetaNeighborUS(
  var_genes = hvgs_keep,
  dat = combo,
  study_id = combo$dataset,
  cell_type = combo$CellType,
  fast_version = TRUE,
  one_vs_best = TRUE, symmetric_output = FALSE
)
message("Finished MetaNeighbor |", Sys.time())

#Save the output
saveRDS(aurocs,file = here("processed-data","MetaNeighbor","All_Species_Modalities_aurocs.Rds"))

sessionInfo()
