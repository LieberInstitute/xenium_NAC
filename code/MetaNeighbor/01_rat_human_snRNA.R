# Use metaneighbor to compare human and rat snRNA-seq
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(MetaNeighbor)
library(orthogene)
library(Seurat)
library(here)

rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

#subset for just neurons
rat_sce <- rat_sce[,rat_sce$ident %in% c("Drd1.1","Drd1.2","Drd2.1","Drd2.2","Drd3","GABA","Chst9",
                                         "Sema5a","Pvalb","Sst","ChAT")]

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

#Load human sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#Remove the neuronal ambiguous population 
sce <- sce[,sce$CellType.Final %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                      "DRD2_MSN_A","DRD2_MSN_B",
                                      "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                      "Excitatory")]

sce

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq

#subset for 121 orthos
sce_sub <- sce[rownames(sce) %in% rownames(rat_sce),]

sce_sub

#subset rat object for genes in the spe
rat_sce <- rat_sce[rownames(rat_sce) %in% rownames(sce_sub),]

#Put rat sce in same order as spe_sub
rat_sce <- rat_sce[rownames(sce_sub),]

identical(rownames(sce_sub),rownames(rat_sce))
stopifnot(identical(rownames(sce_sub),rownames(rat_sce)))

#Add a gene name and gene_id column
rowData(rat_sce)[,c("gene_name","gene_id")] <- rowData(sce_sub)[,c("gene_name","gene_id")]

#Generate a single object containing both rat and human data. 
#Add column data with dataset information.
rat_sce$study_id <- "rat_snRNA"
sce_sub$study_id <- "human_snRNA"
rat_sce$species <- "rat"
sce_sub$species <- "human"

#Remove rowRanges to avoid any issues regarding differences in seqnames between species/objects.
#Before though, save the rownaems to reset after.
genes <- rownames(sce_sub)
stopifnot(identical(genes,rownames(rat_sce)))
rowRanges(rat_sce) <- NULL
rowRanges(sce_sub) <- NULL

#Now subset the column data down to sample and cell type
colData(rat_sce) <- colData(rat_sce)[,c("sample","cellType","study_id","species")]
colnames(colData(rat_sce)) <- c("Sample","CellType.Final","study_id","species")
colData(sce_sub) <- colData(sce_sub)[,c("Sample","CellType.Final","study_id","species")]
sce_sub
rat_sce
#   Combine 
stopifnot(identical(rownames(sce_sub),rownames(rat_sce)))
combo <- cbind(rat_sce, sce_sub)
rownames(combo) <- genes

combo

#   Create a cell type label column that includes dataset of origin,
#   so MetaNeighbor can distinguish rat MSN from human MSN
combo$celltype_label <- paste(
  combo$study_id,
  combo$CellType.Final,
  sep = "|"
  )

#   Pick highly variable genes from the shared ortholog set for the
#   similarity network (MetaNeighbor recommends this over using all genes)
hvgs <- variableGenes(dat = combo, 
                      exp_labels = combo$species)

#Keep top 2000 HVGs
hvgs_keep <- hvgs[1:4000]

#   Run unsupervised MetaNeighbor
aurocs <- MetaNeighborUS(
  var_genes = hvgs_keep,
  dat = combo,
  study_id = combo$species,
  cell_type = combo$CellType.Final,
  fast_version = TRUE
)

#   Plot the AUROC heatmap
pdf(file = here("plots","MetaNeighbor","All_Neurons_both_species.pdf"),width = 10, height = 12)
plotHeatmap(aurocs,cex = .7)
dev.off()

#Just cross-species
rat_types <- grep("^rat", rownames(aurocs), value = TRUE)
human_types <- grep("^human", rownames(aurocs), value = TRUE)

cross_species <- aurocs[rat_types, human_types]
hm1 <- pheatmap::pheatmap(
  cross_species,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  display_numbers = TRUE,
  number_format = "%.2f",
  fontsize_number = 7,
  main = "MetaNeighbor AUROC: Rat vs Human"
)

pdf(file = here("plots","MetaNeighbor","All_Neurons_cross_species.pdf"),width = 10, height = 12)
grid::grid.newpage()
grid::grid.draw(hm1$gtable)
dev.off()

sessionInfo()
