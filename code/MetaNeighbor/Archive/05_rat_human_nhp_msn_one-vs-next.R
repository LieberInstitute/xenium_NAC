# Use metaneighbor to compare human, NHP, rat snRNA-seq
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(MetaNeighbor)
library(orthogene)
library(Seurat)
library(here)

#The NHP object only contains MSNs. Subset all objects to just neurons. 

###########################################
###########     RAT   #####################
###########################################

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

#############################################
###########     Human   #####################
#############################################
#Load human sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

#subset for neurons
sce <- sce[,sce$CellType.Final %in% c("DRD1_MSN_A","DRD1_MSN_B","DRD1_MSN_C","DRD1_MSN_D",
                                      "DRD2_MSN_A","DRD2_MSN_B",
                                      "Inh_A","Inh_B","Inh_C","Inh_D","Inh_E","Inh_F",
                                      "Excitatory")]

sce

# -> To avoid getting into sticky situations, let's 'uniquify' these names:
rowData(sce)$Symbol.uniq <- scuttle::uniquifyFeatureNames(rowData(sce)$gene_id, rowData(sce)$gene_name)
rownames(sce) <- rowData(sce)$Symbol.uniq



#####Find intersection of the genes and subset all objects

#Find the intersection of the mouse, human, and rat
rat_nhp_orthos <- intersect(rownames(rat_sce),rownames(nhp_sce))
all_orthos <- intersect(rat_nhp_orthos,rownames(sce))
length(all_orthos) 

#now subset
#rat
rat_sub <- rat_sce[all_orthos,]
rat_sub

#nhp
nhp_sub <- nhp_sce[all_orthos,]
nhp_sub

#Human
human_sub <- sce[all_orthos,]
human_sub

###############Generate a single object containing both rat and human data. 
#Add column data with dataset information.
human_sub$species <- "Human"
nhp_sub$species <- "NHP"
rat_sub$species <- "Rat"

#Remove rowRanges to avoid any issues regarding differences in seqnames between species/objects.
#Before though, save the rownaems to reset after.
rowRanges(human_sub) <- NULL
rowRanges(rat_sub) <- NULL
rowRanges(nhp_sub) <- NULL

#Now subset the column data down to sample and cell type
#Rat 
colData(rat_sub) <- colData(rat_sub)[,c("sample","cellType","species")]
colnames(colData(rat_sub)) <- c("Sample","CellType","species")
rat_sub

#NHP
colData(nhp_sub) <- colData(nhp_sub)[,c("monkey","MSN_type","species")]
colnames(colData(nhp_sub)) <- c("Sample","CellType","species")
nhp_sub

#Human
colData(human_sub) <- colData(human_sub)[,c("Sample","CellType.Final","species")]
colnames(colData(human_sub)) <- c("Sample","CellType","species")
human_sub

#   Combine 
stopifnot(identical(rownames(human_sub),rownames(nhp_sub)))
stopifnot(identical(rownames(rat_sub),rownames(nhp_sub)))
stopifnot(identical(rownames(human_sub),rownames(rat_sub)))

#Combine all objects
combo <- cbind(rat_sub,nhp_sub,human_sub)
rownames(combo) <- all_orthos

combo

#   Pick highly variable genes from the shared ortholog set for the
#   similarity network 
hvgs <- variableGenes(dat = combo, 
                      exp_labels = combo$species)

#Keep top 4000 HVGs
hvgs_keep <- hvgs[1:4000]

#   Run unsupervised MetaNeighbor
aurocs <- MetaNeighborUS(
  var_genes = hvgs_keep,
  dat = combo,
  study_id = combo$species,
  cell_type = combo$CellType,
  fast_version = TRUE,
  one_vs_best = TRUE, symmetric_output = FALSE
)


# Parse species from row/column names
species_col <- gsub("\\|.*", "", colnames(aurocs))
species_row <- gsub("\\|.*", "", rownames(aurocs))

# Map species to colors
species_colors <- c("Human" = "sienna2", "NHP" = "darkgreen", "Rat" = "dodgerblue3")
col_side <- species_colors[species_col]
row_side <- species_colors[species_row]

pdf(file = here("plots","MetaNeighbor",
                "Human_NHP_Rat_Str_MetaNeighbor_Heatmap.pdf"),
    width = 12, height = 12)

plotHeatmap(aurocs, 
            cex = 0.75,
            ColSideColors = col_side,
            RowSideColors = row_side)

# Add a legend for the species colors
legend("topright", 
       legend = names(species_colors), 
       fill = species_colors, 
       title = "Species", 
       border = NA,
       cex = 0.8)

dev.off()

sessionInfo()
