#module load r_nac
#cd cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/


library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

sce
# class: SingleCellExperiment 
# dim: 36601 103785 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(36601): ENSG00000243485 ENSG00000237613 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(103785): 1_AAACCCAAGACCAACG-1 1_AAACCCACAGTCAGCC-1 ...
# 20_TTTGTTGCAAGATGTA-1 20_TTTGTTGGTACGAAAT-1
# colData names(41): Sample Barcode ... sizeFactor CellType.Final
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#Remove the Neuron_Ambig group
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

dim(sce)
#[1]  36601 103339

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)
# FALSE  TRUE 
# 34977  1624 

#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

dim(sce)
#[1]  34977 103339

#For this each sample needs an anterior/middle/posterior designation. 
#Combine the Anterior and Posterior samples into a single group called "Anterior_Posterior"
#Make a dataframe. 
Ant_Mid_Post <- data.frame(Brain_ID = unique(sce$Brain_ID))

Ant_Mid_Post
# Brain_ID
# 1    Br8325
# 2    Br8492
# 3    Br2720
# 4    Br6423
# 5    Br2743
# 6    Br3942
# 7    Br6432
# 8    Br6471
# 9    Br6522
# 10   Br8667

Ant_Mid_Post <- cbind(Ant_Mid_Post,c("Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior",
                                     "Anterior_Posterior","Anterior_Posterior",
                                     "Anterior_Posterior","Middle",
                                     "Middle","Anterior_Posterior"))

colnames(Ant_Mid_Post)[2] <- "Depth"

Ant_Mid_Post
# Brain_ID              Depth
# 1    Br8325 Anterior_Posterior
# 2    Br8492             Middle
# 3    Br2720             Middle
# 4    Br6423 Anterior_Posterior
# 5    Br2743 Anterior_Posterior
# 6    Br3942 Anterior_Posterior
# 7    Br6432 Anterior_Posterior
# 8    Br6471             Middle
# 9    Br6522             Middle
# 10   Br8667 Anterior_Posterior

#Add depth to the sce object
sce$Depth <- Ant_Mid_Post[match(sce$Brain_ID,Ant_Mid_Post$Brain_ID),"Depth"]

as.data.frame(unique(colData(sce)[,c("Brain_ID","Depth")]))
# Brain_ID              Depth
# 1_AAACCCAAGACCAACG-1    Br8325 Anterior_Posterior
# 3_AAACCCAAGGTGAGCT-1    Br8492             Middle
# 5_AAACCCAGTAATTAGG-1    Br2720             Middle
# 7_AAACCCACACCCTTAC-1    Br6423 Anterior_Posterior
# 9_AAACCCACATTGTAGC-1    Br2743 Anterior_Posterior
# 11_AAACCCAAGACTCCGC-1   Br3942 Anterior_Posterior
# 13_AAACCCAAGGACTGGT-1   Br6432 Anterior_Posterior
# 15_AAACCCAAGGGAGATA-1   Br6471             Middle
# 17_AAACCCAAGATTGCGG-1   Br6522             Middle
# 19_AAACCCACAAGGCCTC-1   Br8667 Anterior_Posterior

#Now subset for DRD1_MSN_B samples only 
sce <- sce[,sce$CellType.Final == "DRD1_MSN_B"]

sce
# class: SingleCellExperiment 
# dim: 34977 6544 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(34977): ENSG00000243485 ENSG00000186092 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames(6544): 1_AAAGGATAGCTCCACG-1 1_AAATGGAAGACGCCCT-1 ...
# 20_TTTCCTCGTGGAAGTC-1 20_TTTGATCGTCGACGCT-1
# colData names(42): Sample Barcode ... CellType.Final Depth
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#Pseudobulk across CellType and Brain_ID
sce_pb <- aggregateAcrossCells(sce,ids = colData(sce)[,c("CellType.Final","Brain_ID")])

sce_pb
# class: SingleCellExperiment 
# dim: 34977 10 
# metadata(1): Samples
# assays(1): counts
# rownames(34977): ENSG00000243485 ENSG00000186092 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames: NULL
# colData names(45): Sample Barcode ... Brain_ID ncells
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

table(sce_pb$CellType.Final)
# DRD1_MSN_B 
# 10 

table(sce_pb$Brain_ID)
# Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667 
# 1      1      1      1      1      1      1      1      1      1 

#Are there any genes that have 0 counts for every gene within the DRD1_MSN_B cluster only. 
table(rowSums(assay(sce_pb, "counts")) == 0)
# FALSE  TRUE 
# 31968  3009 
#Yes let's remove those.

sce_pb <- sce_pb[!rowSums(assay(sce_pb, "counts")) == 0, ]


sce_pb
# class: SingleCellExperiment 
# dim: 31968 10 
# metadata(1): Samples
# assays(1): counts
# rownames(31968): ENSG00000243485 ENSG00000238009 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames: NULL
# colData names(45): Sample Barcode ... Brain_ID ncells
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):
