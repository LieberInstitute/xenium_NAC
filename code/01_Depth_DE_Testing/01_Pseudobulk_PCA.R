#module load r_nac
#cd /dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/
#Goal: Perform pseudobulking and PCA analysis across the Samples that are labeled Anterior, middle, posterior. 


library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(scater)
library(scran)
library(here)

#load the sce object
sce <- readRDS(here("processed-data","12_snRNA","sce_CellType_noresiduals.Rds"))

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

Ant_Mid_Post <- cbind(Ant_Mid_Post,c("Posterior","Middle","Middle","Anterior","Anterior","Posterior","Anterior","Middle","Middle","Posterior"))

colnames(Ant_Mid_Post)[2] <- "Depth"

Ant_Mid_Post
# Brain_ID     Depth
# 1    Br8325 Posterior
# 2    Br8492    Middle
# 3    Br2720    Middle
# 4    Br6423  Anterior
# 5    Br2743  Anterior
# 6    Br3942 Posterior
# 7    Br6432  Anterior
# 8    Br6471    Middle
# 9    Br6522    Middle
# 10   Br8667 Posterior

#Add depth to the sce object
sce$Depth <- Ant_Mid_Post[match(sce$Brain_ID,Ant_Mid_Post$Brain_ID),"Depth"]

as.data.frame(unique(colData(sce)[,c("Brain_ID","Depth")]))
#                       Brain_ID     Depth
# 1_AAACCCAAGACCAACG-1    Br8325 Posterior
# 3_AAACCCAAGGTGAGCT-1    Br8492    Middle
# 5_AAACCCAGTAATTAGG-1    Br2720    Middle
# 7_AAACCCACACCCTTAC-1    Br6423  Anterior
# 9_AAACCCACATTGTAGC-1    Br2743  Anterior
# 11_AAACCCAAGACTCCGC-1   Br3942 Posterior
# 13_AAACCCAAGGACTGGT-1   Br6432  Anterior
# 15_AAACCCAAGGGAGATA-1   Br6471    Middle
# 17_AAACCCAAGATTGCGG-1   Br6522    Middle
# 19_AAACCCACAAGGCCTC-1   Br8667 Posterior

#For pseudobulk DEG testing, you need to test 1 condition vs another. Can't have three designations. Make three columns that
#will be anterior,middle,posterior and the other samples will just say Other. 
sce$Depth_Anterior <- ifelse(sce$Depth == "Anterior",
                             "Anterior",
                             "Other")

sce$Depth_Middle <- ifelse(sce$Depth == "Middle",
                           "Middle",
                           "Other")

sce$Depth_Posterior <- ifelse(sce$Depth == "Posterior",
                              "Posterior",
                              "Other")

#Pseudobulk across CellType and Brain_ID
sce_pb <- aggregateAcrossCells(sce,ids = colData(sce)[,c("CellType.Final","Brain_ID")])


sce_pb
# class: SingleCellExperiment 
# dim: 34977 199 
# metadata(1): Samples
# assays(1): counts
# rownames(34977): ENSG00000243485 ENSG00000186092 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames: NULL
# colData names(48): Sample Barcode ... Brain_ID ncells
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):

#One of the samples is not present in the pseudobulked object

table(sce_pb$CellType.Final)
# Astrocyte_A Astrocyte_B  DRD1_MSN_A  DRD1_MSN_B  DRD1_MSN_C  DRD1_MSN_D 
# 10          10          10          10          10          10 
# DRD2_MSN_A  DRD2_MSN_B Endothelial   Ependymal  Excitatory       Inh_A 
# 10          10          10           9          10          10 
# Inh_B       Inh_C       Inh_D       Inh_E       Inh_F   Microglia 
# 10          10          10          10          10          10 
# Oligo         OPC 
# 10          10 

table(sce_pb$Brain_ID)
# Br2720 Br2743 Br3942 Br6423 Br6432 Br6471 Br6522 Br8325 Br8492 Br8667 
# 20     20     20     20     20     20     20     20     19     20 

#8492 is not represented within the Ependymal population. That's okay. 

###############PCA analysis#####################
#Get library size factors. This function generates "per-cell size factors from library sizes (i.e., total sum of counts per cell)
sce_pb <- computeLibraryFactors(sce_pb)

#Generate log-normalized counts
sce_pb <- logNormCounts(sce_pb)

sce_pb
# class: SingleCellExperiment 
# dim: 34977 199 
# metadata(1): Samples
# assays(2): counts logcounts
# rownames(34977): ENSG00000243485 ENSG00000186092 ... ENSG00000278817
# ENSG00000277196
# rowData names(7): source type ... gene_type binomial_deviance
# colnames: NULL
# colData names(48): Sample Barcode ... Brain_ID ncells
# reducedDimNames(4): GLMPCA_approx tSNE HARMONY tSNE_HARMONY
# mainExpName: NULL
# altExpNames(0):
#logcounts are in the object now. 

#Pull log counts matrix
pb_log_counts <- assay(sce_pb,"logcounts")

#Perform PCA analysis without scaling
pca_no_scale <- prcomp(t(pb_log_counts), center = TRUE, scale. = FALSE)

#Pull PCA results
pca_scores <- as.data.frame(pca_no_scale$x)

#Check out the matrix. 
pca_scores[1:5,1:5]
# PC1        PC2       PC3      PC4      PC5
# 1 -145.7752 -109.21437 -34.36973 86.61906 27.85574
# 2 -162.3320  -99.15418 -17.20643 54.78211 16.94541
# 3 -170.3131  -96.11557 -15.91427 43.95342 16.18650
# 4 -144.8057 -110.36204 -37.43378 82.35241 22.25010
# 5 -159.1286  -99.06644 -15.00062 56.90941 16.61717

#Add metadata to the pca output
pca_scores$Depth <- colData(sce_pb)$Depth
pca_scores$CellType.Final <- colData(sce_pb)$CellType.Final
pca_scores$Brain_ID <- colData(sce_pb)$Brain_ID


#Plot the PCA colored by:
#CellType.Final
pb_celltype_pca <- ggplot(pca_scores,aes(x = PC1, y = PC2, color = CellType.Final)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_celltype_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_CellType.Final_PCA.pdf",
       height = 8,
       width = 8)

#Depth
pb_depth_pca <- ggplot(pca_scores,aes(x = PC1, y = PC2, color = Depth)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_depth_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_Depth_PCA.pdf",
       height = 8,
       width = 8)

#Brain_ID
pb_brainid_pca <- ggplot(pca_scores,aes(x = PC1, y = PC2, color = Brain_ID)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_brainid_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_BrainID_PCA.pdf",
       height = 8,
       width = 8)

######PCA with selection for highly deviant genes. 
#Perform the PCA analysis with selecting for top 4000 highly deviant genes. 
#Pull HDGs
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:4000]
#Pull log counts matrix
pb_log_counts_HDGs <- pb_log_counts[hdgs,]

#Perform PCA analysis without scaling
pca_no_scale_hdgs <- prcomp(t(pb_log_counts_HDGs), center = TRUE, scale. = FALSE)

#Pull PCA results
pca_scores_hdgs <- as.data.frame(pca_no_scale_hdgs$x)

#Add metadata to the pca output
pca_scores_hdgs$Depth <- colData(sce_pb)$Depth
pca_scores_hdgs$CellType.Final <- colData(sce_pb)$CellType.Final
pca_scores_hdgs$Brain_ID <- colData(sce_pb)$Brain_ID


#Plot the PCA colored by:
#CellType.Final
pb_celltype_pca_hdgs <- ggplot(pca_scores_hdgs,aes(x = PC1, y = PC2, color = CellType.Final)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_celltype_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_HDGs_CellType.Final_PCA.pdf",
       height = 8,
       width = 8)

#Depth
pb_depth_pca <-  ggplot(pca_scores_hdgs,aes(x = PC1, y = PC2, color = Depth)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_depth_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_HDGs_Depth_PCA.pdf",
       height = 8,
       width = 8)

#Brain_ID
pb_brainid_pca <-  ggplot(pca_scores_hdgs,aes(x = PC1, y = PC2, color = Brain_ID)) +
  geom_point(size = 4) +
  labs(x = paste("PC1\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC1"]*100,"% Variance Explained")),
       y = paste("PC2\n",paste0(summary(pca_no_scale_hdgs)$importance[2,"PC2"]*100,"% Variance Explained")))

ggsave(plot = pb_brainid_pca,
       filename = "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/All_Clusters_HDGs_BrainID_PCA.pdf",
       height = 8,
       width = 8)

#Split the object by CellType and then perform the PCA. 
#Use the HDGs from above. 
for(i in unique(sce$CellType.Final)){
  print(i)
  
  #Subset for cell type of interest
  sce_sub <- sce_pb[,sce_pb$CellType.Final == i]
  
  #Pull log counts matrix
  pb_log_counts <- assay(sce_sub,"logcounts")
  pb_log_counts <- pb_log_counts[hdgs,]
  
  #Perform PCA analysis without scaling
  pca_no_scale <- prcomp(t(pb_log_counts), center = TRUE, scale. = FALSE)
  
  #Pull PCA results
  pca_scores <- as.data.frame(pca_no_scale$x)
  
  #Add Depth
  pca_scores$Depth <- colData(sce_sub)$Depth
  
  #Plot
  celltype_plot <- ggplot(pca_scores,aes(x = PC1, y = PC2, color = Depth)) +
    geom_point(size = 4) +
    labs(x = paste("PC1\n",paste0(summary(pca_no_scale)$importance[2,"PC1"]*100,"% Variance Explained")),
         y = paste("PC2\n",paste0(summary(pca_no_scale)$importance[2,"PC2"]*100,"% Variance Explained")))
  
  #Save the plot
  ggsave(plot = celltype_plot,
         filename = paste0("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/",i,"_HDGs_Depth_PCA.pdf"),
         height = 8,
         width = 8)
  
}







