#module load r_nac
#Goal: Perform pseudobulking and PCA analysis within a celltype

library(SingleCellExperiment)
library(sessioninfo)
library(ggplot2)
library(ggpubr)
library(scater)
library(scran)
library(scry)
library(here)
library(MASS)

#Pull the iteration 
i <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))

#load the sce object
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

sce

#Remove the neuronal ambiguous population from further analysis. 
sce <- sce[,sce$CellType.Final != "Neuron_Ambig"]

#Do any genes have 0 counts for every cell. 
table(rowSums(assay(sce, "counts")) == 0)

#Remove the genes with 0 counts
sce <- sce[!rowSums(assay(sce, "counts")) == 0, ]

dim(sce)

#For this each sample needs an anterior/middle/posterior designation. 
#Make a dataframe. 
Ant_Mid_Post <- data.frame(Brain_ID = unique(sce$Brain_ID))

Ant_Mid_Post

Ant_Mid_Post <- cbind(Ant_Mid_Post,c("Posterior","Middle","Middle","Anterior","Anterior","Posterior","Anterior","Middle","Middle","Posterior"))

colnames(Ant_Mid_Post)[2] <- "Depth"

Ant_Mid_Post

#Add depth to the sce object
sce$Depth <- Ant_Mid_Post[match(sce$Brain_ID,Ant_Mid_Post$Brain_ID),"Depth"]

table(sce$Depth)

sce$Depth <- factor(x = sce$Depth,levels=c("Anterior","Middle","Posterior"))

class(sce$Depth)

#Get celltype
celltypes <- unique(sce$CellType.Final)
celltype_run <- celltypes[[i]]
message("Running CellType: ", celltype_run)

#Remove the Neuron_Ambig group
sce <- sce[,sce$CellType.Final == celltype_run]

dim(sce)

#compute deviance statistics
# Do for each celltype because overall distribution of gene expression will change bc specific to this celltype
set.seed(1635)
sce <- devianceFeatureSelection(sce,
                                assay = "counts",
                                fam = "binomial",
                                sorted = FALSE,
                                batch = sce$Sample)

#Take top 4000 HDGs
hdgs <- rownames(sce)[order(rowData(sce)$binomial_deviance, decreasing = T)][1:4000]

#Pseudobulk the subsetted object
sce_pb <- aggregateAcrossCells(sce,ids = colData(sce)[,"Brain_ID"])

sce_pb


#Compute logcounts
sce_pb <- computeLibraryFactors(sce_pb)

#Generate log-normalized counts
sce_pb <- logNormCounts(sce_pb)

sce_pb

#Pull log counts matrix
pb_log_counts <- assay(sce_pb,"logcounts")
pb_log_counts <- pb_log_counts[hdgs,]

#Perform PCA analysis without scaling
pca_no_scale <- prcomp(t(pb_log_counts), center = TRUE, scale. = FALSE)
pca_pb <- as.data.frame(pca_no_scale$x[,seq_len(ncol(sce_pb))])


#Add metadata from the pca to the pseudobulked object
#code modified from: https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/08_pseudo_bulking/01_create_pseudobulk_data.R
metadata(sce_pb) <- list("PCA_var_explained" = (summary(pca_no_scale))$importance[2, seq_len(ncol(sce_pb))])
reducedDims(sce_pb) <- list(PCA=pca_pb)

#Add Principal components to the column data
colData(sce_pb) <- cbind(colData(sce_pb),as.data.frame(pca_pb))

celltype_plot <- plotReducedDim(sce_pb,
                                dimred = "PCA",
                                colour_by = "Depth",
                                percentVar = metadata(sce_pb)$PCA_var_explained*100,
                                label_format = c("%s %i", " (%i%%)"),
                                ncomponents = 4,
                                point_size = 3)

#Save the plot
ggsave(plot = celltype_plot,
       filename = paste0("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/PCA_Plots/Cell_Type_Specific/Anterior_Middle_Posterior/",
                         celltype_run,
                         "_HDGs_Depth_3Levels_PCA.pdf"),
       height = 12,
       width = 12)


###########Regression 
x <- colData(sce_pb)

#Force the depth variable to be numeric
x$depth_numeric <- as.numeric(x$Depth)

#Print the depth and depth numeric columns
x[,c("Depth","depth_numeric")]

#Create a directory to put these plots in. 
dir_name <- paste0("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/plots/01_Depth_DE_Testing/Regression/Anterior_Middle_Posterior/",celltype_run)
dir.create(dir_name)

#Make a dataframe to store the R-squared values and the p-values from the regression. 
stat_df <- as.data.frame(matrix(nrow=10,ncol=2))
colnames(stat_df) <- c("R-squared","p")
rownames(stat_df) <- paste0("PC",1:10)

#Add % variance explained to the dataframe
stat_df$Prop_var_explained <- summary(pca_no_scale)$importance["Proportion of Variance", seq_len(ncol(sce_pb))]

for(l in paste0("PC",1:10)){
  print(l)
  #Regression 
  set.seed(1101)
  model_linear <- lm(get(l) ~ depth_numeric, data = x)
  #Pull stats
  p_val <- round(summary(model_linear)$coefficients[2,4],4)
  r_sq  <- round(summary(model_linear)$adj.r.squared,4)
  #Add stats to the dataframe
  stat_df[l,"R-squared"] <- r_sq
  stat_df[l,"p"] <- p_val
  #Make plot
  p <- ggplot(x,aes(y = .data[[l]], x = depth_numeric, color = Depth)) +
    geom_point(size = 4,alpha = 0.5) +
    geom_smooth(method = "lm",color = "black",fill = "lightgrey",se = TRUE) +
    scale_x_continuous(breaks = 1:2, labels = c("anterior_posterior", "middle")) +
    labs(x = "Anatomical Depth",
         y = l,
         title = sprintf("R^2 = %s, p = %s", r_sq, p_val)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  #Save plot
  ggsave(p,filename = paste0(dir_name,"/",l,".pdf"))
}

#Add a column that is the principal component and cell type. 
stat_df$CellType <- celltype_run
stat_df$PCA <- rownames(stat_df)

#Reorder the dataframe
stat_df <- stat_df[,c("PCA","CellType","Prop_var_explained","p","R-squared")]

#Now save the dataframe
saveRDS(stat_df,file = here("processed-data","Regression_stats","Anterior_Middle_Posterior",paste0(celltype_run,"_AMP_Stats.Rds")))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
