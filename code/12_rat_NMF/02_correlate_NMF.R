# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
# Load libraries
library(SingleCellExperiment)
library(pheatmap)
library(reshape2)
library(Seurat)
library(here)


###### Single nucleus RNA-sequencing object
rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

rm(rat_obj)

nmf_res <- readRDS(here("processed-data","rat_NMF","NMF_Results_k53.Rds")) 

###########Correlate sample with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
data <- as.data.frame(rat_sce$sample)
colnames(data) <- "Sample"
onehot_brain <-  dcast(data = data, rownames(data) ~ Sample, length)
#Make the rownames the first column
rownames(onehot_brain) <- onehot_brain[,1]
#Convert first column to numeric 
onehot_brain[,1]<-as.numeric(onehot_brain[,1])
#Reorder based on first column
onehot_brain <- onehot_brain[order(onehot_brain[,1],decreasing=FALSE),]
#Remove the first column
onehot_brain[,1] <- NULL

#Correlate brain ID with NMF patterns. 
pdf(here("plots","rat_nmf","NMF_Sample_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),onehot_brain), fontsize_row = 5)
dev.off()

###########Correlate QC measures with NMF patterns#########
##onehot encoding of the brain variable##
#Make a dataframe of brain IDs 
QC_vars <- colData(rat_sce)[,c("nCount_RNA","nFeature_RNA","percent.mt")]
#Convert to matrix
QC_vars <- as.matrix(QC_vars)

pdf(here("plots","rat_nmf","NMF_QC_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),QC_vars), fontsize_row = 5)
dev.off()

###########Correlate CellType with NMF patterns#########
#Create datafarme of celltype 
ct_data <- as.data.frame(rat_sce$ident)
colnames(ct_data) <- "CellType" 

#One hot encode the cell type
onehot_CellType <-  dcast(data = ct_data, rownames(ct_data) ~ CellType, length)
rownames(onehot_CellType) <- onehot_CellType[,1]
onehot_CellType[,1] <- as.numeric(onehot_CellType[,1])
onehot_CellType <- onehot_CellType[order(onehot_CellType[,1],decreasing=FALSE),]
onehot_CellType[,1]<-NULL


###correlate with nmf patterns
pdf(here("plots","rat_nmf","NMF_CellType_correlation_heatmap.pdf"))
pheatmap(cor(t(nmf_res@h),onehot_CellType), fontsize_row = 5)
dev.off()


########### Aggregate NMF patterns #########
# create dataframe
aggr_data <- data.frame(colData(rat_sce), t(nmf_res@h))

# aggregate NMF patterns across cell types
# grep "NMF" to get all NMF patterns
aggr_data2 <- aggregate(x = aggr_data[,grep("nmf", colnames(aggr_data))],
                        by = list(aggr_data$ident),
                        FUN = mean)

aggr_data2[1:5,1:5]


# move Group.1 to row names, then drop
rownames(aggr_data2) <- aggr_data2$Group.1
aggr_data2 <- aggr_data2[,-1]

pdf(here("plots","rat_nmf","NMF_CellType_correlation_aggregated_heatmap.pdf"))
pheatmap(aggr_data2,
         color=colorRampPalette(c("blue","white","red"))(100),
         cluster_cols=T,
         cluster_rows=T,
         scale="column",
         fontsize_col = 5)
dev.off()



#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
