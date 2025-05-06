#Goal: Additional exploratory analysis of singleR results
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.4.x
library(SpatialExperiment)
library(ComplexHeatmap)
library(sessioninfo)
library(ggplot2)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe
# class: SpatialExperiment 
# dim: 541 2375058 
# metadata(11): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
# DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2375058): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(56): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe
# class: SpatialExperiment 
# dim: 366 2375058 
# metadata(11): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(2375058): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(56): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

#Load all of the SingleR results
filenames <- list.files(here("processed-data","06_label_transfer"),full.names = TRUE)
file_list <- lapply(filenames, function(x) get(load(x)))

res <- do.call(what = rbind,file_list)
dim(res)
#[1] 2375058       4

#Create a key column 
res$key <- rownames(res)
colData(spe)$key <- rownames(colData(spe))

#Put res in the order of the spe object
res <- res[match(colData(spe)$key,res$key),]

#Some sanity checks that it worked. 
identical(res$key,spe$key)
#[1] TRUE

identical(res$key,rownames(colData(spe)))
#[1] TRUE

#Add the labels and pruned labels to the object
colData(spe) <- cbind(colData(spe),res[,c("labels","pruned.labels")])

#Calculate percentage of all cells in each sample that are a certain cell type. 
celltype_by_sample <- as.data.frame.matrix(table(spe$pruned.labels,spe$Sample))
colnames(celltype_by_sample)
# [1] "Br6660_NAc1_580"   "Br6660_Nac10_4080" "Br6660_Nac11_5580"
# [4] "Br6660_NAc2_1090"  "Br6660_NAc3_1580"  "Br6660_NAc4_2080" 
# [7] "Br6660_NAc5_2580"  "Br6660_NAc6_3080"  "Br6660_NAc7_3580" 
# [10] "Br6660_NAc8_4580"  "Br6660_NAc9_5080" 

celltype_by_sample <- celltype_by_sample[,c(1,4:11,2,3)]
colnames(celltype_by_sample)
# [1] "Br6660_NAc1_580"   "Br6660_NAc2_1090"  "Br6660_NAc3_1580" 
# [4] "Br6660_NAc4_2080"  "Br6660_NAc5_2580"  "Br6660_NAc6_3080" 
# [7] "Br6660_NAc7_3580"  "Br6660_NAc8_4580"  "Br6660_NAc9_5080" 
# [10] "Br6660_Nac10_4080" "Br6660_Nac11_5580"

#Make a ComplexHeatmap that contains all of the CellTypes 
celltype_by_sample_pct

#Calculate percentages 
celltype_by_sample_pct <- sweep(celltype_by_sample,MARGIN = 2,colSums(celltype_by_sample),"/") * 100
celltype_by_sample_pct$CellType <- rownames(celltype_by_sample_pct)


#Make a heatmap of celltype_by_sample_pct with ComplexHeatmap
ct_mat <- celltype_by_sample_pct[,-12]
ct_mat <- as.matrix(ct_mat)

#Create a row annotation that will allow for consistent colors between snRNA-seq and this. 

# Use rownames of matrix as class labels
class <- factor(rownames(ct_mat))

# Subset color vector to just those classes
cluster_cols_list <- list(class = cluster_cols[levels(class)])

# Row annotation
row_ha <- rowAnnotation(class = class, col = cluster_cols_list)

min(ct_mat)
#[1] 0.03725362

max(ct_mat)
#[1] 42.25915


# Draw heatmap
hm1 <- Heatmap(ct_mat, 
               name = "Percent of Sample",
               right_annotation = row_ha,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "black", lwd = 1))

pdf(file = here("plots","06_label_transfer","SingleR_Predictions_all.pdf"),
    width = 8,
    height = 8)
draw(hm1)
dev.off()


#Melt dataframe and plot
celltype_by_sample_pct_melt <- reshape2::melt(celltype_by_sample_pct)

#Load cluster colors
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE)
# Loading objects:
#   cluster_cols

#Remove Neuron_Ambig
cluster_cols <- cluster_cols[-14]

p <- ggplot(data = celltype_by_sample_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(stat="identity",position = "stack") +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  labs(x = "Sample", y = "% of Sample",fill = "Predicted\nCell Type") 
ggsave(plot = p,filename = here("plots","06_label_transfer","Pct_Predicted_CellType_Per_Sample.pdf"))

#Calculate the percentage of D1 cells split by population
D1s <- celltype_by_sample[grep("DRD1",rownames(celltype_by_sample)),]
D1s_pct <- sweep(D1s,MARGIN = 2,colSums(D1s),"/") * 100
D1s_pct$CellType <- rownames(D1s_pct)
D1s_pct_melt <- reshape2::melt(D1s_pct)

p2 <- ggplot(data = D1s_pct_melt,aes(x = variable,y = value,fill = CellType)) +
  geom_bar(stat="identity",position = "stack") +
  scale_fill_manual(values = cluster_cols) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
  labs(x = "Sample", y = "% of DRD1-MSNs",fill = "Predicted\nCell Type") 
ggsave(plot = p2,filename = here("plots","06_label_transfer","D1s_Predicted_Pct_Per_Sample.pdf"))


for(i in unique(D1s_pct_melt$CellType)){
  x <- subset(D1s_pct_melt,subset=(CellType == i))
  p3 <- ggplot(data = x,aes(x = variable,y = value,fill = CellType)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    labs(x = "Sample", y = "% of DRD1-MSNs",fill = "Predicted\nCell Type") 
  ggsave(plot = p3,filename = here("plots","06_label_transfer",paste0(i,"_Predicted_PctofD1s_Per_Sample.pdf")))
}


#Oligos dominate the heatmap above to Per to the point you can't really tell any differences. 
#What if we did the percentage of all found 
celltype_pct <- sweep(celltype_by_sample,MARGIN = 1,rowSums(celltype_by_sample),"/") * 100
celltype_pct <- as.matrix(celltype_pct)

range(celltype_pct)
#[1]  0.2720976 28.4759719


class <- factor(rownames(celltype_pct))
cols_list <- list(class = cluster_cols[levels(class)])

col_fun <- circlize::colorRamp2(c(0,10,20),c("#2166AC", "white", "#B2182B"))

# Row annotation
row_ha_neurons <- rowAnnotation(class = class, col = cols_list)
hm2 <- Heatmap(celltype_pct, 
               name = "Percent of Cell Type",
               right_annotation = row_ha_neurons,
               cluster_columns = FALSE,
               rect_gp = gpar(col = "black", lwd = 1),
               col = col_fun)

pdf(file = here("plots","06_label_transfer","SingleR_Predictions_all_PercentofcellType.pdf"),
    width = 8,
    height = 8)
draw(hm2)
dev.off()

#Like the D1 plots above,make a bargraph for each cell type
for(i in unique(celltype_by_sample_pct_melt$CellType)){
  x <- subset(celltype_by_sample_pct_melt,subset=(CellType == i))
  p <- ggplot(data = x,aes(x = variable,y = value,fill = CellType)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    labs(x = "Sample", y = "% of Sample",fill = "Predicted\nCell Type")
  ggsave(plot = p,filename = here("plots","06_label_transfer","Percent_Sample",paste0(i,"_Predicted_Pct_Per_Sample.pdf")))
}


#Make sample plots but as percent of each cell type represented 
celltype_pct_melt <- reshape2::melt(celltype_pct)
for(i in unique(celltype_pct_melt$Var1)){
  x <- subset(celltype_pct_melt,subset=(Var1 == i))
  p <- ggplot(data = x,aes(x = Var2,y = value,fill = Var1)) +
    geom_bar(stat="identity") +
    scale_fill_manual(values = cluster_cols) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45,hjust = 1)) +
    labs(x = "Sample", y = "% of Cell Type",fill = "Predicted\nCell Type")
  ggsave(plot = p,filename = here("plots","06_label_transfer","Percent_CellType",paste0(i,"_Predicted_Pct_Per_CellType.pdf")))
}


###Reproduciblity
sessionInfo()
# R version 4.4.3 Patched (2025-02-28 r88189)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.4.x/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.0
# 
# locale:
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] grid      stats4    stats     graphics  grDevices datasets  utils     methods   base     
# 
# other attached packages:
#   [1] ComplexHeatmap_2.22.0       here_1.0.1                  ggplot2_3.5.2               sessioninfo_1.2.3          
# [5] SpatialExperiment_1.16.0    SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0 Biobase_2.66.0             
# [9] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3         IRanges_2.40.1              S4Vectors_0.44.0           
# [13] BiocGenerics_0.52.0         MatrixGenerics_1.18.1       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] shape_1.4.6.1           circlize_0.4.16         gtable_0.3.6            rjson_0.2.23           
# [5] GlobalOptions_0.1.2     lattice_0.22-6          Cairo_1.6-2             vctrs_0.6.5            
# [9] tools_4.4.3             generics_0.1.3          parallel_4.4.3          tibble_3.2.1           
# [13] cluster_2.1.8           pkgconfig_2.0.3         Matrix_1.7-2            RColorBrewer_1.1-3     
# [17] lifecycle_1.0.4         GenomeInfoDbData_1.2.13 compiler_4.4.3          farver_2.1.2           
# [21] stringr_1.5.1           textshaping_1.0.1       codetools_0.2-20        clue_0.3-66            
# [25] pillar_1.10.2           crayon_1.5.3            DelayedArray_0.32.0     magick_2.8.6           
# [29] iterators_1.0.14        abind_1.4-8             foreach_1.5.2           tidyselect_1.2.1       
# [33] digest_0.6.37           stringi_1.8.7           dplyr_1.1.4             reshape2_1.4.4         
# [37] labeling_0.4.3          rprojroot_2.0.4         colorspace_2.1-1        cli_3.6.5              
# [41] SparseArray_1.6.2       magrittr_2.0.3          S4Arrays_1.6.0          dichromat_2.0-0.1      
# [45] withr_3.0.2             scales_1.4.0            UCSC.utils_1.2.0        XVector_0.46.0         
# [49] httr_1.4.7              ragg_1.4.0              png_0.1-8               GetoptLong_1.0.5       
# [53] doParallel_1.0.17       rlang_1.1.6             Rcpp_1.0.14             glue_1.8.0             
# [57] jsonlite_2.0.0          R6_2.6.1                plyr_1.8.9              systemfonts_1.2.3      
# [61] zlibbioc_1.52.0        

