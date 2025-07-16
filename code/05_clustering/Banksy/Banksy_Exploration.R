#Goal: Load Banksy results and explore
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.4.x
library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(escheR)
library(here)

#Load spe object containing normalized coutns
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))

spe
# class: SpatialExperiment 
# dim: 541 4884175 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
#   DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4884175): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
#   Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Subset for only probes
gene_expression_idx <- which(rowData(spe)$Type == "Gene Expression")
spe <- spe[gene_expression_idx,]

spe
# class: SpatialExperiment 
# dim: 366 4884175 
# metadata(22): Samples Samples ... Samples Samples
# assays(3): counts nucleus_normcounts cell_normcounts
# rownames(366): ABCC9 ADAMTS12 ... ZBBX ZDHHC23
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4884175): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
#   Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
# mainExpName: NULL
# altExpNames(0):
# spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Add Banksy clusters to the object. 
Banksy_clusters <- read.csv(here("processed-data","05_Clustering","Banksy_clusters.csv"))
head(Banksy_clusters)

table(Banksy_clusters$V1)
      1       2       3       4       5       6       7       8       9      10 
1245255  645179  624322  616903  345006  327998  295049  219861  192707  115127 
     11      12      13      14      15 
 114991   60695   56504   24434     144 

#X is just rownames, remove it. 
Banksy_clusters <- Banksy_clusters[,-1]

#Rename the columns
colnames(Banksy_clusters) <- c("Banksy_Cluster","key")

identical(rownames(colData(spe)),Banksy_clusters$key)
#[1] TRUE

#Just cbind the clusters into the colData
colData(spe)$Banksy_Cluster <- as.factor(Banksy_clusters$Banksy_Cluster)

###########Complext heatmap of basic markers to classify Banksy_Cluster
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$Banksy_Cluster)

############set up columns for heatmaps. 
#Set marker genes to be included on the heatmap.  
Panel <- readxl::read_excel(here("processed-data","Files_For_Upload","NAc_Xenium_Panel_Final_withNotes.xlsx"))
markers_all <- Panel$Gene[order(Panel$Cell_Type)]

marker_labels <-  Panel$Cell_Type[order(Panel$Cell_Type)]
marker_labels <- factor(x = marker_labels,
                        levels =  unique(marker_labels))

colors_markers <- list(marker = c(`D1-MSN` = "#332288",
                                  `D2-MSN` = "#117733",
                                  `D1-Islands` = "#44AA99",
                                  MSNs = "#88CCEE",
                                  `Inhibitory-Neuron` = "#DDCC77",
                                  `Spatial-Marker` = "#8C564B",
                                  Astrocytes = "#D81B60",
                                  Endothelial = "#CC6677",
                                  Ependymal = "#AA4499",
                                  Microglia = "#882255",
                                  Other = "#FE6100",
                                  Excitatory = "#FFB000",
                                  Neurons = "#1E88E5",
                                  OT = "gray",
                                  `Lateral-Septum` = "black",
                                  `Diagonal-Band` = "white",
                                  `OUD-GWAS` = "green",
                                  `SCZ-GWAS` = "blue"))

correct_order <- unique(Panel$Cell_Type[order(Panel$Cell_Type)])

colors_markers$marker <- colors_markers$marker[correct_order]

col_ha <- ComplexHeatmap::columnAnnotation(marker = marker_labels,
                                           show_annotation_name = FALSE,
                                           show_legend = TRUE,
                                           col = colors_markers)

###########set up rows for heatmap. 
dat <- assay(spe,"nucleus_normcounts")
#rownames(dat) <- rowData(sce)$gene_name
dim(dat)

dat <- dat[markers_all,]
dim(dat)
dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)
max(hm_mat)
#[1] 3.321265

min(hm_mat)
#[1] -1.944437

col_fun <- circlize::colorRamp2(c(-2,0,3.5),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              bottom_annotation = col_ha,
                              row_title = NULL,
                              col = col_fun)


pdf(file = here("plots","05_clustering","Banksy_Clusters_with100genepanel_columns_clustered.pdf"),
    width = 18,
    height = 12)
ComplexHeatmap::draw(hm)
dev.off()

hm2 <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = TRUE,
                              cluster_columns = FALSE,
                              bottom_annotation = col_ha,
                              row_title = NULL,
                              col = col_fun)


pdf(file = here("plots","05_clustering","Banksy_Clusters_with100genepanel.pdf"),
    width = 18,
    height = 12)
ComplexHeatmap::draw(hm2)
dev.off()


#Name the clusters
cluster_annotation_df <- data.frame(Banksy_cluster = levels(spe$Banksy_Cluster),
                                    Annotation = c("Astrocyte-A","WM","Astrocyte-B","MSN-1","Excitatory","Endothelial",
                                                   "Microglia","D1_Island_A","Inh_A","Other","Ependymal",
                                                   "D1_Island_B","Chat_Inh"))

spe$Annotation <- cluster_annotation_df$Annotation[match(spe$Banksy_Cluster,
                                                         cluster_annotation_df$Banksy_cluster)]


#Read in the singleR results. Maybe that will provide some information as to what "Other is"
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
#Just use the labels column 
colData(spe) <- cbind(colData(spe),res[,c("labels","pruned.labels")])

labels_annotation <- as.data.frame.matrix(table(spe$pruned.labels,spe$Annotation))

#Calculate percentages 
labels_annotation_pct <- sweep(labels_annotation,MARGIN = 2,colSums(labels_annotation),"/") * 100
labels_annotation_pct$snRNAseq_cluster <- rownames(labels_annotation_pct)
labels_annotation_pct_melt <- reshape2::melt(labels_annotation_pct)
#Using snRNAseq_cluster as id variables

load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",
     verbose = TRUE)
# Loading objects:
#   cluster_cols

cluster_cols
# Oligo  DRD1_MSN_A  DRD2_MSN_A         OPC   Microglia   Ependymal 
# "#4F4753"   "#ECA31C"   "#58B6ED"   "#0D9F72"   "#F2E642"   "#0077B9" 
# Astrocyte_A  DRD1_MSN_B Endothelial       Inh_A  DRD2_MSN_B Astrocyte_B 
# "#D95F00"   "#D079AA"   "#D00DFF"   "#35FB00"   "#F80091"   "#FF0016" 
# DRD1_MSN_C Neuro_Ambig  DRD1_MSN_D       Inh_B       Inh_C       Inh_D 
# "#2A4BF9"      "grey"   "#FB3DD9"   "#7A0096"   "#854222"   "#A7F281" 
# Inh_E  Excitatory       Inh_F 
# "#0DFBFA"   "#5C6300"     "black" 

cluster_cols <- cluster_cols[-14]


p1 <- ggplot(data = labels_annotation_pct_melt,aes(x = variable,y=value,fill =snRNAseq_cluster)) +
  geom_bar(position = "stack",stat = "identity") +
  scale_fill_manual(values = cluster_cols) +
  labs(x = "Annotation",
       y = "% of Annotation",
       fill = "snRNA-seq\nCluster") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust = 1))
ggsave(filename = here("plots","05_clustering","Annotation_percentage_snRNA-seq_cluster_bar.pdf"),plot = p1)


#Change other to Inh_B
spe$Annotation <- ifelse(spe$Annotation == "Other",
                         "Inh_B",
                         spe$Annotation)

#Create a color palette that makes sense. 
#make some new brain colors
cluster_cols <- Polychrome::createPalette(length(unique(spe$Annotation)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(cluster_cols) <- unique(spe$Annotation)

#Switch colors for Excitatory and D1_Island_B
cluster_cols[1] <- "#B2E0FF" 
cluster_cols[13] <- "#D70D5D"

#Plot with escheR 
for (i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[, spe$Sample == i]
  
  p <-make_escheR(sub_spe) |>
    add_fill("Annotation")+
    ggtitle(i) +
    scale_fill_manual(values = cluster_cols) +
    theme(element_text(hjust = 0.5))

  png(filename = here("plots","05_clustering","Banksy",paste0(i,".png")),
      units = "in",height = 20,width = 20,res = 300)
  print(p)
  dev.off()
  
}

#Also remake the heatmap with domains named. 
###########Complext heatmap of basic markers to classify Banksy_Cluster
#Code from https://github.com/LieberInstitute/septum_lateral/blob/main/snRNAseq_mouse/code/02_analyses/Complex%20Heatmap.R
splitit <- function(x) split(seq(along = x), x)

cell_idx <- splitit(spe$Annotation)

##column annotation can be used from above. Now set up the rows. 
###########set up rows for heatmap. 
population <- unique(spe$Annotation)

pop_markers <- list(population = c(Excitatory = "#B2E0FF" ,
                                   WM = "#0084E3",
                                   `Astrocyte-B` = "#FFCD0D",
                                   `Astrocyte-A` = "#F5700D",
                                   `MSN-1` =  "#1C5342",
                                   Microglia = "#00FB16",
                                   Inh_A = "#FD00FB",
                                   Endothelial = "#FFB0E4",
                                   Inh_B = "#00FBC8",
                                   Chat_Inh = "#834522",
                                   Ependymal =  "#86B235",
                                   D1_Island_A =   "#9916FF",
                                   D1_Island_B = "#D70D5D"))

row_ha <- ComplexHeatmap::rowAnnotation(population = population,
                                        show_annotation_name = FALSE,
                                        show_legend = TRUE,
                                        col = pop_markers)

###########set up rows for heatmap. 
dat <- assay(spe,"nucleus_normcounts")
#rownames(dat) <- rowData(sce)$gene_name
dim(dat)

dat <- dat[markers_all,]
dim(dat)
dat <- as.matrix(dat)

hm_mat <- scale(t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i])))),
                center = TRUE,
                scale = TRUE)
max(hm_mat)
#[1] 3.321265

min(hm_mat)
#[1] -1.944437
#hm_mat <- t(do.call(cbind, lapply(cell_idx, function(i) rowMeans(dat[markers_all, i]))))

hm_mat <- hm_mat[population,]

col_fun <- circlize::colorRamp2(c(-2,0,3.5),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              bottom_annotation = col_ha,
                              right_annotation = row_ha,
                              row_title = NULL,
                              col = col_fun)


pdf(file = here("plots","05_clustering","Annotated_Banksy_Clusters_with100genepanel_columns_clustered.pdf"),
    width = 18,
    height = 12)
ComplexHeatmap::draw(hm)
dev.off()

hm2 <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                               name = "centered,scaled",
                               cluster_rows = TRUE,
                               cluster_columns = FALSE,
                               bottom_annotation = col_ha,
                               row_title = NULL,
                               col = col_fun)


pdf(file = here("plots","05_clustering","Annotated_Banksy_Clusters_with100genepanel.pdf"),
    width = 18,
    height = 12)
ComplexHeatmap::draw(hm2)
dev.off()


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
#   [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
# [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
# [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
# [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
# [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: US/Eastern
# tzcode source: system (glibc)
# 
# attached base packages:
#   [1] stats4    stats     graphics  grDevices datasets  utils     methods  
# [8] base     
# 
# other attached packages:
#   [1] here_1.0.1                  escheR_1.6.0               
# [3] ggplot2_3.5.2               sessioninfo_1.2.3          
# [5] SpatialExperiment_1.16.0    SingleCellExperiment_1.28.1
# [7] SummarizedExperiment_1.36.0 Biobase_2.66.0             
# [9] GenomicRanges_1.58.0        GenomeInfoDb_1.42.3        
# [11] IRanges_2.40.1              S4Vectors_0.44.0           
# [13] BiocGenerics_0.52.0         MatrixGenerics_1.18.1      
# [15] matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] gtable_0.3.6            circlize_0.4.16         shape_1.4.6.1          
# [4] rjson_0.2.23            GlobalOptions_0.1.2     lattice_0.22-6         
# [7] Cairo_1.6-2             vctrs_0.6.5             tools_4.4.3            
# [10] generics_0.1.3          parallel_4.4.3          Polychrome_1.5.1       
# [13] tibble_3.2.1            cluster_2.1.8           pkgconfig_2.0.3        
# [16] Matrix_1.7-2            RColorBrewer_1.1-3      scatterplot3d_0.3-44   
# [19] readxl_1.4.5            lifecycle_1.0.4         GenomeInfoDbData_1.2.13
# [22] stringr_1.5.1           compiler_4.4.3          farver_2.1.2           
# [25] textshaping_1.0.1       codetools_0.2-20        ComplexHeatmap_2.22.0  
# [28] clue_0.3-66             pillar_1.10.2           crayon_1.5.3           
# [31] DelayedArray_0.32.0     magick_2.8.6            iterators_1.0.14       
# [34] abind_1.4-8             foreach_1.5.2           tidyselect_1.2.1       
# [37] digest_0.6.37           stringi_1.8.7           reshape2_1.4.4         
# [40] dplyr_1.1.4             labeling_0.4.3          rprojroot_2.0.4        
# [43] grid_4.4.3              colorspace_2.1-1        cli_3.6.5              
# [46] SparseArray_1.6.2       magrittr_2.0.3          S4Arrays_1.6.0         
# [49] dichromat_2.0-0.1       withr_3.0.2             scales_1.4.0           
# [52] UCSC.utils_1.2.0        XVector_0.46.0          httr_1.4.7             
# [55] cellranger_1.1.0        ragg_1.4.0              png_0.1-8              
# [58] GetoptLong_1.0.5        doParallel_1.0.17       rlang_1.1.6            
# [61] Rcpp_1.0.14             glue_1.8.0              jsonlite_2.0.0         
# [64] plyr_1.8.9              R6_2.6.1                systemfonts_1.2.3      
# [67] zlibbioc_1.52.0   
