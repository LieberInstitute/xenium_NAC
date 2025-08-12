#Goal: Load Banksy results and explore
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
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
# DeprecatedCodeword_0393
# rowData names(8): ID Symbol ... subsets_any_neg subsets_GEX
# colnames(4884175): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
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
# Br6436_Nac_11_5650_68180 Br6436_Nac_11_5650_68181
# colData names(59): Sample Barcode ... cell_area.sf nucleus_area.sf
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_final y_final
# imgData names(1): sample_id

#Add Banksy clusters to the object. 
Banksy_clusters <- read.csv(here("processed-data","05_Clustering","Banksy_clusters.csv"))
head(Banksy_clusters)
# X V1                V2
# 1 1  6 Br6660_NAc1_580_1
# 2 2  6 Br6660_NAc1_580_2
# 3 3  6 Br6660_NAc1_580_3
# 4 4  6 Br6660_NAc1_580_4
# 5 5  6 Br6660_NAc1_580_5
# 6 6  6 Br6660_NAc1_580_6

table(Banksy_clusters$V1)
# 1       2       3       4       5       6       7       8       9      10 
# 1245255  645179  624322  616903  345006  327998  295049  219861  192707  115127 
# 11      12      13      14      15 
# 114991   60695   56504   24434     144 

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
#[1] 3.605704

min(hm_mat)
#-2.331799

col_fun <- circlize::colorRamp2(c(-2.5,0,3.75),c("blue","white","red"))

hm <- ComplexHeatmap::Heatmap(matrix = hm_mat,
                              name = "centered,scaled",
                              cluster_rows = TRUE,
                              cluster_columns = TRUE,
                              bottom_annotation = col_ha,
                              row_title = NULL,
                              col = col_fun)


pdf(file = here("plots","05_clustering","Banksy","Banksy_Clusters_with100genepanel_columns_clustered.pdf"),
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


pdf(file = here("plots","05_clustering","Banksy","Banksy_Clusters_with100genepanel.pdf"),
    width = 18,
    height = 12)
ComplexHeatmap::draw(hm2)
dev.off()


#Create a color palette that makes sense. 
#Plot Banksy clusters on the tissue sections to help with annotation. 
cluster_cols <- Polychrome::createPalette(length(unique(spe$Banksy_Cluster)),
                                          c("#D81B60", "#1E88E5","#ffcd14","#ff7814","#004D40"))
names(cluster_cols) <- unique(spe$Banksy_Cluster)

#Plot with escheR 
for (i in unique(spe$Sample)){
  print(i)
  sub_spe <- spe[, spe$Sample == i]
  
  p <-make_escheR(sub_spe) |>
    add_fill("Banksy_Cluster")+
    ggtitle(i) +
    scale_fill_manual(values = cluster_cols) +
    theme(element_text(hjust = 0.5))
  
  png(filename = here("plots","05_clustering","Banksy","Pre_Annotation",paste0(i,"_PreAnnotation.png")),
      units = "in",height = 20,width = 20,res = 300)
  print(p)
  dev.off()
  
}

table(spe$Sample,spe$Banksy_Cluster)
# 1     2     3     4     5     6     7     8     9
# Br6660_NAc1_580    45864 14909 26384 42426 12775  8271  8282  9316  1365
# Br6660_NAc2_1090   29867 12250 17411 25588  8433 11212  5416  5977  1588
# Br6660_NAc3_1580   50536 25039 28480 34854 13301 34183 10192  8439  4872
# Br6660_NAc4_2080   50709 22763 27865 28181 14618 24549  4597  8976 11898
# Br6660_NAc5_2580   70382 19535 38875 30332 16856 27914 10138 11169 10436
# Br6660_NAc6_3080   67554 22573 36065 24392 17490 16999 16255 10228  8846
# Br6660_NAc7_3580   64676 22934 37024 20519 17985 23220 15605  9844  9097
# Br6660_NAc8_4580   53188 46764 28659 21805 20644 23645  6789 10553 15760
# Br6660_NAc9_5080   53336 50278 24604 26900 18944 21285 10412 11134 20054
# Br6660_Nac10_4080  62849 42419 36544 22515 19489 15116  9185 11335 14701
# Br6660_Nac11_5580  30674 45837 12866 28113 17372 22092  4274  8451 13400
# Br6436_Nac1_650    82827 27056 40964 41614 17655 16431 25988 15295  9226
# Br6436_Nac2_1150   71182 21827 33908 31246 15036  7195 16656 12313  9328
# Br6436_Nac3_1650   72166 27645 35330 32406 17837 10386 22866 13279 10487
# Br6436_Nac_4_2150  64678 31690 29967 29388 14137 12129 19473 10225  6774
# Br6436_Nac_5_2650  56653 28667 26172 27150 15617 13066 15446  8458  5968
# Br6436_Nac_6_3150  84005 40680 33750 40739 20812  9421    21  9842  7650
# Br6436_Nac_7_3660  69738 41488 33910 37178 19041 14621 27114 13503  7884
# Br6436_Nac_8_4150  62574 32466 29865 25558 15519  6881 23586 11260  7543
# Br6436_Nac_9_4650  36991 21727 17425 12246  9850  3751 16105  6847  6025
# Br6436_Nac_10_5150 50565 32321 22285 25598 16265  4090 20299  9879  6657
# Br6436_Nac_11_5650 14241 14311  5969  8155  5330  1541  6350  3538  3148
# 
# 10    11    12    13    14    15
# Br6660_NAc1_580     4944  2905     4  6033    31     0
# Br6660_NAc2_1090    3182  1695    14  1896    21     0
# Br6660_NAc3_1580    4995  3275    22  3577   231     0
# Br6660_NAc4_2080    5573  3231   657   370    16     0
# Br6660_NAc5_2580    7293  4772  1474  2880    22     0
# Br6660_NAc6_3080    6990  7566  3534   778    60     0
# Br6660_NAc7_3580    6542  8642  3083   726    33     0
# Br6660_NAc8_4580    5776  8936  6442  2354    39     0
# Br6660_NAc9_5080    6055  9796  2520  3213    50     0
# Br6660_Nac10_4080   6655 11009  9368  7138    23     0
# Br6660_Nac11_5580   4304  8447  3599   928    65   144
# Br6436_Nac1_650     6284  4137  3244  7374    86     0
# Br6436_Nac2_1150    5349  3933   746 12321    87     0
# Br6436_Nac3_1650    5639  4472  2149  5427   141     0
# Br6436_Nac_4_2150   4932  3793  1577   138   574     0
# Br6436_Nac_5_2650   4244  3708  2851   366  2309     0
# Br6436_Nac_6_3150   5867  4349  3528   229 19335     0
# Br6436_Nac_7_3660   5942  5079  1859   348   426     0
# Br6436_Nac_8_4150   4680  4603  1263    95   168     0
# Br6436_Nac_9_4650   3350  4104  6151    17   112     0
# Br6436_Nac_10_5150  4908  4509  5204   221   453     0
# Br6436_Nac_11_5650  1623  2030  1406    75   152     0

#Use AggregateAcrossCells to Pseudobulk the data and generate some boxplots. 
ids_df <- S4Vectors::DataFrame(
  Sample         = colData(spe)$Sample,
  Banksy_Cluster = colData(spe)$Banksy_Cluster
)

spe_pb <- scuttle::aggregateAcrossCells(spe, ids = ids_df, use.assay.type = "counts")

##AggregateAcrossCells removes the column names. However, it should be the Sample name and banksy cluster
#which is in the column data
#Let's check if that is correct
nrow(colData(spe_pb)) == ncol(spe_pb)
#[1] TRUE

#Spot check counts
# Pick a pseudobulk column
test_col <- 1
test_sample  <- colData(sce_pb)$Sample[test_col]
test_cluster <- colData(sce_pb)$Banksy_Cluster[test_col]

# Cells from original spe that belong to this combination
sel <- colData(spe)$Sample == test_sample &
  colData(spe)$Banksy_Cluster == test_cluster

# Sum their counts
orig_sum <- Matrix::rowSums(assay(spe, "counts")[, sel, drop=FALSE])

# Compare to pseudobulk column
all(orig_sum == assay(sce_pb, "counts")[, test_col])


#Pull counts matrix from spe
spe_counts <- as.matrix(assay(spe,"counts"))

output <- vector()
for(i in 1:ncol(spe_pb)){
  print(i)
  # Pick a pseudobulk column
  test_col <- i
  test_sample  <- colData(spe_pb)$Sample[test_col]
  test_cluster <- colData(spe_pb)$Banksy_Cluster[test_col]
  
  # Cells from original spe that belong to this combination
  cells <- which(colData(spe)$Sample == test_sample & colData(spe)$Banksy_Cluster == test_cluster)
  
  # Sum their counts
  orig_sum <- Matrix::rowSums(assay(spe, "counts")[, cells])
  
  # Compare to pseudobulk column
 output[i] <- identical(orig_sum,assay(spe_pb, "counts")[, test_col])
  
}

table(output)
# output
# TRUE 
#  309 

#Add the column names back. 
colnames(spe_pb) <- paste(colData(spe_pb)$Banksy_Cluster,colData(spe_pb)$Sample,sep = ".")
rownames(colData(spe_pb)) <- paste(colData(spe_pb)$Banksy_Cluster,colData(spe_pb)$Sample,sep = ".")
identical(colnames(spe_pb),rownames(colData(spe_pb)))
#[1] TRUE

# log-normalize
spe_pb <- scuttle::computeLibraryFactors(spe_pb)
spe_pb <- scuttle::logNormCounts(spe_pb)

#Remove the hyphen from a couple of the gene names. 
rownames(spe_pb) <- gsub("-","",rownames(spe_pb))


for(i in rownames(spe_pb)){
  print(i)
  colData(spe_pb)[[i]] <- assay(spe_pb,"logcounts")[i,]
}

#make boxplots to help with annotation 
for(i in rownames(spe_pb)){
  print(i)
  p <- ggplot(as.data.frame(colData(spe_pb)),aes(x = Banksy_Cluster,y = .data[[i]],fill = Banksy_Cluster)) +
    geom_boxplot(outlier.shape = NA) + #Don't plot outlier twice
    geom_jitter(alpha = 0.6) +
    ggtitle(i) +
    theme_bw() +
    theme(legend.position = "none")
  ggsave(p,filename = here("plots","05_clustering","Banksy","Pre_Annotation",
                           "Pseudobulked_Boxplots",paste0(i,".pdf")),
         height = 5, width = 6.5)
}

#Also ran SingleR which is an unbiased way to use snRNA-seq for annotation of single cell data. 
#Load all of the SingleR results
filenames <- list.files(here("processed-data","06_label_transfer"),full.names = TRUE)
file_list <- lapply(filenames, function(x) get(load(x)))

res <- do.call(what = rbind,file_list)
dim(res)
#[1] 4774688       4

#Create a key column 
res$key <- rownames(res)
colData(spe)$key <- rownames(colData(spe))

#Put res in the order of the spe object
res <- res[match(colData(spe)$key,res$key),]

#Some sanity checks that it worked. 
identical(res$key,colData(spe)$key)
#[1] TRUE

all(res$key == rownames(colData(spe)))
#[1] TRUE

#Add the labels and pruned labels to the object
colData(spe) <- cbind(colData(spe),res[,c("labels","pruned.labels")])

# Make a stacked bargraph where x is the banksy cluster and y is the percent of taht banksy cluster 
# mapping to each snRNA-seq cluster
df <- as.data.frame(colData(spe))

x <- as.matrix(table(df$labels,df$Banksy_Cluster))
x_pct <- sweep(x = x,MARGIN = 2,STATS = colSums(x),"/")*100
x_pct <- as.data.frame(x_pct)
colnames(x_pct) <- c("snRNA_CellType","Banksy_Cluster","Percent_Banksy_Cluster")

#Load snRNA-seq cell type colors
load("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/070924_21colors_celltypeFinal.rda",verbose = TRUE)
# Loading objects:
#   cluster_cols

cluster_cols <- cluster_cols[-14]

p <- ggplot(x_pct,aes(x = Banksy_Cluster,y = Percent_Banksy_Cluster,fill = snRNA_CellType)) +
  geom_bar(stat = "identity") + 
  scale_fill_manual(values = cluster_cols)

ggsave(p,filename = here("plots","05_clustering","Banksy","Pre_Annotation","SingleR_PreAnnotation_StackedBar.pdf"))

sessionInfo()
# R version 4.5.0 Patched (2025-05-21 r88220)
# Platform: x86_64-conda-linux-gnu
# Running under: Rocky Linux 9.4 (Blue Onyx)
# 
# Matrix products: default
# BLAS:   /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/lib/libRblas.so 
# LAPACK: /jhpce/shared/community/core/conda_R/4.5/R/lib64/R/lib/libRlapack.so;  LAPACK version 3.12.1
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
#   [1] here_1.0.1                  escheR_1.8.0               
# [3] ggplot2_3.5.2               sessioninfo_1.2.3          
# [5] SpatialExperiment_1.18.1    SingleCellExperiment_1.30.1
# [7] SummarizedExperiment_1.38.1 Biobase_2.68.0             
# [9] GenomicRanges_1.60.0        GenomeInfoDb_1.44.0        
# [11] IRanges_2.42.0              S4Vectors_0.46.0           
# [13] BiocGenerics_0.54.0         generics_0.1.4             
# [15] MatrixGenerics_1.20.0       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] tidyselect_1.2.1        viridisLite_0.4.2       dplyr_1.1.4            
# [4] vipor_0.4.7             farver_2.1.2            viridis_0.6.5          
# [7] digest_0.6.37           rsvd_1.0.5              lifecycle_1.0.4        
# [10] cluster_2.1.8.1         Cairo_1.6-2             magrittr_2.0.3         
# [13] compiler_4.5.0          rlang_1.1.6             tools_4.5.0            
# [16] S4Arrays_1.8.0          labeling_0.4.3          scatterplot3d_0.3-44   
# [19] DelayedArray_0.34.1     RColorBrewer_1.1-3      abind_1.4-8            
# [22] BiocParallel_1.42.0     withr_3.0.2             grid_4.5.0             
# [25] beachmat_2.24.0         colorspace_2.1-1        scales_1.4.0           
# [28] iterators_1.0.14        dichromat_2.0-0.1       cli_3.6.5              
# [31] crayon_1.5.3            ragg_1.4.0              httr_1.4.7             
# [34] rjson_0.2.23            readxl_1.4.5            scuttle_1.18.0         
# [37] ggbeeswarm_0.7.2        parallel_4.5.0          cellranger_1.1.0       
# [40] XVector_0.48.0          vctrs_0.6.5             Matrix_1.7-3           
# [43] jsonlite_2.0.0          BiocSingular_1.24.0     GetoptLong_1.0.5       
# [46] BiocNeighbors_2.2.0     ggrepel_0.9.6           irlba_2.3.5.1          
# [49] clue_0.3-66             beeswarm_0.4.0          scater_1.36.0          
# [52] systemfonts_1.2.3       magick_2.8.6            foreach_1.5.2          
# [55] glue_1.8.0              codetools_0.2-20        cowplot_1.1.3          
# [58] Polychrome_1.5.4        shape_1.4.6.1           gtable_0.3.6           
# [61] UCSC.utils_1.4.0        ScaledMatrix_1.16.0     ComplexHeatmap_2.24.0  
# [64] tibble_3.2.1            pillar_1.10.2           GenomeInfoDbData_1.2.14
# [67] circlize_0.4.16         R6_2.6.1                textshaping_1.0.1      
# [70] doParallel_1.0.17       rprojroot_2.0.4         lattice_0.22-7         
# [73] png_0.1-8               Rcpp_1.0.14             gridExtra_2.3          
# [76] SparseArray_1.8.0       pkgconfig_2.0.3         GlobalOptions_0.1.2 
