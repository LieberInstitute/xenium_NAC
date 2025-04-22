#Goal: Calculate QC metrics and explore. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.4.x
#Modified from:
 # https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/f3b87d7a0356e716468b304196f5cd362fd91515/code/analysis/02_xenium_qc/01_threshold_outliers.R
 # https://github.com/LieberInstitute/spatialAmygdala/blob/32108cdc145bb83aa73822c862dc28299a29b47c/code/Xenium/03_quality_control/01_perCellQC.R
library(SpatialExperiment)
library(scattermore)
library(scuttle)
library(scater)
library(escheR)
library(scran)
library(dplyr)
library(here)

#Load function to plot column data on top of tissue.
plot_outliers_on_tissue <- function(x, outliers){
  plist <- vector(mode = "list",length = length(outliers))
  for(i in seq_along(outliers)){
    outlier_col <- outliers[i]
    #Make the plot
    p <- make_escheR(x) |>
      add_fill(var = outlier_col) +
      ggtitle(paste0(
        sum(colData(x)[[outlier_col]], na.rm = TRUE), " outliers"
      )) +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_scattermore()
    
    plist[[i]] <- p
  }
  return(plist)
}
  
#Read in the RDS file from 01_build_spe. At the moment, this contains only sections from Br6660 (runs 1-5)
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))

spe


#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel

#Add QC metrics 
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        GEX = is_GEX))

colnames(colData(spe))

# Remove empty cells
empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]
spe <- spe[, colSums(counts(spe)) > 0]

all_outlier_ids <- c(empty_cells)

#How do nuclei/cell area correlate with counts + detected? 
#Nuc area by UMIs
png(here("plots","03_qc","Nucleus_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "sum")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","Nucleus_area_by_detected.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "detected")
dev.off()

#Nuc area by total_counts
png(here("plots","03_qc","Nucleus_area_by_totalcounts.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "total_counts")
dev.off()

#Cell area by UMIs
png(here("plots","03_qc","cell_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "sum")
dev.off()

#Cell area by nGenes
png(here("plots","03_qc","cell_area_by_nGenes.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "detected")
dev.off()

#Cell area by total_counts
png(here("plots","03_qc","Cell_area_by_totalcounts.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "total_counts")
dev.off()


#Split these plots by sample using facet_wrap
coldata_df <- as.data.frame(colData(spe))

#Nuc area by sum
png(here("plots","03_qc","nucleus_area_by_sum_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = nucleus_area, y = sum)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","nucleus_area_by_nGenes_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = nucleus_area, y = detected)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Cell area by sum
png(here("plots","03_qc","cell_area_by_sum_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = cell_area, y = sum)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","cell_area_by_nGenes_split.png"),height = 2500, width = 2500,res = 300)
ggplot(coldata_df, aes(x = cell_area, y = detected)) +
  geom_point(alpha = 0.5, 
             size = 0.7) +
  facet_wrap(~Sample, 
             scales = "free")
dev.off()

#There are NA values in unassigned, negative controls, and neg codewords ?
#For each sample identify which cells are in top 1% of these three metrics
#Input the data into a table
outlier_mat <- as.data.frame(matrix(nrow = 5,
                                    ncol = length(unique(spe$Sample))))

colnames(outlier_mat) <- unique(spe$Sample)
rownames(outlier_mat) <- c("unassigned_outlier",
                           "negProbe_outlier",
                           "negCodeword_outlier",
                           "detected_outlier",
                           "total_counts_outlier") 


for(i in unique(spe$Sample)){
  print(sprintf("------------%s------------", unique(spe$Sample)[i]))
  spe_sub <- spe[,spe$Sample == i]
  
  #unassigned, negProbe/Codeword will be calculated with quantiles
  unassigned_value <- quantile(x = spe_sub$subsets_unassigned_percent,.99,na.rm = TRUE)
  spe_sub$unassigned_outlier <- spe_sub$subsets_unassigned_percent >= unassigned_value
  print(unassigned_value)
  outlier_mat["unassigned_outlier",i] <- sum(spe_sub$unassigned_outlier,na.rm = TRUE)
  
  negProbe_value <- quantile(x = spe_sub$subsets_negProbe_percent,.99,na.rm = TRUE)
  spe_sub$negProbe_outlier <- spe_sub$subsets_negProbe_percent >= negProbe_value
  print(negProbe_value)
  outlier_mat["negProbe_outlier",i] <- sum(spe_sub$negProbe_outlier,na.rm = TRUE)
  
  negCodeword_value <- quantile(x = spe_sub$subsets_negCodeword_percent,.99,na.rm = TRUE)
  spe_sub$negCodeword_outlier <- spe_sub$subsets_negCodeword_percent >= negCodeword_value
  print(negCodeword_value)
  outlier_mat["negCodeword_outlier",i] <- sum(spe_sub$negCodeword_outlier,na.rm = TRUE)
  
  #Calculate sum and detected outliers with isOutlier 
  spe_sub$detected_outlier <- isOutlier(spe_sub$detected,type = "lower")
  outlier_mat["detected_outlier",i] <- sum(spe_sub$detected_outlier)
  
  spe_sub$total_counts_outlier <- isOutlier(spe_sub$total_counts,type = "lower")
  outlier_mat["total_counts_outlier",i] <- sum(spe_sub$total_counts_outlier)
  
  pdf(here("plots","03_qc","Sample_Specific_outliers",paste0(i,"_outliers.pdf")))
  tissue_plots <- plot_outliers_on_tissue(x = spe_sub,
                                          outliers = c("unassigned_outlier",
                                                       "negProbe_outlier",
                                                       "negCodeword_outlier",
                                                       "detected_outlier",
                                                       "total_counts_outlier"))
  lapply(tissue_plots, print)
  dev.off()
  
  # aggregate all outliers by taking the union
  spe_sub$is_outlier <- spe_sub$unassigned_outlier | spe_sub$negProbe_outlier | spe_sub$negCodeword_outlier | spe_sub$detected_outlier | spe_sub$total_counts_outlier
  outlier_ids <- ccolnames(spe_sub)[spe_sub$is_outlier]
  all_outlier_ids <- c(all_outlier_ids,outlier_ids)
}

#How many total outliers? 
print(length(all_outlier_ids))


#Save the outliers
save(all_outlier_ids,file = here("processed-data","Xenium_QC","all_outlier_ids.rda"))

#save the mat for the number of outliers for each metric for each sample
save(outlier_mat,file = here("processed-data","Xenium_QC","outlier_matrix.rda"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
