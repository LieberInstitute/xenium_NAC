#Goal: Calculate QC metrics and explore. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.5
#Modified from:
 # https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/f3b87d7a0356e716468b304196f5cd362fd91515/code/analysis/02_xenium_qc/01_threshold_outliers.R
 # https://github.com/LieberInstitute/spatialAmygdala/blob/32108cdc145bb83aa73822c862dc28299a29b47c/code/Xenium/03_quality_control/01_perCellQC.R
 # https://pachterlab.github.io/voyager/articles/vig5_xenium.html#quality-control
library(SpatialExperiment)
library(scattermore)
library(tidyverse)
#library(Voyager)
library(scuttle)
library(scater)
library(escheR)
library(scran)
library(dplyr)
library(here)

#Plot colDAta on tissue via escheR. Thanks Cindy! 
plot_coldata_on_tissue <- function(x, column_name){
  plist <- list()
  sample <- unique(x$Sample)
  for (i in 1:length(sample)){
    x_sub <- x[, x$Sample == sample[i]]
    p <- make_escheR(x_sub) %>%
      add_fill(column_name)+
      ggtitle(sample[i])+
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_scattermore()
    plist[[sample[i]]] <- p
  }
  return(plist)
}


  
#Read in the RDS file from 01_build_spe. At the moment, this contains only sections from Br6660 (runs 1-5)
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))

spe

#Make the sample column a factor
spe$Sample <- factor(x = spe$Sample, levels = unique(spe$Sample))

levels(spe$Sample)

#Find genes for scuttle subsets 
is_neg <- stringr::str_detect(rownames(spe), "^NegControlProbe")
is_neg2 <- stringr::str_detect(rownames(spe), "^NegControlCodeword")
is_unassigned <- stringr::str_detect(rownames(spe), "^Unassigned")
is_anyneg <- is_neg | is_neg2 | is_unassigned
is_GEX <- rowData(spe)$Type == "Gene Expression" #This will help identify QC based on just the genes in panel


#Per: https://genomics.uci.edu/wp-content/uploads/sites/30/PDF_UCI_GRT_Hub_Xenium_Workshop_20240118-120fec8a1ce5134f.pdf
# "Negative Control Codewords are
# a random subset of codewords,
# with identical properties to gene
# codewords" 
#and
# "By definition, a call made to negative control codeword is an error"

#Unassigned probes help measure noise and off-target activity because they do not match any gene codewords

#Negative control probes are sequences that should not bind anything and therefore measure false positive rates and specificity

#Add QC metrics 
spe <- scuttle::addPerCellQCMetrics(spe, subsets = list(negProbe = is_neg,
                                                        negCodeword = is_neg2,
                                                        unassigned = is_unassigned,
                                                        any_neg = is_anyneg,
                                                        GEX = is_GEX))

colnames(colData(spe))

# Remove empty cells
empty_cells <- colnames(spe)[colSums(counts(spe)) == 0]
spe <- spe[, colSums(counts(spe)) > 0]
dim(spe)

#There are some cells within this object in which 100% of the reads are negative controls - Remove them
spe <- spe[,spe$subsets_any_neg_percent < 100]
dim(spe)


##########################################
#Will use the concatenated any_neg QC metric as preliminary analysis demonstrated that most true 
#outliers were odentified by percent negative probe while there were a few outliers identified by the 
#other metrics. This seems to work as a catch all. 
#Calculate quantiles for any_neg
any_neg99_cutoffs <- tapply(colData(spe)$subsets_any_neg_percent, colData(spe)$Sample, function(x){
  quantile(x,.99,na.rm = TRUE)
})

#Above generates any array. Convert to named numeric vector. 
any_neg99_cutoffs_num <- as.numeric(any_neg99_cutoffs)
names(any_neg99_cutoffs_num) <- names(any_neg99_cutoffs)

ref_values <- any_neg99_cutoffs_num[as.character(spe$Sample)]

spe$any_neg_99 <- spe$subsets_any_neg_percent > ref_values  # or any other condition
table(spe$any_neg_99)


#Write a csv with number of retained and cutoffs included. 
cbind(as.matrix(table(spe$Sample,spe$any_neg_99)),any_neg99_cutoffs) %>% 
  as.data.frame() %>% 
  rename(percentile_99 = any_neg99_cutoffs) %>% 
  mutate(Sample = rownames(.),
         CellsRetained = `FALSE`,
         Outliers = `TRUE`) %>%
  select(Sample,everything()) %>%
  `rownames<-`(NULL) %>%
  write.csv(file = here("processed-data","Xenium_QC","99th_Percentile_Outliers_Any_Neg.csv"),
            quote = FALSE)
  
#Plot the percentage of any negative probe by detected and sum
#detected
for(i in unique(spe$Sample)){
  spe_sub <- spe[,spe$Sample == i]
  p <- plotColData(spe_sub,x = "detected",y="subsets_any_neg_percent",color_by = "any_neg_99")
  ggsave(p, filename = here("plots","03_qc","Sample_Specific_Outliers",paste0(i,"_anyneg_99_detected.png")))
}

#Sum 
for(i in unique(spe$Sample)){
  spe_sub <- spe[,spe$Sample == i]
  p <- plotColData(spe_sub,x = "sum",y="subsets_any_neg_percent",color_by = "any_neg_99")
  ggsave(p, filename = here("plots","03_qc","Sample_Specific_Outliers",paste0(i,"_anyneg_99_sum.png")))
}


#########
##Use adaptive thresholds to identify problematic cells based on # of detected genes + total umis

#isoutlier for detected and sum
#detected
spe$detected_out <- isOutlier(
  spe$detected,
  log = TRUE,
  batch = as.factor(spe$Sample),
  nmads = 3,
  type = "lower"
)

table(spe$Sample,spe$detected_out)

p <- plotColData(spe, x = "Sample", y = "detected", color_by = "detected_out") +
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(p, filename = here("plots","03_qc","detected_outliers_violin.png"))

table(spe$Sample,spe$detected_out) %>% 
  as.data.frame.matrix() %>%
  mutate(Sample = rownames(.),
         CellsRetained = `FALSE`,
         Outliers = `TRUE`) %>%
  select(Sample,CellsRetained,Outliers) %>%
  `rownames<-`(NULL) %>%
  write.csv(file = here("processed-data","Xenium_QC","Detected_Outliers.csv"),
            quote = FALSE)

#sum
spe$sum_out <- isOutlier(
  spe$sum,
  log = TRUE,
  batch = as.factor(spe$Sample),
  nmads = 3,
  type = "lower"
)

p <- plotColData(spe, x = "Sample", y = "sum", color_by = "sum_out") +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
  ggsave(p, filename = here("plots","03_qc","sum_outliers_violin.png"))

table(spe$Sample,spe$sum_out) %>% 
  as.data.frame.matrix() %>%
  mutate(Sample = rownames(.),
         CellsRetained = `FALSE`,
         Outliers = `TRUE`) %>%
  select(Sample,CellsRetained,Outliers) %>%
  `rownames<-`(NULL) %>%
  write.csv(file = here("processed-data","Xenium_QC","Sum_Outliers.csv"),
            quote = FALSE)

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

#Concert spe to sfe to use some of the colData  plotting functions
#sfe <- toSpatialFeatureExperiment(spe)

#Pull coldata columns that include percent counts mapping to negProbes/Codewords/unassigned/GEX
#cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_percent$")]
#cols_use

#p <- plotColDataHistogram(sfe, cols_use, bins = 100)
#ggsave(filename = here("plots","03_qc","subsets_histogram.pdf"),plot = p)

#Just as in the Voyager Xenium QC workflow, the histrograms for negative control probes and unassigned are 
#centered on 0. Additionally, most of the cells have 100% of their counts derived from actual genes.
#Do for all but the GEX subset which is not centered on 0. 
#p2 <- plotColDataHistogram(sfe, cols_use[-5], bins = 100) + 
#  scale_x_log10() 
#ggsave(filename = here("plots","03_qc","subsets_histogram_logscale.pdf"),plot = p2)

#This actually looks pretty good with most cells falling under 1%. 
#Let's make violin plots of the subsets by Sample
#for(i in cols_use){
#  print(i)
#  if(i == "subsets_GEX_percent"){
#    p_i <- plotColData(spe, 
#                       x = "Sample", 
#                       y = i,
#                       color_by = "Sample") +
#      theme(axis.text.x = element_text(angle = 45, hjust = 1),
#            legend.position = "none")  +
#      stat_summary(fun = median, 
#                   fun.min = median, 
#                   fun.max = median,
#                   geom = "crossbar", 
#                   width = 0.3)
#    ggsave(filename = here("plots","03_qc",paste0(i,"_violin_by_sample.png")),plot = p_i)
#  }else{
#    p_i <- plotColData(spe, 
#                       x = "Sample", 
#                       y = i,
#                      color_by = "Sample") +
#      scale_y_log10() +
#      theme(axis.text.x = element_text(angle = 45, hjust = 1),
#            legend.position = "none")  +
#      stat_summary(fun = median, 
#                   fun.min = median, 
#                   fun.max = median,
#                   geom = "crossbar", 
#                   width = 0.3)
#    ggsave(filename = here("plots","03_qc",paste0(i,"_violin_by_sample.png")),plot = p_i) 
#  }
#}

#Several warnings generated above from cells that have 0 values for these QC measures. 

#Do the same but for sum or the total of each count. 
#cols_use <- names(colData(sfe))[str_detect(names(colData(sfe)), "_sum$")]
#cols_use

#p_sums <- plotColDataHistogram(sfe, cols_use, bins = 100, ncol = 2) +
#  scale_y_log10()
#ggsave(filename = here("plots","03_qc","sums_histograms.png"),plot = p_sums)

#Generate violins as above. 
#for(i in cols_use){
#  print(i)
#  if(i == "subsets_GEX_sum"){
#    p_i <- plotColData(spe, 
#                       x = "Sample", 
#                       y = i,
#                       color_by = "Sample") +
#      scale_y_log10() +
#      theme(axis.text.x = element_text(angle = 45, hjust = 1),
#            legend.position = "none")  +
#      stat_summary(fun = median, 
#                   fun.min = median, 
#                   fun.max = median,
#                   geom = "crossbar", 
#                   width = 0.3)
#    ggsave(filename = here("plots","03_qc",paste0(i,"_violin_by_sample.png")),plot = p_i)
#  }else{
#    p_i <- plotColData(spe, 
#                       x = "Sample", 
#                       y = i,
#                       color_by = "Sample") +
#      theme(axis.text.x = element_text(angle = 45, hjust = 1),
#            legend.position = "none")  +
#      stat_summary(fun = median, 
#                  fun.min = median, 
#                   fun.max = median,
#                   geom = "crossbar", 
#                   width = 0.3)
#    ggsave(filename = here("plots","03_qc",paste0(i,"_violin_by_sample.png")),plot = p_i)
#  }
#}

#What about cell size and nucleus area? Can these be used as a QC for the segmentation? 
#p <- plotColDataHistogram(sfe, c("cell_area", "nucleus_area"), scales = "free_x")
#ggsave(filename = here("plots","03_qc","cell_nucleus_area_freex.png"),plot = p)


p <- plotColData(spe,x = "Sample", y = "cell_area") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(filename = here("plots","03_qc","cell_area_by_Sample_violin.png"),plot = p)

p <- plotColData(spe,x = "Sample", y = "nucleus_area") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")  +
  stat_summary(fun = median, 
               fun.min = median, 
               fun.max = median,
               geom = "crossbar", 
               width = 0.3)
ggsave(filename = here("plots","03_qc","nucleus_area_by_Sample_violin.png"),plot = p)

##########
##Remove cells currently identified as low-quality. 
#Using detected and sum outliers, as well as any_neg 99th percentile. 
spe$discard <- spe$any_neg_99 | spe$detected_out | spe$sum_out


#Save discard metrics as a csv
table(spe$Sample,spe$discard) %>% 
  as.data.frame.matrix() %>%
  mutate(Sample = rownames(.),
         CellsRetained = `FALSE`,
         Outliers = `TRUE`,
         PercentRemoved = Outliers/(Outliers+CellsRetained)*100) %>%
  select(Sample,CellsRetained,Outliers,PercentRemoved) %>%
  `rownames<-`(NULL) %>%
  write.csv(file = here("processed-data","Xenium_QC","All_Outliers.csv"),
            quote = FALSE)


discard_tissue_plot <- plot_coldata_on_tissue(x = spe,column_name = "discard")
for(i in names(discard_tissue_plot)){
  print(i)
  ggsave(plot = discard_tissue_plot[[i]],
         filename = here("plots","03_qc","Sample_Specific_Outliers",
                         paste0(i,"_discard_outliers.png")),
         height = 16,
         width = 16)
}

#Save spe before removing low quality cells
message(paste0("Saving spe object with QC - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_withQC.Rds"))


#Remove all of the low quality cells
spe <- spe[,!spe$discard]
dim(spe)
spe

#Save cleaned SPE
message(paste0("Saving cleaned SPE object - ",Sys.time()))
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_clean.Rds"))


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
