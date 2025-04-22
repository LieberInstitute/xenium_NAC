#Goal: Calculate QC metrics and explore. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
#module load conda_R/4.4.x
#Modified from:
 # https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/f3b87d7a0356e716468b304196f5cd362fd91515/code/analysis/02_xenium_qc/01_threshold_outliers.R
 # https://github.com/LieberInstitute/spatialAmygdala/blob/32108cdc145bb83aa73822c862dc28299a29b47c/code/Xenium/03_quality_control/01_perCellQC.R
library(SpatialExperiment)
library(scuttle)
library(scater)
library(scran)
library(here)

#Read in the RDS file from 01_build_spe. At the moment, this contains only sections from Br6660 (runs 1-5)
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))

spe
# class: SpatialExperiment 
# dim: 541 2452295 
# metadata(11): Samples Samples ... Samples Samples
# assays(1): counts
# rownames(541): ABCC9 ADAMTS12 ... DeprecatedCodeword_0381
# DeprecatedCodeword_0393
# rowData names(3): ID Symbol Type
# colnames(2452295): Br6660_NAc1_580_1 Br6660_NAc1_580_2 ...
# Br6660_Nac11_5580_201084 Br6660_Nac11_5580_201085
# colData names(32): Sample Barcode ... Race PrimaryDx
# reducedDimNames(0):
#   mainExpName: NULL
# altExpNames(0):
#   spatialCoords names(2) : x_centroid y_centroid
# imgData names(1): sample_id

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
# [1] "Sample"                       "Barcode"                     
# [3] "cell_id"                      "transcript_counts"           
# [5] "control_probe_counts"         "genomic_control_counts"      
# [7] "control_codeword_counts"      "unassigned_codeword_counts"  
# [9] "deprecated_codeword_counts"   "total_counts"                
# [11] "cell_area"                    "nucleus_area"                
# [13] "nucleus_count"                "segmentation_method"         
# [15] "sample_id"                    "Xenium_Run_ID"               
# [17] "Slide_ID"                     "Donor"                       
# [19] "CryoSection_Date"             "Slide_Sample_Reagent_Lot"    
# [21] "Decoding_Reagent_B_Lot"       "Decoding_Reagent_A_Lot"      
# [23] "Human_Brain_Add_On_Lot"       "Custom_Panel_Lot"            
# [25] "Day_1_Start"                  "Date_Bench_Assay_Completed"  
# [27] "Xenium_Instrument"            "Xenium_Start_Date"           
# [29] "Age"                          "Sex"                         
# [31] "Race"                         "PrimaryDx"                   
# [33] "sum"                          "detected"                    
# [35] "subsets_negProbe_sum"         "subsets_negProbe_detected"   
# [37] "subsets_negProbe_percent"     "subsets_negCodeword_sum"     
# [39] "subsets_negCodeword_detected" "subsets_negCodeword_percent" 
# [41] "subsets_unassigned_sum"       "subsets_unassigned_detected" 
# [43] "subsets_unassigned_percent"   "subsets_GEX_sum"             
# [45] "subsets_GEX_detected"         "subsets_GEX_percent"         
# [47] "total"  

#How do nuclei/cell area correlate with counts + detected? 
#Nuc area by UMIs
png(here("plots","03_qc","Nucleus_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "sum")
dev.off()

#Nuc area by nGenes
png(here("plots","03_qc","Nucleus_area_by_detected.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "nucleus_area", y = "detected")
dev.off()

#Cell area by UMIs
png(here("plots","03_qc","cell_area_by_sum.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "sum")
dev.off()

#Cell area by nGenes
png(here("plots","03_qc","cell_area_by_nGenes.png"),height = 1000, width = 1000,res = 300)
plotColData(spe,x = "cell_area", y = "detected")
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





#What are the distributions of negative probes, codewords, and unassigned 



