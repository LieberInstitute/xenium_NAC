#Goal: Calculate size factors based on nucleus and cell area
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5
#Code modified from: 
  #https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/Xenium/03_quality_control/01_perCellQC.R
  #https://github.com/LylaAtta123/normalization-analyses/blob/main/R/xenium.ipynb
library(SpatialExperiment)
library(sessioninfo)
library(ggplot2)
library(here)

spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_clean.Rds"))

spe

#Based on Atta et al. 2024, generate non-count based scaling factors
#Calculate both cell area and nucleus area values
spe$cell_area.sf <- spe$cell_area / median(spe$cell_area)
spe$nucleus_area.sf <- spe$nucleus_area / median(spe$nucleus_area)

#Generate histograms 
cell_area_hist <- ggplot(colData(spe),aes(x = cell_area.sf)) +
  geom_histogram(bins = 50,fill = "#70bbfe",color = "#FFF") +
  labs(x = "Scaling Factor",
       y = "Frequency",
       title = "Cell area") +
  theme_bw() +
  theme(plot.title = element_text(hjust= 0.5))
ggsave(plot = cell_area_hist,filename = here("plots","03_qc","cell_area_scaling_hist.pdf"))

nuc_area_hist <- ggplot(colData(spe),aes(x = nucleus_area.sf)) +
  geom_histogram(bins = 50,fill = "#70bbfe",color = "#FFF") +
  labs(x = "Scaling Factor",
       y = "Frequency",
       title = "Nucleus area") +
  theme_bw() +
  theme(plot.title = element_text(hjust= 0.5))
ggsave(plot = nuc_area_hist,filename = here("plots","03_qc","nucleus_area_scaling_hist.pdf"))

#Lines 58-61 straight from: https://github.com/LieberInstitute/spatialAmygdala/blob/devel/code/Xenium/03_quality_control/01_perCellQC.R
# normalize the counts by the nucleus and cell area scaling factors
assay(spe, "nucleus_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$nucleus_area.sf, transform="log", assay.type="counts")
assay(spe, "cell_normcounts") <- scuttle::normalizeCounts(spe, size.factors=spe$cell_area.sf, transform="log", assay.type="counts")


#Save spe with normalized counts
saveRDS(spe,here("processed-data","02_build_spe","SPEs","spe_NormCounts.Rds"))


###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()

