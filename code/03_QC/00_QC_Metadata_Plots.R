#Goal: Plot metrics as violins + on top of tissue. 
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.4.x
#code modified from https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/02_xenium_qc/00_plot_metrics_on_tissue.R

library(SpatialExperiment)
library(sessioninfo)
library(scattermore)
library(patchwork)
library(tidyverse)
library(escheR)
library(scater)
library(scran)
library(here)


#Read in the RDS file from 01_build_spe. At the moment, this contains only sections from Br6660 (runs 1-5)
spe <- readRDS(here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))

spe

###DEFINE FUNCTIONS
#Modify the plotColData function 
plot_colData_nac <- function(object,x,y,color){
  plotColData(object = object, x = x , y = y, colour_by = color) +
    ggtitle(y) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5),
          legend.position = "none") +
    stat_summary(fun = median, 
                 fun.min = median, 
                 fun.max = median,
                 geom = "crossbar", 
                 width = 0.3)
}

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
    plist[[i]] <- p
  }
  return(plist)
}

###

#Make sure the plot variables that could be numeric, are actually factors
spe$Slide_ID <- as.factor(spe$Slide_ID)
spe$Slide_Sample_Reagent_Lot <- as.factor(spe$Slide_Sample_Reagent_Lot)
spe$Decoding_Reagent_B_Lot <- as.factor(spe$Decoding_Reagent_B_Lot)
spe$Decoding_Reagent_A_Lot <- as.factor(spe$Decoding_Reagent_A_Lot)
spe$Human_Brain_Add_On_Lot <- as.factor(spe$Human_Brain_Add_On_Lot)


#Define some metricsa adn variables to plot. 
metrics_to_plot <- c("total_counts","unassigned_codeword_counts",
                     "cell_area","nucleus_area","transcript_counts")

plot_vars <- c("Sample", "Slide_ID", "Slide_Sample_Reagent_Lot",
               "Decoding_Reagent_B_Lot", "Decoding_Reagent_A_Lot",
               "Human_Brain_Add_On_Lot", "Custom_Panel_Lot", "Xenium_Instrument")

for(i in metrics_to_plot){
  message(paste("Plotting", i, "-", Sys.time()))
  
  plot_list <- lapply(plot_vars, function(j) {
    plot_colData_nac(object = spe, y = i, x = j, color = j)
  })
  
  # Combine with patchwork
  combined_plot <- wrap_plots(plot_list, ncol = 4) + 
    plot_annotation(title = paste(i, "QC Violin Plots"))
  
  # Save to file
  ggsave(filename = here("plots", 
                         "03_qc", 
                         paste0(i, 
                                "_violins_combined.png")),
         plot = combined_plot,
         width = 16, 
         height = 16, 
         dpi = 300)
}

message(paste("Moving to plot metrics on tissue -",Sys.time())) 

#on tissue
for(i in metrics_to_plot){
  message(paste("Plotting",i,"-",Sys.time()))
  pdf(here("plots", "03_qc", paste0(i,"_ontissue.pdf")))
  tissue_plots <- plot_coldata_on_tissue(spe, i)
  lapply(tissue_plots, print)
  dev.off()
}

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
