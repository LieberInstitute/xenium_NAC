#Goal: Perform xenium QC
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.4.x
#code modified from https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/02_xenium_qc/

library(SpatialExperiment)
library(scattermore)
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
