#Goal: Build raw SPE
#cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
#module load conda_R/4.5
#code modified from https://github.com/LieberInstitute/spatialDLPFC_SCZ_XENIUM/blob/devel/code/analysis/01_build_spe/01_build_spe.R

library(SingleCellExperiment)
library(SpatialExperiment)
library(sessioninfo)
library(tidyverse)
library(here)
library(escheR)


#Load sample_info dataframe that contains information about each sample including depth and directories. 
sample_info <- read.csv(here("processed-data","Sample_Info_NAcXenium_All.csv"))

#Add the full sample path to the sample_info dataframe
sample_info$full_data_path <- file.path(paste(here("raw-data","xenium"),sample_info$Output_Directory,sep = "/"))

#Add  a column of the depth, which is the last aspect of the Sample_ID
sample_info$Depth <- as.numeric(sub(".*_(\\d{3,4})$", "\\1", sample_info$Sample_ID))
  
all_spes <- vector(mode = "list",length = nrow(sample_info))
names(all_spes) <- sample_info$Sample_ID
for(i in 1:nrow(sample_info)){
  print(i)
  #Pull path to data and sample_ID
  sample_path <- sample_info[i,"full_data_path"]
  Sample <- sample_info[i,"Sample_ID"]
  
  print(sprintf("----------%s----------", Sample))
  
  sample_info_use <- sample_info[i,]
  
  #Paths to counts matrix and cell csv file
  counts_path <- here(sample_path,"cell_feature_matrix.h5")
  cell_info_path <- here(sample_path,"cells.csv.gz")
  
  #Creae sce object
  sce <- DropletUtils::read10xCounts(samples = counts_path, sample.names = Sample)
  counts(sce) <- methods::as(DelayedArray::realize(counts(sce)), "dgCMatrix") # Convert to delayed array
  
  cell_info <- vroom::vroom(cell_info_path) #reads in cell csv file much faster than standard functions
  
  #Add cell info to sce
  colData(sce) <- cbind(colData(sce),cell_info)
  spe <- toSpatialExperiment(sce,spatialCoordsNames = c("x_centroid","y_centroid"))
  rownames(spe) <- rowData(spe)$Symbol #No need to uniquify because Ensembl speicifed during panel generation
  
  #Need to assign donor-specific column names because we will merge everything after loop is finished. 
  colnames(spe) <- paste(Sample,rownames(cell_info),sep = "_")
  
  #Add in relevant metadata. 
  #First technical variables. 
  spe$Xenium_Run_ID <- sample_info_use$Xenium_Run_ID
  spe$Slide_ID <- sample_info_use$Slide_ID
  spe$Donor <- sample_info_use$Donor
  spe$CryoSection_Date <- sample_info_use$Date_of_Cryosection
  spe$Slide_Sample_Reagent_Lot <- sample_info_use$PN.1000460_Lot_Slides_Sample_Reagents
  spe$Decoding_Reagent_B_Lot <- sample_info_use$PN.1000625_Lot_Decoding_Reagent_Module_B
  spe$Decoding_Reagent_A_Lot <- sample_info_use$PN.1000624_Lot_Decoding_Reagent_Module_A
  spe$Human_Brain_Add_On_Lot <- sample_info_use$PN.1000599_Lot_Human_Brain_Add.On
  spe$Custom_Panel_Lot <- sample_info_use$PN.1000561_Lot_Add.On_Custom_Panel
  spe$Day_1_Start <- sample_info_use$Day_1_Start
  spe$Date_Bench_Assay_Completed <- sample_info_use$Date_Bench_Assay_Completed
  spe$Xenium_Instrument <- sample_info_use$Xenium_Instrument_Used
  spe$Xenium_Start_Date <- sample_info_use$Start_Date_Xenium_Run
  
  #Now biological variables
  #spe$sample_id <- Sample
  spe$Donor <- sample_info_use$Donor
  spe$Age <- sample_info_use$Age
  spe$Sex <- sample_info_use$Sex
  spe$Race <- sample_info_use$Race
  spe$PrimaryDx <- sample_info_use$PrimaryDx
  
  #Next step is to rotate the samples. 
  coords <- spatialCoords(spe)
  #metadata(spe)$original_coords <- coords
  
  if(Sample %in% c("Br6660_NAc1_580", "Br6660_NAc3_1580", "Br6660_NAc4_2080", 
                   "Br6660_NAc5_2580", "Br6660_NAc6_3080", "Br6660_NAc7_3580", 
                   "Br6660_NAc8_4580", "Br6660_NAc9_5080", "Br6660_Nac10_4080", 
                   "Br6660_Nac11_5580","Br6436_Nac1_650", "Br6436_Nac2_1150", 
                   "Br6436_Nac3_1650", "Br6436_Nac_4_2150", "Br6436_Nac_5_2650", 
                   "Br6436_Nac_6_3150", "Br6436_Nac_7_3660", "Br6436_Nac_8_4150", 
                   "Br6436_Nac_9_4650", "Br6436_Nac_10_5150")){
    #rotate 90 degrees COUNTER clockwise
    #To do this, first find the center point
    center_x <- mean(range(coords[,1])) 
    center_y <- mean(range(coords[,2]))
    
    #Calculate distance of points from the origin. 
    x0 <- coords[,1] - center_x
    y0 <- coords[,2] - center_y
    
    #Apply the rotation by making x = -y and y = x
    x_rot <- -y0
    y_rot <- x0
    
    #Now add x0+x_rot and y0+y_rot to translate everything back to original space. 
    x_final <- x_rot + center_x
    y_final <- y_rot + center_y
    
    #Convert the coordinates
    spatialCoords(spe) <- cbind(x_final,y_final)
  }else{
    if(Sample == "Br6436_Nac_11_5650"){
      #rotate 90 degrees clockwise
      #To do this, first find the center point
      center_x <- mean(range(coords[,1])) 
      center_y <- mean(range(coords[,2]))
      
      #Calculate distance of points from the origin. 
      x0 <- coords[,1] - center_x
      y0 <- coords[,2] - center_y
      
      #Apply the rotation by making x = -y and y = x
      x_rot <- y0
      y_rot <- -x0
      
      #Now add x0+x_rot and y0+y_rot to translate everything back to original space. 
      x_final <- x_rot + center_x
      y_final <- y_rot + center_y
      
      #Convert the coordinates
      spatialCoords(spe) <- cbind(x_final,y_final)
    }
  }
  
  #Plot total_counts on top of the tissue to visualize rotation. 
  x <- escheR::make_escheR(spe) |>
    escheR::add_fill("total_counts")
  ggsave(plot = x,
         filename = here("plots","02_build_spe","Rotation_check",
                         paste0(Sample,"_post_rotation.png")),
         height = 16, width = 18)
  
  
  all_spes[[Sample]] <- spe
  
  #clear memory
  rm(sce)
  rm(cell_info)
  gc() #garbage collection
}

#Combine all of the spes 
message(paste0("Combining the spes - ",Sys.time()))
all_spes <- do.call(cbind,all_spes)

message(paste0("Saving SPE - ",Sys.time()))
saveRDS(all_spes,here("processed-data","02_build_spe","SPEs","spe_raw.Rds"))

###Reproduciblity
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()



























