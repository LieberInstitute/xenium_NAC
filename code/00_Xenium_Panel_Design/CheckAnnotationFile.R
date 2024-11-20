#Goal for this script: Check that annotations file has cell type assigned to the correct barcode 
#cd ~/NAc_Xenium_Panel
#module load r_nac


library(SingleCellExperiment)
library(scran)
library(scater)
library(here)

#load the sce object. 
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")


#read in the annotations file 
annotations <- read.csv(here("Panel_Design_Files","annotations.csv"))
head(annotations)

#Check if barcodes in annotations file are in the sce object
message("Check if barcodes in annotations file are in the sce object")
all(annotations$barcode %in% rownames(colData(sce)))

stopifnot(all(annotations$barcode %in% rownames(colData(sce))))

all(annotations$barcode %in% colnames(sce))

stopifnot(all(annotations$barcode %in% colnames(sce)))

#Checking that barcodes all have the correct annotation. 
#Create an empty vector. 
empty_vec <- vector(length = nrow(annotations))

#Check 
for(i in 1:nrow(annotations)){
  #Pull barcode and annotation 
  barcode <- annotations[i,"barcode"]
  annotation <- annotations[i,"annotation"]
  
  #Use barcode to find the cell within the metadata of the sce object
  ground_truth <- colData(sce)[barcode,"CellType.Final"]

  #Does the annotation within the annotations dataframe = the celltype in the coldata of the object? 
  empty_vec[i] <- annotation == ground_truth
}

message("Are all values within the vector TRUE?")
table(empty_vec)

all(empty_vec)

message("Are any of the barcodes duplicated in the annotations file?")
any(duplicated(annotations$barcode))

message("Are the number of unique barcodes the same as the number of rows in the annotation dataframe?")
length(unique(annotations$barcode)) == nrow(annotations)


sessioninfo::session_info()
