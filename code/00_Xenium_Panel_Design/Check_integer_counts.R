#Xenium Panel Designer requires un-normalized counts. The guide suggests doublechecking this. 
#While potentially overkill, probably good to do the check. Running in job because the matrices are large. 

library(SingleCellExperiment)
library(sessioninfo)
library(scater)
library(scran)
library(here)


#load the sce object. 
sce <- readRDS("/dcs04/lieber/marmaypag/spatialNac_LIBD4125/spatial_NAc/processed-data/12_snRNA/sce_CellType_noresiduals.Rds")

sce

stopifnot(identical(rownames(colData(sce)),colnames(sce)))


#Pull counts matrix. 
matrix.1 <- counts(sce)

matrix.1 <- c(as.matrix(matrix.1))

#Convert with as.integer and create second matrix 
matrix.2 <- as.integer(counts(sce))

#check the class. matrix.1 should be numeric and matrix.2 should be integer
class(matrix.1)

class(matrix.2)

#check dims, head,and tail of matrix.1
length(matrix.1)

head(matrix.1)

tail(matrix.1) 

#check dims, head,and tail of matrix.2 
length(matrix.2)

head(matrix.2)

tail(matrix.2)


#Check with all.equal according to https://www.10xgenomics.com/analysis-guides/creating-single-cell-references-for-xenium-custom-panel-design-from-seurat-or-anndata
message("Are all values within counts matrix integers?")
all.equal(matrix.1,matrix.2) 


session_info()
