# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5
#Load libraries
library(SingleCellExperiment)
library(here)
library(scuttle)
library(Matrix)
library(sessioninfo)
library("RcppML",lib.loc = "/users/rphillip/R/4.3.x") #Need development version for singlet. 
library("singlet",lib.loc = "/users/rphillip/R/4.3.x")

###### Single nucleus RNA-sequencing object
rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

#Run cross validation
message("Running NMF - ", Sys.time())
Sys.time()
options(RcppML.threads=4)
x <- RcppML::nmf(assay(rat_sce,"logcounts"),
                 k=53,
                 tol = 1e-06,
                 maxit = 1000,
                 verbose = T,
                 L1 = 0.1,
                 seed = 908,
                 mask_zeros = FALSE,
                 diag = TRUE,
                 nonneg = TRUE)


#Save the cross validation results
saveRDS(x, file = here("processed-data","rat_NMF","NMF_Results_k53.Rds"))

#sesion info
print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
