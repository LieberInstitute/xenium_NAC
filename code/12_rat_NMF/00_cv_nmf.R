# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
# module load conda_R/4.5

library(SingleCellExperiment)
library(sessioninfo)
library(orthogene)
library(scuttle)
library(Matrix)
library(Seurat)
library("RcppML",lib.loc = "/users/rphillip/R/4.3.x") #Need development version for singlet. 
library("singlet",lib.loc = "/users/rphillip/R/4.3.x")
library(here)

###### Single nucleus RNA-sequencing object
rat_obj <- readRDS(file = here("processed-data","Day_Lab_Panel","5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

#Run cross validation
message("Running cross validation - ", Sys.time())
cvnmf <- cross_validate_nmf(
  logcounts(rat_sce),
  ranks=c(30,40,50,60,70,80,90,100,125),
  n_replicates = 3,
  tol = 1e-03,
  maxit = 100,
  verbose = 3,
  L1 = 0.1,
  L2 = 0,
  threads = 0,
  test_density = 0.2
)

#Save the cross validation results
saveRDS(cvnmf, file = here("processed-data","rat_NMF","nmf_cv_results.Rds"))

#plot the results
pdf(here("plots","rat_nmf","cross_validation_results.pdf"))
plot(cvnmf)
dev.off()

print("Reproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
