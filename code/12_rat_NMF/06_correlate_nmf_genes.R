library(SingleCellExperiment)
library(pheatmap)
library(reshape2)
library(Seurat)
library(Matrix)
library(here)

###### Single nucleus RNA-sequencing object
rat_obj <- readRDS(file = here("processed-data", "Day_Lab_Panel", "5_NAc0069_integrated_clean_NAcOnly.rds"))
rat_sce <- as.SingleCellExperiment(x = rat_obj)
rat_sce

rm(rat_obj)

nmf_res <- readRDS(here("processed-data", "rat_NMF", "NMF_Results_k53.Rds"))

W <- nmf_res@w  # genes x patterns
H <- nmf_res@h  # patterns x cells

# ---- CONFIG --------------------------------------------------
top_n    <- 20
plot_dir <- here("plots", "rat_nmf")

# ---- Heatmap: top genes x chosen patterns (loading weight) -----------------
chosen_patterns <- c("nmf10","nmf11","nmf18","nmf24")
message("Building gene x pattern loading heatmap for: ",
        paste(chosen_patterns, collapse = ", "))
 
chosen_genes <- unique(unlist(lapply(chosen_patterns, function(pat) {
  names(sort(W[, pat], decreasing = TRUE))[seq_len(top_n)]
})))
 
heat_mat <- W[chosen_genes, chosen_patterns, drop = FALSE]
 
pdf(file.path(plot_dir,paste0(
  "NMF_top_genes_heatmap_", paste(chosen_patterns, collapse = "_"), ".pdf"
)))
pheatmap(
  heat_mat,
  color        = colorRampPalette(c("white", "firebrick"))(100),
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  fontsize_row = 7,
  main         = paste("Top genes -", paste(chosen_patterns, collapse = " vs "))
)
dev.off()
 
#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
 
