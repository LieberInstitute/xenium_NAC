library("RcppML", lib.loc = "/users/rphillip/R/4.3.x")
library(Seurat)
library(orthogene)
library(Matrix)
library(ggplot2)
library(here)

# ---- CONFIG --------------------------------------------------
seu_path      <- "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/xenium_NAC_collab/compiled_area_norm_L0.0_L0.9_k30_d30_r1.0_banksy_SeuObj.rds"
loadings_path <- here("processed-data", "rat_NMF", "NMF_Results_k53.Rds")
seurat_assay  <- "Xenium"
seurat_layer  <- "data"  
out_path      <- here("processed-data", "mouse_xenium", "proj_rat_nmf_on_mouse_xenium.rds")
seu_out_path  <- here("processed-data", "mouse_xenium", "compiled_area_norm_L0.0_L0.9_k30_d30_r1.0_banksy_SeuObj_NMFproj.rds")
plot_base_dir <- here("plots", "nmf_patterns_rat_mouse") 
# ---------------------------------------------------------------------------

message("Loading mouse Xenium Seurat object - ", Sys.time())
seu <- readRDS(seu_path)
seu

DefaultAssay(seu) <- seurat_assay
seu[[seurat_assay]] <- JoinLayers(seu[[seurat_assay]])
seu <- NormalizeData(seu, assay = seurat_assay, verbose = TRUE)

message("Loading rat NMF results - ", Sys.time())
nmf_res <- readRDS(loadings_path)
W <- nmf_res@w

# ---- Rat -> mouse 1:1 ortholog mapping (cached) -----------------------------
gene_map <- convert_orthologs(
  gene_df         = rownames(W),
  gene_input      = "listA",
  gene_output     = "dict",
  input_species   = "rat",
  output_species  = "mouse",
  non121_strategy = "drop_both_species",
  method          = "gprofiler"
)

keep <- rownames(W) %in% names(gene_map)
W <- W[keep, , drop = FALSE]
rownames(W) <- unlist(gene_map[rownames(W)])

# ---- Align genes between (now mouse-named) loadings and the Xenium panel ---
common_genes <- intersect(rownames(W), rownames(seu))
W_sub <- W[common_genes, , drop = FALSE]

# Pattern (NMF factor) names - fall back to generic names if W has none
pattern_names <- colnames(W_sub)
if (is.null(pattern_names)) pattern_names <- paste0("nmf", seq_len(ncol(W_sub)))
colnames(W_sub) <- pattern_names

message("Extracting log-normalized expression from Seurat object - ", Sys.time())
expr <- GetAssayData(seu, assay = seurat_assay, layer = seurat_layer)
expr <- expr[common_genes, , drop = FALSE]
expr <- as(expr, "CsparseMatrix")
stopifnot(identical(rownames(expr), rownames(W_sub)))

# patterns x cells
message("Projecting with RcppML::project (NNLS, non-negative) - ", Sys.time())
proj <- project(expr, w = W_sub, L1 = 0)
rownames(proj) <- pattern_names
colnames(proj) <- colnames(expr)

# remove rowSums == 0
proj1 <- proj[rowSums(proj) == 0, , drop = FALSE]
proj2 <- proj[rowSums(proj) != 0, , drop = FALSE]

# normalize
proj2 <- apply(proj2, 1, function(x){x / sum(x)})   # becomes cells x patterns
proj1 <- t(proj1)                                    # becomes cells x patterns

proj_final <- cbind(proj2, proj1)
proj_final <- proj_final[, match(rownames(proj), colnames(proj_final))]

message("Saving raw and normalized projections - ", Sys.time())
saveRDS(proj,       out_path)
saveRDS(proj_final, sub("\\.rds$", "_final.rds", out_path))

# ---- Project into the Seurat object ----------------------------------------
message("Adding NMF projection to Seurat object - ", Sys.time())
proj_t <- t(proj)[colnames(seu), , drop = FALSE]       
proj_final <- proj_final[colnames(seu), , drop = FALSE] 

# meta.data columns: nmf{i} (raw) and nmf{i}_norm (normalized)
meta_raw  <- as.data.frame(proj_t)
colnames(meta_raw) <- pattern_names
meta_norm <- as.data.frame(proj_final)
colnames(meta_norm) <- paste0(pattern_names, "_norm")

seu <- AddMetaData(seu, metadata = meta_raw)
seu <- AddMetaData(seu, metadata = meta_norm)

# reducedDims: NMF_Proj (raw) and NMF_Proj_norm (normalized)
emb_raw <- as.matrix(proj_t)
colnames(emb_raw) <- paste0("NMF_", seq_len(ncol(emb_raw)))
seu[["NMF_Proj"]] <- CreateDimReducObject(
  embeddings = emb_raw,
  loadings   = as.matrix(W_sub),
  key        = "NMF_",
  assay      = seurat_assay
)

emb_norm <- as.matrix(proj_final)
colnames(emb_norm) <- paste0("NMFnorm_", seq_len(ncol(emb_norm)))
seu[["NMF_Proj_norm"]] <- CreateDimReducObject(
  embeddings = emb_norm,
  loadings   = as.matrix(W_sub),
  key        = "NMFnorm_",
  assay      = seurat_assay
)


# ---- Plot each NMF pattern on each tissue section ---------------------------
dir.create(plot_base_dir, showWarnings = FALSE, recursive = TRUE)

fov_names <- Images(seu)

message("Plotting ", length(pattern_names), " NMF patterns across ",
        length(fov_names), " tissue sections - ", Sys.time())

for (pat in pattern_names) {
  feature <- paste0(pat, "_norm")
  pat_dir <- file.path(plot_base_dir, pat)
  dir.create(pat_dir, showWarnings = FALSE, recursive = TRUE)

  for (fov in fov_names) {
    p <- ImageFeaturePlot(
      seu,
      fov      = fov,
      features = feature,
      cols     = c("lightgrey", "darkred")
    ) 

    ggsave(
      filename = file.path(pat_dir, paste0(pat, "_", fov, ".png")),
      plot     = p,
      width    = 12,
      height   = 12,
      dpi      = 300
    )
  }
  message(">>> Finished ", pat, " (", length(fov_names), " sections) - ", Sys.time())
}

message("Saving Seurat object with NMF projection - ", Sys.time())
saveRDS(seu, seu_out_path)

#-------------------------- Reproducibility information -----------------------#
message("\nReproducibility information:")
Sys.time()
proc.time()
options(width = 120)
sessioninfo::session_info()
