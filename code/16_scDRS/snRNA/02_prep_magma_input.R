# Prepare MAGMA inputs for snRNA-seq scDRS pipeline.
#
# 1. Builds a gene location file (hg19) keyed by the snRNA h5ad gene names.
#    With full transcriptome (~20k+ genes) the MAGMA gene sets will be richer
#    than the 366-gene Xenium panel version.
# 2. Symlinks shared resources from the parent scDRS pipeline (p-value files,
#    1000G reference, SNP location).
#
# Usage: sbatch 02_prep_magma_input.sh

library(here)
library(data.table)

out_dir <- here("processed-data", "16_scDRS", "snRNA", "magma")
parent_dir <- here("processed-data", "16_scDRS", "magma")

## ---- 1. Gene location file (hg19 coords, keyed by h5ad gene names) ------
h5ad <- here("processed-data", "16_scDRS", "snRNA", "h5ad", "snRNA_counts.h5ad")
var_name <- as.character(rhdf5::h5read(h5ad, "var/_index"))

# Try to read Ensembl IDs if available
var_id <- tryCatch(
  as.character(rhdf5::h5read(h5ad, "var/gene_id")),
  error = function(e) rep(NA_character_, length(var_name))
)

gene_meta <- fread("/dcs04/lieber/shared/statsgen/LDSC/base/gene_meta_hg19.txt")
setnames(gene_meta, c(
  "ensembl", "ensembl_v", "chr", "start", "end",
  "symbol", "gc", "type", "entrez"
))
gene_meta <- gene_meta[chr %in% as.character(1:22)]
gene_meta[, span := end - start]
setorder(gene_meta, -span)
by_id <- gene_meta[!duplicated(ensembl)]
by_symbol <- gene_meta[!is.na(symbol) & symbol != ""][!duplicated(symbol)]

idx_id <- match(var_id, by_id$ensembl)
idx_sym <- match(var_name, by_symbol$symbol)

gene_loc <- data.table(
  name  = var_name,
  chr   = fifelse(!is.na(idx_id), by_id$chr[idx_id], by_symbol$chr[idx_sym]),
  start = fifelse(!is.na(idx_id), by_id$start[idx_id], by_symbol$start[idx_sym]),
  end   = fifelse(!is.na(idx_id), by_id$end[idx_id], by_symbol$end[idx_sym])
)
message(
  "h5ad genes: ", length(var_name),
  "; matched by Ensembl ID: ", sum(!is.na(idx_id)),
  "; by symbol fallback: ", sum(is.na(idx_id) & !is.na(idx_sym)),
  "; unmatched (dropped): ", sum(is.na(gene_loc$chr))
)
gene_loc <- gene_loc[!is.na(chr)][!duplicated(name)]

fwrite(
  gene_loc,
  file.path(out_dir, "ref", "genes_hg19_symbol.gene.loc"),
  sep = "\t", col.names = FALSE
)
message("Gene loc: ", nrow(gene_loc), " genes")

## ---- 2. Symlink shared resources from parent pipeline -------------------
# P-value files
pval_link <- file.path(out_dir, "pval")
if (!file.exists(pval_link)) {
  file.symlink(file.path(parent_dir, "pval"), pval_link)
  message("Symlinked pval -> ", file.path(parent_dir, "pval"))
}

# 1000G merged bfile + SNP loc
ref_files <- c(
  "1000G.EUR.QC.merged.bed", "1000G.EUR.QC.merged.bim",
  "1000G.EUR.QC.merged.fam", "merge_list.txt"
)
for (f in ref_files) {
  target <- file.path(out_dir, "ref", f)
  source <- file.path(parent_dir, "ref", f)
  if (!file.exists(target) && file.exists(source)) {
    file.symlink(source, target)
    message("Symlinked ", f)
  }
}

sessioninfo::session_info()
