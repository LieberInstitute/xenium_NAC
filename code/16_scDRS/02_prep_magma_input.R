# Prepare MAGMA inputs from LDSC-munged GWAS summary statistics.
#
# 1. Builds a MAGMA gene location file keyed by the Visium-HD h5ad var names.
#    NOTE on genome builds: the Visium-HD object is GRCh38, but MAGMA must run
#    in hg19 to match the GWAS sumstats and the 1000G reference. The build
#    only affects SNP->gene assignment (all hg19 here); the handoff to the
#    expression object is by gene name. hg19 coordinates are looked up by
#    Ensembl ID (stable across builds) first, then by symbol, covering ~95%
#    of h5ad genes. Restricting to h5ad genes also keeps scDRS top-1000 gene
#    sets from being diluted by genes absent from the object.
# 2. Writes per-trait p-value files (SNP, P, N). The LDSC-munged files carry
#    Z rather than P, so P is derived as 2 * pnorm(-|Z|).
#
# Usage: sbatch 02_prep_magma_input.sh

library(here)
library(data.table)

manifest <- fread(here("processed-data", "16_scDRS", "trait_manifest.tsv"))
out_dir <- here("processed-data", "16_scDRS", "magma")

## ---- 1. Gene location file (hg19 coords, keyed by h5ad gene names) ------
h5ad <- here("processed-data", "HD_Full_Analysis", "h5ad", "VHD_sfe_counts.h5ad")
var_name <- as.character(rhdf5::h5read(h5ad, "var/_index"))
var_id <- as.character(rhdf5::h5read(h5ad, "var/ID"))

gene_meta <- fread("/dcs04/lieber/shared/statsgen/LDSC/base/gene_meta_hg19.txt")
setnames(gene_meta, c(
  "ensembl", "ensembl_v", "chr", "start", "end",
  "symbol", "gc", "type", "entrez"
))
gene_meta <- gene_meta[chr %in% as.character(1:22)]
# Longest span first, for de-duplication of IDs/symbols
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

## ---- 2. Per-trait p-value files ------------------------------------------
# NOTE: no SNP location file is written here. MAGMA's --annotate reads SNP
# locations directly from the merged 1000G PLINK .bim in step 03, so a
# separate .snp.loc would be unused.
for (i in seq_len(nrow(manifest))) {
  key <- manifest$trait_key[i]
  # Some files have a space-separated header over tab-separated data, which
  # breaks fread's auto-detection; parse the header separately.
  hdr <- strsplit(readLines(manifest$file_path[i], n = 1), "\\s+")[[1]]
  ss <- fread(manifest$file_path[i], skip = 1, header = FALSE)
  stopifnot(ncol(ss) == length(hdr))
  setnames(ss, hdr)

  stopifnot(all(c("SNP", "N", "Z") %in% names(ss)))
  ss <- ss[!is.na(Z) & !is.na(N) & !is.na(SNP)]
  # MAGMA requires unique SNP IDs; drop duplicates, keeping the first.
  n_dup <- sum(duplicated(ss$SNP))
  if (n_dup > 0) {
    ss <- ss[!duplicated(SNP)]
  }
  ss[, P := 2 * pnorm(-abs(Z))]
  # MAGMA dislikes P == 0; floor at smallest representable double
  ss[P < .Machine$double.xmin, P := .Machine$double.xmin]
  ss[, N := round(N)]

  fwrite(
    ss[, .(SNP, P, N)],
    file.path(out_dir, "pval", paste0(key, ".pval.txt")),
    sep = "\t"
  )
  message(key, ": ", nrow(ss), " SNPs", if (n_dup > 0) paste0(" (", n_dup, " duplicate SNPs dropped)") else "")
}

sessioninfo::session_info()
