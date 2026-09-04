# Score curated drug-target gene sets (drugs of abuse + psychiatric drug
# classes) on Visium-HD cells and test enrichment across spatial domains and
# cell types.
#
# Scoring: per-cell mean logcounts over the target set (HDF5-friendly),
# z-scaled across cells. Testing: Kruskal-Wallis across groups + per-group
# Wilcoxon (group vs rest) with BH FDR.
#
# Usage: sbatch 03_curated_drug_sets.sh

library(SpatialFeatureExperiment)
library(SingleCellExperiment)
library(HDF5Array)
library(DelayedArray)
library(data.table)
library(here)

out_dir <- here("processed-data", "17_drug_analysis", "curated_sets")

## ---- Curated target sets ---------------------------------------------------
# Receptor/transporter targets of major drugs of abuse and psychiatric drug
# classes (brain-expressed direct targets; sources: IUPHAR/ChEMBL mechanisms).
drug_sets <- list(
  opioids            = c("OPRM1", "OPRD1", "OPRK1", "OPRL1"),
  psychostimulants   = c("SLC6A3", "SLC6A4", "SLC6A2", "SLC18A2", "TAAR1"),
  nicotine           = c("CHRNA4", "CHRNB2", "CHRNA5", "CHRNA3", "CHRNB4",
                         "CHRNA6", "CHRNA7"),
  alcohol            = c("GABRA1", "GABRA2", "GABRB1", "GABRB2", "GABRB3",
                         "GABRG2", "GRIN1", "GRIN2A", "GRIN2B", "GLRA1"),
  cannabinoids       = c("CNR1", "CNR2", "FAAH", "MGLL"),
  benzodiazepines    = c("GABRA1", "GABRA2", "GABRA3", "GABRA5", "GABRG2"),
  dissociatives      = c("GRIN1", "GRIN2A", "GRIN2B", "GRIN2D"),
  psychedelics       = c("HTR2A", "HTR2C", "HTR1A"),
  antipsychotics     = c("DRD2", "DRD3", "DRD4", "HTR2A", "HTR1A", "ADRA1A",
                         "HRH1", "CHRM1"),
  antidepressants    = c("SLC6A4", "SLC6A2", "HTR1A", "HTR2A", "HTR3A"),
  mood_stabilizers   = c("GSK3B", "IMPA1", "INPP1")
)

## ---- Load object -----------------------------------------------------------
sfe <- loadHDF5SummarizedExperiment(
  here("processed-data", "HD_Full_Analysis", "sfe_spatial_annotated")
)
stopifnot(all(c("Spatial_Domain", "labels") %in% colnames(colData(sfe))))

# Remove cells outside the NAc (Hypo and Excitatory spatial domains)
# and Excitatory cell-type labels from any remaining spatial domains
sd_exclude <- c("Hypo", "Excitatory")
label_exclude <- c("Excitatory")
keep <- !sfe$Spatial_Domain %in% sd_exclude & !sfe$labels %in% label_exclude
message(sprintf("Excluding %d cells (Spatial_Domain in {%s} or labels in {%s})",
                sum(!keep), paste(sd_exclude, collapse = ", "),
                paste(label_exclude, collapse = ", ")))
sfe <- sfe[, keep]

lc <- assay(sfe, "logcounts")

## ---- Per-cell scores -------------------------------------------------------
missing_genes <- lapply(drug_sets, setdiff, y = rownames(sfe))
for (s in names(missing_genes)) {
  if (length(missing_genes[[s]]) > 0) {
    message(s, ": missing from object: ",
            paste(missing_genes[[s]], collapse = ", "))
  }
}

scores <- sapply(drug_sets, function(genes) {
  genes <- intersect(genes, rownames(sfe))
  as.numeric(colMeans(lc[genes, , drop = FALSE]))
})
scores_z <- scale(scores,center = TRUE,scale = TRUE)
rownames(scores_z) <- colnames(sfe)

score_dt <- data.table(
  cell_id = colnames(sfe),
  Sample = sfe$Sample,
  Spatial_Domain = as.character(sfe$Spatial_Domain),
  labels = as.character(sfe$labels),
  scores_z
)
fwrite(score_dt, file.path(out_dir, "curated_set_scores_per_cell.tsv.gz"))

## ---- Group tests -----------------------------------------------------------
test_grouping <- function(group_col) {
  groups <- score_dt[[group_col]]
  res <- rbindlist(lapply(names(drug_sets), function(set) {
    x <- score_dt[[set]]
    kw_p <- kruskal.test(x, factor(groups))$p.value
    rbindlist(lapply(unique(groups), function(g) {
      in_g <- groups == g
      wt <- wilcox.test(x[in_g], x[!in_g])
      data.table(
        drug_set = set, group = g, kw_p = kw_p,
        median_in = median(x[in_g]), median_out = median(x[!in_g]),
        delta_median = median(x[in_g]) - median(x[!in_g]),
        wilcox_p = wt$p.value, n_cells = sum(in_g)
      )
    }))
  }))
  res[, wilcox_fdr := p.adjust(wilcox_p, method = "BH")]
  fwrite(res, file.path(out_dir, sprintf("curated_set_enrichment_%s.csv", group_col)))
  res
}

res_domain <- test_grouping("Spatial_Domain")
res_labels <- test_grouping("labels")

message("Top enriched (Spatial_Domain):")
print(head(res_domain[order(-delta_median)], 15))

sessioninfo::session_info()
