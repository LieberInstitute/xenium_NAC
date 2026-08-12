#!/usr/bin/env Rscript
argv <- commandArgs(FALSE)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1])))

source(file.path(script_dir, "common.R"))

args <- commandArgs(TRUE)

val <- function(n) {
  i <- match(n, args)
  if (is.na(i) || i == length(args)) {
    stop("Missing ", n)
  }
  args[i + 1]
}

out <- val("--output-dir")

prep <- readRDS(file.path(out, "checkpoints", "prepared_data.rds"))

if (!file.exists(file.path(out, "checkpoints", "PILOT_ACCEPTED.txt"))) stop("Pilot has not been accepted")

cts <- sort(unique(as.character(prep$support_decisions$CellType[prep$support_decisions$grid_variant ==
  "primary"])))

if (length(prep$genes) != 366 || length(cts) != 20) stop("Frozen production universe must be 366 genes x 20 CellTypes")

manifest <- data.frame(task_index = seq_along(prep$genes), Gene = prep$genes, checkpoint_stem = sprintf(
  "%03d_%s",
  seq_along(prep$genes), gsub("[^A-Za-z0-9._-]", "_", prep$genes)
), stringsAsFactors = FALSE)

sources <- c(file.path(script_dir, "common.R"), normalizePath(sub("^--file=", "", grep("^--file=", argv,
  value = TRUE
)[1])), file.path(script_dir, "05_fit_production_gene.R"), file.path(script_dir, "06_summarize_production.R"))

identity <- axyd_hash_text(c(prep$identity$identity, vapply(sources, axyd_hash_file, character(1))))

plan <- list(
  identity = identity, prepared_identity = prep$identity$identity, manifest = manifest, celltypes = cts,
  n_genes = 366L, n_celltypes = 20L, n_planned_per_model = 7320L, n_planned_combined = 14640L, approved_utc = format(Sys.time(),
    tz = "UTC", usetz = TRUE
  ), approval = "explicit user approval after accepted D1 marker pilot",
  sources = sources
)

dir.create(file.path(out, "checkpoints", "production"), FALSE, TRUE)

axyd_csv(manifest, file.path(out, "checkpoints", "production_gene_manifest.csv"))

axyd_rds(plan, file.path(out, "checkpoints", "production_plan.rds"))

writeLines("PASS", file.path(out, "checkpoints", "PRODUCTION_APPROVED.txt"))

cat("Production plan frozen: 366 genes x 20 CellTypes x M0/M1\n")
