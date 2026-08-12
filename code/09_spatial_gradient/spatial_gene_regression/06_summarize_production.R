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

plots <- val("--plot-dir")

prod <- file.path(out, "production")

pplots <- file.path(plots, "production")

dir.create(prod, FALSE, TRUE)

dir.create(pplots, FALSE, TRUE)

plan <- readRDS(file.path(out, "checkpoints", "production_plan.rds"))

paths <- file.path(out, "checkpoints", "production", paste0(plan$manifest$checkpoint_stem, ".rds"))

missing <- paths[!file.exists(paths)]

if (length(missing)) stop("Production checkpoints missing: ", length(missing))

obj <- lapply(paths, readRDS)

if (any(vapply(
  obj, function(z) !identical(z$identity, plan$identity) || !identical(z$status, "complete"),
  logical(1)
))) {
  stop("Stale/incomplete production checkpoint")
}

r <- axyd_bind(lapply(obj, `[[`, "results"))

if (nrow(r) != plan$n_planned_combined || anyDuplicated(r[c("Gene", "CellType", "model_id")]) || !setequal(
  r$Gene,
  plan$manifest$Gene
) || !setequal(r$CellType, plan$celltypes)) {
  stop("Incomplete 366 x 20 x 2 production universe")
}

r$p_BH_tested_model_family <- NA_real_

r$p_Bonferroni_planned_model_family_7320 <- NA_real_

r$p_Bonferroni_planned_combined_family_14640 <- NA_real_

for (m in c("M0", "M1")) {
  ix <- r$model_id == m & r$model_status == "tested" & is.finite(r$p_AP_nominal)
  r$p_BH_tested_model_family[ix] <- p.adjust(r$p_AP_nominal[ix], "BH")
  r$p_Bonferroni_planned_model_family_7320[ix] <- pmin(1, r$p_AP_nominal[ix] * plan$n_planned_per_model)
  r$p_Bonferroni_planned_combined_family_14640[ix] <- pmin(1, r$p_AP_nominal[ix] * plan$n_planned_combined)
}

r$significant_BH_model_family <- r$p_BH_tested_model_family < 0.05

r$significant_Bonferroni_7320 <- r$p_Bonferroni_planned_model_family_7320 < 0.05

r$significant_Bonferroni_14640 <- r$p_Bonferroni_planned_combined_family_14640 < 0.05

r <- r[order(r$model_id, r$planned_test_index), ]

m0 <- r[r$model_id == "M0", ]

m1 <- r[r$model_id == "M1", ]

a <- m0

b <- m1

names(a)[!names(a) %in% c("Gene", "CellType")] <- paste0("M0_", names(a)[!names(a) %in% c("Gene", "CellType")])

names(b)[!names(b) %in% c("Gene", "CellType")] <- paste0("M1_", names(b)[!names(b) %in% c("Gene", "CellType")])

cmp <- merge(a, b, by = c("Gene", "CellType"), all = TRUE, sort = FALSE)

cmp$delta_beta_AP <- cmp$M1_beta_AP_per_mm - cmp$M0_beta_AP_per_mm

cmp$SE_inflation_ratio <- cmp$M1_SE_AP / cmp$M0_SE_AP

cmp$direction_agreement <- sign(cmp$M0_beta_AP_per_mm) == sign(cmp$M1_beta_AP_per_mm)

runtime <- data.frame(Gene = vapply(obj, `[[`, character(1), "Gene"), elapsed_seconds = vapply(
  obj, `[[`,
  numeric(1), "elapsed_seconds"
), peak_recorded_R_memory_mb = vapply(obj, `[[`, numeric(1), "peak_recorded_R_memory_mb"))

status <- as.data.frame(with(r, table(model_id, model_status, useNA = "ifany")))

names(status) <- c("model_id", "model_status", "n")

sig <- do.call(rbind, lapply(c("M0", "M1"), function(m) {
  z <- r[r$model_id == m, ]
  data.frame(
    model_id = m, n_planned = nrow(z), n_tested = sum(z$model_status == "tested"), n_not_tested = sum(z$model_status ==
      "not_tested"), n_model_failed = sum(z$model_status == "model_failed"), BH_significant = sum(z$significant_BH_model_family,
      na.rm = TRUE
    ), Bonferroni_7320_significant = sum(z$significant_Bonferroni_7320, na.rm = TRUE),
    Bonferroni_14640_significant = sum(z$significant_Bonferroni_14640, na.rm = TRUE)
  )
}))

axyd_csv(plan$manifest, file.path(prod, "00_production_gene_manifest.csv"))

axyd_csv(r, file.path(prod, "01_all_M0_M1_results.csv"))

axyd_csv(m0, file.path(prod, "02_M0_global_results.csv"))

axyd_csv(m1, file.path(prod, "03_M1_global_results.csv"))

axyd_csv(cmp, file.path(prod, "04_M0_vs_M1_comparison.csv"))

axyd_csv(status, file.path(prod, "05_model_status_summary.csv"))

axyd_csv(runtime, file.path(prod, "06_runtime_and_memory.csv"))

axyd_csv(sig, file.path(prod, "07_significance_summary.csv"))

png(file.path(pplots, "01_production_significance_counts.png"), 1800, 1200, res = 200)

mat <- t(as.matrix(sig[c("BH_significant", "Bonferroni_7320_significant", "Bonferroni_14640_significant")]))

barplot(mat,
  beside = TRUE, names.arg = sig$model_id, col = c("#0072B2", "#D55E00", "#009E73"), ylab = "Number significant",
  main = "Production multiplicity sensitivity"
)

legend("topright", rownames(mat), fill = c("#0072B2", "#D55E00", "#009E73"), bty = "n")

dev.off()

png(file.path(pplots, "02_M0_vs_M1_AP_coefficients.png"), 1600, 1500, res = 200)

z <- cmp[cmp$M0_model_status == "tested" & cmp$M1_model_status == "tested", ]

plot(z$M0_beta_AP_per_mm, z$M1_beta_AP_per_mm,
  pch = 16, cex = 0.45, col = rgb(0, 0.45, 0.7, 0.35), xlab = "M0 AP beta/mm",
  ylab = "M1 AP beta/mm", main = "All tested Gene x CellType combinations"
)

abline(0, 1, lty = 2, col = "red")

grid()

dev.off()

png(file.path(pplots, "03_production_nominal_pvalue_QQ.png"), 1800, 900, res = 180)

par(mfrow = c(1, 2))

for (m in c("M0", "M1")) {
  p <- r$p_AP_nominal[r$model_id == m & r$model_status == "tested"]
  qqplot(-log10(ppoints(length(p))), -log10(sort(p)),
    pch = 16, cex = 0.4, xlab = "Expected -log10(p)",
    ylab = "Observed -log10(p)", main = paste(m, "nominal AP p-values")
  )
  abline(0, 1, lty = 2, col = "red")
}

dev.off()

plot_status <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(
    file.path(script_dir, "07_plot_D1_marker_effect_sizes.R"),
    "--output-dir",
    out,
    "--plot-dir",
    plots
  )
)

if (plot_status != 0) {
  stop("D1 marker effect-size plotting failed")
}

lines <- c(
  "# Full production summary", "", paste0("Generated: ", format(Sys.time(), tz = "UTC", usetz = TRUE)),
  "", "Frozen universe: 366 genes x 20 CellTypes = 7,320 planned rows per model; 14,640 combined M0/M1 rows.",
  "M0 is primary; M1 is the density-adjusted sensitivity model. Both use identical rows from the primary 2000-um unbalanced repeated-absolute-tile panel; every retained tile occurs in at least four selected physical slices.",
  "BH is calculated separately over successfully tested M0 and M1 p-values as specified. Conservative Bonferroni values use all 7,320 planned tests per model and all 14,640 planned tests combined.",
  "", paste0(
    "- ", sig$model_id, ": tested=", sig$n_tested, ", not_tested=", sig$n_not_tested, ", model_failed=",
    sig$n_model_failed, ", BH significant=", sig$BH_significant, ", Bonferroni(7320) significant=",
    sig$Bonferroni_7320_significant, ", Bonferroni(14640) significant=", sig$Bonferroni_14640_significant
  )
)

writeLines(lines, file.path(prod, "08_production_summary.md"))

writeLines("PASS", file.path(out, "checkpoints", "PRODUCTION_COMPLETE.txt"))

unlink(file.path(out, "RESULTS_STALE_AFTER_UNBALANCED_TILE_REDESIGN.md"))

cat("Full 366-gene x 20-CellType production summary complete\n")
