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

gene <- val("--gene")

resume <- "--resume" %in% args

prep <- readRDS(file.path(out, "checkpoints", "prepared_data.rds"))

plan <- readRDS(file.path(out, "checkpoints", "production_plan.rds"))

if (!identical(plan$prepared_identity, prep$identity$identity)) stop("Production plan/prepared identity mismatch")

j <- match(gene, plan$manifest$Gene)

if (is.na(j)) stop("Unknown production gene: ", gene)

checkpoint <- file.path(out, "checkpoints", "production", paste0(plan$manifest$checkpoint_stem[j], ".rds"))

if (resume && file.exists(checkpoint)) {
  q <- readRDS(checkpoint)
  if (identical(q$identity, plan$identity) && identical(q$status, "complete")) {
    cat("Reusing production checkpoint", gene, "\n")
    quit(save = "no")
  }
}

started <- Sys.time()

cc <- read.csv(gzfile(file.path(prep$input_dir, "primary_counts.csv.gz")), check.names = FALSE)[c(
  "tile_group_id",
  gene
)]

results <- list()

ri <- 0L

empty <- function(ct, m, status, reason, n = 0, nt = 0, ns = 0) {
  data.frame(
    planned_test_index = (j -
      1L) * plan$n_celltypes + match(ct, plan$celltypes), Gene = gene, CellType = ct, model_id = m,
    grid_origin_id = "primary", support_design = prep$config$support$design,
    minimum_slices_per_absolute_tile = as.integer(prep$config$support$minimum_slices_per_absolute_tile),
    n_slices = ns, n_tiles = nt, n_rows = n, AP_min_mm = NA, AP_max_mm = NA, AP_span_mm = NA, beta_AP_per_mm = NA,
    SE_AP = NA, CI_AP_lower = NA, CI_AP_upper = NA, test_statistic_AP = NA, df_AP = NA, p_AP_nominal = NA,
    AP_rate_ratio_per_mm = NA, AP_percent_rate_change_per_mm = NA, covariance_method = "CR2 clustered by physical slice; Satterthwaite df",
    inference_status = status, support_status = status, model_status = status, failure_reason = reason,
    beta_density_within = NA, SE_density_within = NA, beta_density_between = NA, SE_density_between = NA,
    pearson_AP_density_between = NA, spearman_AP_density_between = NA, VIF_AP = NA, condition_number = NA,
    collinearity_flag = "", AP_density_between_coefficient_correlation = NA, stringsAsFactors = FALSE
  )
}

for (ct in plan$celltypes) {
  dec <- prep$support_decisions[prep$support_decisions$grid_variant == "primary" & prep$support_decisions$CellType ==
    ct, ]
  d <- prep$rows[prep$rows$grid_variant == "primary" & prep$rows$CellType == ct, ]
  if (nrow(dec) != 1 || dec$status != "tested" || !nrow(d)) {
    for (m in c("M0", "M1")) {
      ri <- ri + 1L
      results[[ri]] <- empty(ct, m, "not_tested", if (nrow(dec)) {
        dec$failure_reason
      } else {
        "support_decision_missing"
      })
    }
    next
  }
  d$raw_count <- cc[[gene]][match(d$tile_group_id, cc$tile_group_id)]
  d <- d[order(d$AP_um, d$absolute_tile), ]
  d$absolute_tile <- droplevels(factor(d$absolute_tile))
  d$Sample <- factor(d$Sample, levels = unique(d$Sample[order(d$AP_um)]))
  if (anyNA(d$raw_count) || any(d$raw_count < 0) || any(d$exposure <= 0) || any(d$n_cells <= 0)) {
    stop("Production row reconciliation failed: ", gene, " / ", ct)
  }
  reason <- character()
  if (sum(d$raw_count > 0) < as.integer(prep$config$gene_eligibility$minimum_positive_rows)) {
    reason <- c(reason, "insufficient_positive_rows")
  }
  if (sum(d$raw_count) < as.numeric(prep$config$gene_eligibility$minimum_total_raw_counts)) {
    reason <- c(reason, "insufficient_total_raw_counts")
  }
  if (length(reason)) {
    for (m in c("M0", "M1")) {
      ri <- ri + 1L
      results[[ri]] <- empty(
        ct, m, "not_tested", paste(reason, collapse = ";"), nrow(d), nlevels(d$absolute_tile),
        nlevels(d$Sample)
      )
    }
    next
  }
  if (!all(c("density_within", "density_between", "valid_tissue_area_mm2", "density_cells_per_mm2") %in%
    names(d))) {
    stop("Production M1 density fields missing")
  }
  sl <- unique(d[c("Sample", "AP_centered_mm", "density_between")])
  pr <- suppressWarnings(cor(sl$AP_centered_mm, sl$density_between, method = "pearson"))
  sr <- suppressWarnings(cor(sl$AP_centered_mm, sl$density_between, method = "spearman"))
  r2 <- tryCatch(summary(lm(AP_centered_mm ~ density_between, data = sl))$r.squared, error = function(e) NA_real_)
  vif <- if (is.finite(r2) && r2 < 1) {
    1 / (1 - r2)
  } else {
    Inf
  }
  cond <- tryCatch(axyd_design_condition(d), error = function(e) Inf)
  flag <- if (vif >= as.numeric(prep$config$collinearity$VIF_severe)) {
    "severe_collinearity"
  } else if (vif >= as.numeric(prep$config$collinearity$VIF_warning)) {
    "collinearity_warning"
  } else if (cond > as.numeric(prep$config$collinearity$condition_warning)) {
    "condition_warning"
  } else {
    "none"
  }
  for (m in c("M0", "M1")) {
    f <- axyd_cr2_fit(d, m)
    if (f$status != "tested") {
      ri <- ri + 1L
      results[[ri]] <- empty(
        ct, m, "model_failed", f$reason, nrow(d), nlevels(d$absolute_tile),
        nlevels(d$Sample)
      )
      next
    }
    ri <- ri + 1L
    results[[ri]] <- data.frame(
      planned_test_index = (j - 1L) * plan$n_celltypes + match(ct, plan$celltypes),
      Gene = gene, CellType = ct, model_id = m, grid_origin_id = "primary",
      support_design = prep$config$support$design,
      minimum_slices_per_absolute_tile = as.integer(prep$config$support$minimum_slices_per_absolute_tile),
      n_slices = nlevels(d$Sample),
      n_tiles = nlevels(d$absolute_tile), n_rows = nrow(d), AP_min_mm = min(d$AP_um) / 1000, AP_max_mm = max(d$AP_um) / 1000,
      AP_span_mm = diff(range(d$AP_um)) / 1000, beta_AP_per_mm = f$beta, SE_AP = f$se, CI_AP_lower = f$lo,
      CI_AP_upper = f$hi, test_statistic_AP = f$stat, df_AP = f$df, p_AP_nominal = f$p, AP_rate_ratio_per_mm = exp(f$beta),
      AP_percent_rate_change_per_mm = 100 * (exp(f$beta) - 1), covariance_method = "CR2 clustered by physical slice; Satterthwaite df",
      inference_status = "success", support_status = "tested", model_status = "tested", failure_reason = "",
      beta_density_within = if (m == "M1") {
        f$dw["beta"]
      } else {
        NA
      }, SE_density_within = if (m == "M1") {
        f$dw["se"]
      } else {
        NA
      }, beta_density_between = if (m == "M1") {
        f$db["beta"]
      } else {
        NA
      }, SE_density_between = if (m == "M1") {
        f$db["se"]
      } else {
        NA
      }, pearson_AP_density_between = if (m == "M1") {
        pr
      } else {
        NA
      }, spearman_AP_density_between = if (m == "M1") {
        sr
      } else {
        NA
      }, VIF_AP = if (m == "M1") {
        vif
      } else {
        NA
      }, condition_number = if (m == "M1") {
        cond
      } else {
        NA
      }, collinearity_flag = if (m == "M1") {
        flag
      } else {
        "not_applicable"
      }, AP_density_between_coefficient_correlation = if (m == "M1") {
        f$coef_cor
      } else {
        NA
      }, stringsAsFactors = FALSE
    )
  }
}

obj <- list(
  identity = plan$identity, status = "complete", Gene = gene, results = axyd_bind(results),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), peak_recorded_R_memory_mb = round(max(gc()[
    ,
    6
  ]), 1), session = capture.output(sessionInfo())
)

if (nrow(obj$results) != 40 || anyDuplicated(obj$results[c("Gene", "CellType", "model_id")])) {
  stop(
    "Invalid production checkpoint universe for ",
    gene
  )
}

axyd_rds(obj, checkpoint)

cat("Production M0/M1 complete for", gene, "\n")
