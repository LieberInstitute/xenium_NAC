#!/usr/bin/env Rscript
argv <- commandArgs(FALSE)

script_dir <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1])))

source(file.path(script_dir, "common.R"))

args <- commandArgs(TRUE)

val <- function(n, default = NULL) {
  i <- match(n, args)
  if (is.na(i)) {
    return(default)
  }
  if (i == length(args)) {
    stop("Missing value for ", n)
  }
  args[i + 1]
}

out <- val("--output-dir")

gene <- val("--gene")

resume <- "--resume" %in% args

prep <- readRDS(file.path(out, "checkpoints", "prepared_data.rds"))

if (!file.exists(file.path(out, "checkpoints", "PREPARE_M0_READY.txt")) || !isTRUE(prep$cr2_validated)) stop("M0 preparation/CR2 hard gates have not passed")

if (!gene %in% prep$manifest$Gene) stop("Unknown pilot gene: ", gene)

j <- match(gene, prep$manifest$Gene)

checkpoint <- file.path(out, "checkpoints", "pilot", paste0(prep$manifest$checkpoint_stem[j], ".rds"))

dir.create(dirname(checkpoint), FALSE, TRUE)

if (resume && file.exists(checkpoint)) {
  q <- readRDS(checkpoint)
  if (identical(q$identity$identity, prep$identity$identity)) {
    cat("Reusing pilot checkpoint", gene, "\n")
    quit(save = "no")
  }
}

started <- Sys.time()

counts <- list(primary = read.csv(gzfile(file.path(prep$input_dir, "primary_counts.csv.gz")), check.names = FALSE)[c(
  "tile_group_id",
  gene
)], shifted_xy_plus1000um = read.csv(gzfile(file.path(prep$input_dir, "shifted_xy_plus1000um_counts.csv.gz")),
  check.names = FALSE
)[c("tile_group_id", gene)])

panel_labels <- prep$panel$source_label[prep$panel$Gene == gene]

label <- paste(panel_labels, collapse = ";")

empty_result <- function(ct, v, m, status, reason, n = 0, nt = 0, ns = 0) {
  data.frame(
    Gene = gene, CellType = ct,
    marker_source_label = label, marker_own_or_cross_label = if (ct %in% panel_labels) "own_label" else "cross_label",
    model_id = m, grid_origin_id = v, support_design = prep$config$support$design,
    minimum_slices_per_absolute_tile = as.integer(prep$config$support$minimum_slices_per_absolute_tile),
    n_slices = ns, n_tiles = nt, n_rows = n, AP_min_mm = NA, AP_max_mm = NA,
    AP_span_mm = NA, beta_AP_per_mm = NA, SE_AP = NA, CI_AP_lower = NA, CI_AP_upper = NA, test_statistic_AP = NA,
    df_AP = NA, p_AP_nominal = NA, p_AP_adjusted_family = NA, AP_rate_ratio_per_mm = NA, AP_percent_rate_change_per_mm = NA,
    covariance_method = "CR2 clustered by physical slice", inference_status = status, support_status = status,
    model_status = status, failure_reason = reason, beta_density_within = NA, SE_density_within = NA,
    beta_density_between = NA, SE_density_between = NA, pearson_AP_density_between = NA, spearman_AP_density_between = NA,
    VIF_AP = NA, condition_number = NA, collinearity_flag = "", AP_density_between_coefficient_correlation = NA,
    stringsAsFactors = FALSE
  )
}

results <- list()

loso <- list()

obsout <- list()

masksens <- list()

ri <- li <- oi <- mi <- 0L

for (v in c("primary", "shifted_xy_plus1000um")) {
  for (ct in unlist(prep$config$pilot$celltypes)) {
    d <- prep$rows[prep$rows$grid_variant == v & prep$rows$CellType == ct, ]
    decision <- prep$support_decisions[prep$support_decisions$grid_variant == v & prep$support_decisions$CellType ==
      ct, ]
    if (!nrow(decision) || decision$status != "tested" || !nrow(d)) {
      for (m in c("M0", "M1")) {
        ri <- ri + 1L
        results[[ri]] <- empty_result(ct, v, m, "not_tested", if (nrow(decision)) {
          decision$failure_reason
        } else {
          "support_decision_missing"
        })
      }
      next
    }
    cc <- counts[[v]]
    d$raw_count <- cc[[gene]][match(d$tile_group_id, cc$tile_group_id)]
    d <- d[order(d$AP_um, d$absolute_tile), ]
    d$absolute_tile <- droplevels(factor(d$absolute_tile))
    d$Sample <- factor(d$Sample, levels = unique(d$Sample[order(d$AP_um)]))
    if (anyNA(d$raw_count) || any(d$raw_count < 0) || any(d$exposure <= 0) || any(d$n_cells <= 0)) {
      stop("Count/exposure row reconciliation failed for ", gene, " / ", ct, " / ", v)
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
        results[[ri]] <- empty_result(
          ct, v, m, "not_tested", paste(reason, collapse = ";"),
          nrow(d), nlevels(d$absolute_tile), nlevels(d$Sample)
        )
      }
      next
    }
    has_density <- isTRUE(prep$area_available) && all(c(
      "density_within", "density_between", "valid_tissue_area_mm2",
      "density_cells_per_mm2"
    ) %in% names(d))
    if (has_density) {
      sl <- unique(d[c("Sample", "AP_um", "AP_centered_mm", "density_between")])
      pr <- cor(sl$AP_centered_mm, sl$density_between, method = "pearson")
      sr <- cor(sl$AP_centered_mm, sl$density_between, method = "spearman")
      r2 <- summary(lm(AP_centered_mm ~ density_between, data = sl))$r.squared
      vif <- if (is.finite(r2) && r2 < 1) {
        1 / (1 - r2)
      } else {
        Inf
      }
      cond <- axyd_design_condition(d)
      flag <- if (vif >= as.numeric(prep$config$collinearity$VIF_severe)) {
        "severe_collinearity"
      } else if (vif >= as.numeric(prep$config$collinearity$VIF_warning)) {
        "collinearity_warning"
      } else if (cond > as.numeric(prep$config$collinearity$condition_warning)) {
        "condition_warning"
      } else {
        "none"
      }
    } else {
      pr <- sr <- vif <- cond <- NA_real_
      flag <- "density_area_unavailable"
    }
    for (m in c("M0", "M1")) {
      if (m == "M1" && !has_density) {
        ri <- ri + 1L
        q <- empty_result(
          ct, v, m, "not_tested", "cell_coordinate_support_area_unavailable",
          nrow(d), nlevels(d$absolute_tile), nlevels(d$Sample)
        )
        q$support_status <- "tested"
        results[[ri]] <- q
        next
      }
      f <- axyd_cr2_fit(d, m)
      if (f$status != "tested") {
        ri <- ri + 1L
        results[[ri]] <- empty_result(
          ct, v, m, "model_failed", f$reason, nrow(d), nlevels(d$absolute_tile),
          nlevels(d$Sample)
        )
        next
      }
      ri <- ri + 1L
      results[[ri]] <- data.frame(
        Gene = gene, CellType = ct, marker_source_label = label, marker_own_or_cross_label = if (ct %in%
          panel_labels) {
          "own_label"
        } else {
          "cross_label"
        }, model_id = m, grid_origin_id = v, support_design = prep$config$support$design,
        minimum_slices_per_absolute_tile = as.integer(prep$config$support$minimum_slices_per_absolute_tile),
        n_slices = nlevels(d$Sample), n_tiles = nlevels(d$absolute_tile),
        n_rows = nrow(d), AP_min_mm = min(d$AP_um) / 1000, AP_max_mm = max(d$AP_um) / 1000, AP_span_mm = diff(range(d$AP_um)) / 1000,
        beta_AP_per_mm = f$beta, SE_AP = f$se, CI_AP_lower = f$lo, CI_AP_upper = f$hi, test_statistic_AP = f$stat,
        df_AP = f$df, p_AP_nominal = f$p, p_AP_adjusted_family = NA, AP_rate_ratio_per_mm = exp(f$beta),
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
      keep <- intersect(c(
        "Sample", "AP_um", "CellType", "grid_variant", "absolute_tile", "bx",
        "by", "n_cells", "exposure", "valid_tissue_area_mm2", "density_cells_per_mm2", "density_within",
        "density_between", "raw_count"
      ), names(d))
      oo <- d[keep]
      oo$Gene <- gene
      oo$model_id <- m
      oo$fitted_count <- fitted(f$fit)
      oo$fitted_rate <- oo$fitted_count / oo$exposure
      oo$observed_rate <- oo$raw_count / oo$exposure
      oi <- oi + 1L
      obsout[[oi]] <- oo
      for (s in levels(d$Sample)) {
        dd <- droplevels(d[d$Sample != s, ])
        lf <- axyd_cr2_fit(dd, m)
        li <- li + 1L
        loso[[li]] <- data.frame(
          Gene = gene, CellType = ct, model_id = m, grid_origin_id = v,
          omitted_Sample = s, n_slices_remaining = nlevels(dd$Sample), full_beta_AP_per_mm = f$beta,
          LOSO_beta_AP_per_mm = if (lf$status == "tested") {
            lf$beta
          } else {
            NA
          }, LOSO_SE_AP = if (lf$status == "tested") {
            lf$se
          } else {
            NA
          }, LOSO_df_AP = if (lf$status == "tested") {
            lf$df
          } else {
            NA
          }, LOSO_p_AP = if (lf$status == "tested") {
            lf$p
          } else {
            NA
          }, direction_agrees = if (lf$status == "tested") {
            sign(lf$beta) == sign(f$beta)
          } else {
            NA
          }, relative_beta_change = if (lf$status == "tested") {
            abs(lf$beta - f$beta) / max(abs(f$beta), as.numeric(prep$config$collinearity$beta_zero_tolerance))
          } else {
            NA
          }, status = lf$status, failure_reason = lf$reason, stringsAsFactors = FALSE
        )
      }
    }
    if (has_density && nrow(prep$mask_sensitivity_rows)) {
      for (ps in unique(prep$mask_sensitivity_rows$mask_parameter_set)) {
        ds <- prep$mask_sensitivity_rows[prep$mask_sensitivity_rows$grid_variant == v & prep$mask_sensitivity_rows$CellType ==
          ct & prep$mask_sensitivity_rows$mask_parameter_set == ps, ]
        key <- paste(d$Sample, d$tile_group_id)
        skey <- paste(ds$Sample, ds$tile_group_id)
        ds <- ds[match(key, skey), ]
        if (nrow(ds) != nrow(d) || anyNA(ds$tile_group_id) || !identical(
          as.character(ds$Sample),
          as.character(d$Sample)
        ) || !identical(as.character(ds$absolute_tile), as.character(d$absolute_tile))) {
          stop(
            "M0/M1 mask-sensitivity row contract failed for ", gene, " / ", ct, " / ", v,
            " / ", ps
          )
        }
        ds$raw_count <- d$raw_count
        ds$Sample <- d$Sample
        ds$absolute_tile <- d$absolute_tile
        sf <- axyd_cr2_fit(ds, "M1")
        mi <- mi + 1L
        masksens[[mi]] <- data.frame(
          Gene = gene, CellType = ct, grid_origin_id = v, mask_parameter_set = ps,
          n_slices = nlevels(ds$Sample), n_tiles = nlevels(ds$absolute_tile), n_rows = nrow(ds),
          model_status = sf$status, failure_reason = if (sf$status == "tested") {
            ""
          } else {
            sf$reason
          }, beta_AP_per_mm = if (sf$status == "tested") {
            sf$beta
          } else {
            NA
          }, SE_AP = if (sf$status == "tested") {
            sf$se
          } else {
            NA
          }, CI_AP_lower = if (sf$status == "tested") {
            sf$lo
          } else {
            NA
          }, CI_AP_upper = if (sf$status == "tested") {
            sf$hi
          } else {
            NA
          }, df_AP = if (sf$status == "tested") {
            sf$df
          } else {
            NA
          }, p_AP_nominal = if (sf$status == "tested") {
            sf$p
          } else {
            NA
          }, AP_rate_ratio_per_mm = if (sf$status == "tested") {
            exp(sf$beta)
          } else {
            NA
          }, stringsAsFactors = FALSE
        )
      }
    }
  }
}

obj <- list(
  identity = prep$identity, status = "complete", Gene = gene, results = axyd_bind(results),
  leave_one_slice_out = axyd_bind(loso), observations = axyd_bind(obsout), mask_parameter_sensitivity = axyd_bind(masksens),
  elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")), peak_recorded_R_memory_mb = round(max(gc()[
    ,
    6
  ]), 1), session = capture.output(sessionInfo())
)

axyd_rds(obj, checkpoint)

cat("2000-um absolute-XY M0/M1 pilot plus mask-parameter sensitivity complete for", gene, "\n")
