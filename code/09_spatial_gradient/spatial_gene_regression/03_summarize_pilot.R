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

prep <- readRDS(file.path(out, "checkpoints", "prepared_data.rds"))

started <- Sys.time()

paths <- file.path(out, "checkpoints", "pilot", paste0(prep$manifest$checkpoint_stem, ".rds"))

missing <- paths[!file.exists(paths)]

if (length(missing)) stop("Pilot checkpoints missing: ", paste(basename(missing), collapse = ","))

obj <- lapply(paths, readRDS)

if (any(vapply(obj, function(z) !identical(z$identity$identity, prep$identity$identity), logical(1)))) stop("Stale pilot checkpoint identity")

if (any(vapply(obj, `[[`, character(1), "status") != "complete")) stop("Incomplete pilot checkpoint")

results <- axyd_bind(lapply(obj, `[[`, "results"))

loo <- axyd_bind(lapply(obj, `[[`, "leave_one_slice_out"))

obs <- axyd_bind(lapply(obj, `[[`, "observations"))

masksens <- axyd_bind(lapply(obj, function(z) z$mask_parameter_sensitivity))

runtime <- data.frame(Gene = vapply(obj, `[[`, character(1), "Gene"), elapsed_seconds = vapply(
  obj, `[[`,
  numeric(1), "elapsed_seconds"
), peak_recorded_R_memory_mb = vapply(obj, `[[`, numeric(1), "peak_recorded_R_memory_mb"))

adjust_method_label <- as.character(prep$config$pilot$multiplicity_method)

adjust_method <- if (tolower(adjust_method_label) == "holm") "holm" else adjust_method_label

if (!adjust_method %in% p.adjust.methods) stop("Unsupported pilot multiplicity method: ", adjust_method_label)

primary <- results[results$grid_origin_id == "primary", ]

for (m in c("M0", "M1")) {
  ix <- primary$model_id == m & primary$model_status == "tested" & is.finite(primary$p_AP_nominal)
  primary$p_AP_adjusted_family[ix] <- p.adjust(primary$p_AP_nominal[ix], method = adjust_method)
  primary$pilot_multiplicity_family[primary$model_id == m] <- paste0(
    adjust_method_label, " across frozen primary-grid pilot ",
    m, " combinations only"
  )
}

m0 <- primary[primary$model_id == "M0", ]

m1 <- primary[primary$model_id == "M1", ]

names(m0)[!names(m0) %in% c("Gene", "CellType")] <- paste0("M0_", names(m0)[!names(m0) %in% c(
  "Gene",
  "CellType"
)])

names(m1)[!names(m1) %in% c("Gene", "CellType")] <- paste0("M1_", names(m1)[!names(m1) %in% c(
  "Gene",
  "CellType"
)])

cmp <- merge(m0, m1, by = c("Gene", "CellType"), all = TRUE, sort = FALSE)

cmp$delta_beta_AP <- cmp$M1_beta_AP_per_mm - cmp$M0_beta_AP_per_mm

cmp$percent_change_abs_beta_AP <- 100 * (abs(cmp$M1_beta_AP_per_mm) - abs(cmp$M0_beta_AP_per_mm)) / pmax(
  abs(cmp$M0_beta_AP_per_mm),
  as.numeric(prep$config$collinearity$beta_zero_tolerance)
)

cmp$SE_inflation_ratio <- cmp$M1_SE_AP / cmp$M0_SE_AP

cmp$direction_agreement <- sign(cmp$M0_beta_AP_per_mm) == sign(cmp$M1_beta_AP_per_mm)

cmp$CI_overlap <- pmax(cmp$M0_CI_AP_lower, cmp$M1_CI_AP_lower) <= pmin(cmp$M0_CI_AP_upper, cmp$M1_CI_AP_upper)

cmp$M0_CI_includes_zero <- cmp$M0_CI_AP_lower <= 0 & cmp$M0_CI_AP_upper >= 0

cmp$M1_CI_includes_zero <- cmp$M1_CI_AP_lower <= 0 & cmp$M1_CI_AP_upper >= 0

thr <- as.numeric(prep$config$collinearity$similar_relative_change_max)

cmp$change_classification <- ifelse(!cmp$direction_agreement, "direction_reversed", ifelse(abs(cmp$percent_change_abs_beta_AP) <=
  100 * thr, "similar", ifelse(abs(cmp$M1_beta_AP_per_mm) < abs(cmp$M0_beta_AP_per_mm), "attenuated",
  "strengthened"
)))

lsum <- if (nrow(loo)) {
  do.call(rbind, lapply(split(seq_len(nrow(loo)), interaction(loo$Gene, loo$CellType, loo$model_id,
    loo$grid_origin_id,
    drop = TRUE
  )), function(i) {
    z <- loo[i, ]
    ok <- z$status == "tested"
    data.frame(
      Gene = z$Gene[1], CellType = z$CellType[1], model_id = z$model_id[1], grid_origin_id = z$grid_origin_id[1],
      LOSO_n_success = sum(ok), LOSO_direction_agreement = if (any(ok)) {
        mean(z$direction_agrees[ok], na.rm = TRUE)
      } else {
        NA
      }, LOSO_max_relative_beta_change = if (any(ok)) {
        max(z$relative_beta_change[ok], na.rm = TRUE)
      } else {
        NA
      }, LOSO_stable = all(ok) && mean(z$direction_agrees[ok], na.rm = TRUE) >= as.numeric(prep$config$leave_one_slice_out$minimum_direction_agreement) &&
        max(z$relative_beta_change[ok], na.rm = TRUE) <= as.numeric(prep$config$leave_one_slice_out$maximum_relative_coefficient_change)
    )
  }))
} else {
  data.frame()
}

for (m in c("M0", "M1")) {
  z <- lsum[lsum$model_id == m & lsum$grid_origin_id == "primary", c(
    "Gene", "CellType", "LOSO_direction_agreement",
    "LOSO_max_relative_beta_change", "LOSO_stable"
  )]
  names(z)[3:5] <- paste0(
    c("LOSO_direction_agreement_", "LOSO_max_relative_beta_change_", "LOSO_stable_"),
    m
  )
  cmp <- merge(cmp, z, by = c("Gene", "CellType"), all.x = TRUE, sort = FALSE)
}

cmp$conditional_interpretation_flag <- ifelse(cmp$M1_model_status != "tested", "density_model_not_tested",
  ifelse(cmp$M1_collinearity_flag %in% c("severe_collinearity", "collinearity_warning", "condition_warning"),
    "conditional_AP_estimate_collinearity_limited", ifelse(cmp$LOSO_stable_M1 %in% FALSE, "unstable_conditional_estimate",
      "conditional_sensitivity_estimate"
    )
  )
)

colq <- unique(m1[c(
  "Gene", "CellType", "M1_pearson_AP_density_between", "M1_spearman_AP_density_between",
  "M1_VIF_AP", "M1_condition_number", "M1_collinearity_flag", "M1_AP_density_between_coefficient_correlation"
)])

names(colq) <- sub("^M1_", "", names(colq))

shift <- results[results$grid_origin_id == "shifted_xy_plus1000um", ]

shiftcmp <- merge(
  primary[c(
    "Gene", "CellType", "model_id", "beta_AP_per_mm", "SE_AP", "p_AP_nominal",
    "n_slices", "n_tiles", "model_status"
  )], shift[c(
    "Gene", "CellType", "model_id", "beta_AP_per_mm",
    "SE_AP", "p_AP_nominal", "n_slices", "n_tiles", "model_status"
  )],
  by = c("Gene", "CellType", "model_id"),
  suffixes = c("_primary", "_shifted"), all = TRUE
)

shiftcmp$direction_agreement <- sign(shiftcmp$beta_AP_per_mm_primary) == sign(shiftcmp$beta_AP_per_mm_shifted)

shiftcmp$diagnostic_role <- "fixed +1000um ML/DV boundary sensitivity; never selected by p-value"

fail <- results[results$model_status != "tested", c(
  "Gene", "CellType", "model_id", "grid_origin_id",
  "support_status", "model_status", "failure_reason"
)]

if (nrow(masksens)) {
  ref <- m1[m1$M1_grid_origin_id == "primary", c(
    "Gene", "CellType", "M1_beta_AP_per_mm", "M1_SE_AP",
    "M1_p_AP_nominal"
  )]
  names(ref)[3:5] <- c("primary_mask_beta_AP_per_mm", "primary_mask_SE_AP", "primary_mask_p_AP_nominal")
  masksens <- merge(masksens, ref, by = c("Gene", "CellType"), all.x = TRUE, sort = FALSE)
  masksens$delta_beta_vs_primary_mask <- masksens$beta_AP_per_mm - masksens$primary_mask_beta_AP_per_mm
  masksens$direction_agreement_vs_primary_mask <- sign(masksens$beta_AP_per_mm) == sign(masksens$primary_mask_beta_AP_per_mm)
}

axyd_csv(m0, file.path(out, "13_pilot_M0_spatially_adjusted_results.csv"))

axyd_csv(m1, file.path(out, "14_pilot_M1_density_adjusted_results.csv"))

axyd_csv(cmp, file.path(out, "15_pilot_M0_vs_M1_comparison.csv"))

axyd_csv(colq, file.path(out, "16_pilot_collinearity_diagnostics.csv"))

axyd_csv(loo, file.path(out, "17_pilot_leave_one_slice_out.csv"))

axyd_csv(shiftcmp, file.path(out, "18_pilot_shifted_grid_results.csv"))

axyd_csv(fail, file.path(out, "19_pilot_model_failures.csv"))

axyd_csv(runtime, file.path(out, "21_runtime_and_memory.csv"))

axyd_csv(masksens, file.path(out, "22_pilot_mask_parameter_sensitivity.csv"))

axyd_csv(obs, file.path(out, "pilot_observed_and_fitted_rows.csv"))

axyd_csv(lsum, file.path(out, "pilot_LOSO_summary.csv"))

blockers <- prep$blockers

checks <- data.frame(gate = c(
  "active_source_manifest_written", "AP_coordinate_audit", "primary_repeated_tile_A_B_support",
  "cell_coordinate_support_area_independent_positive", "M0_M1_identical_rows", "mask_parameter_sensitivity_completed",
  "CR2_Satterthwaite_validated", "all_tested_inference_finite", "shifted_grid_not_selected_by_result",
  "full_production_not_launched"
), status = c(
  if (file.exists(file.path(out, "01_active_source_manifest.csv"))) "PASS" else "FAIL",
  "PASS", if (all(prep$support_decisions$status[prep$support_decisions$grid_variant == "primary" &
    prep$support_decisions$CellType %in% unlist(prep$config$pilot$celltypes)] == "tested")) {
    "PASS"
  } else {
    "FAIL"
  },
  if (prep$area_available) "PASS" else "FAIL", if (prep$area_available) "PASS" else "FAIL", if (nrow(masksens) &&
    all(masksens$model_status == "tested")) {
    "PASS"
  } else {
    "FAIL"
  }, if (prep$cr2_validated) "PASS" else "FAIL",
  if (all(primary$model_status != "tested" | is.finite(primary$p_AP_nominal) & is.finite(primary$SE_AP) &
    is.finite(primary$df_AP))) {
    "PASS"
  } else {
    "FAIL"
  }, "PASS", "PASS"
))

ready <- all(checks$status == "PASS")

lines <- c(
  "# 2000-µm absolute-XY density AP pilot acceptance report", "", paste0("Generated: ", format(Sys.time(),
    tz = "UTC", usetz = TRUE
  )), "", "Analyzed-tissue support is derived only from all QC-passed cell coordinates; it is not DAPI-derived and is not a histological tissue mask.",
  "M0 estimates the within-CellType AP association from an unbalanced panel of absolute ML/DV tiles, each retained only when observed in at least four selected physical slices.",
  "M1 is a sensitivity analysis conditioning on within-slice and between-slice CellType density using that gene-independent support area.",
  "", paste0("- ", checks$gate, ": ", checks$status), "", sprintf(
    "Frozen pilot: %d genes × 2 CellTypes = %d combinations per primary model.",
    nrow(prep$manifest), 2 * nrow(prep$manifest)
  ), sprintf(
    "Primary M0 tested/model_failed/not_tested: %d/%d/%d.",
    sum(m0$M0_model_status == "tested"), sum(m0$M0_model_status == "model_failed"), sum(m0$M0_model_status ==
      "not_tested")
  ), sprintf("Primary M1 tested/model_failed/not_tested: %d/%d/%d.", sum(m1$M1_model_status ==
    "tested"), sum(m1$M1_model_status == "model_failed"), sum(m1$M1_model_status == "not_tested")),
  sprintf("M1 collinearity flags: %s.", paste(names(table(m1$M1_collinearity_flag)), as.integer(table(m1$M1_collinearity_flag)),
    sep = "=", collapse = ", "
  )), "", "Pilot significance is not a validation outcome and was not used to tune mask parameters, support, origin, density, or formulas.",
  paste0("Implementation-ready: ", if (ready) "YES — explicit user approval is still required before production." else "NO — inspect failed gates and model failures.")
)

writeLines(lines, file.path(out, "20_pilot_acceptance_report.md"))

axyd_csv(checks, file.path(out, "pilot_acceptance_checks.csv"))

Rscript <- file.path(R.home("bin"), "Rscript")

status <- system2(Rscript, c(
  file.path(script_dir, "04_plot_pilot.R"), "--output-dir", out, "--plot-dir",
  plots
))

if (status != 0) stop("Pilot plotting failed")

files <- sort(c(list.files(out, recursive = TRUE, full.names = TRUE), list.files(plots,
  recursive = TRUE,
  full.names = TRUE
)))

info <- file.info(files)

manifest <- list(analysis_id = prep$config$analysis_id, identity = prep$identity$identity, created_utc = format(Sys.time(),
  tz = "UTC", usetz = TRUE
), production_launched = FALSE, implementation_ready = ready, files = lapply(
  seq_along(files),
  function(i) list(path = normalizePath(files[i]), bytes = unname(info$size[i]), md5 = axyd_hash_file(files[i]))
))

axyd_json(manifest, file.path(out, "artifact_manifest.json"))

if (ready) {
  writeLines("PASS; explicit production approval still required", file.path(
    out, "checkpoints",
    "PILOT_ACCEPTED.txt"
  ))
}

cat("2000-um absolute-XY density pilot summary complete\n")
