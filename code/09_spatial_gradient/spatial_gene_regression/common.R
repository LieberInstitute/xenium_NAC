suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
})

axyd_hash_file <- function(p) unname(tools::md5sum(p))

axyd_hash_text <- function(x) {
  p <- tempfile()
  on.exit(unlink(p))
  writeLines(x, p)
  axyd_hash_file(p)
}

axyd_bool <- function(x) tolower(as.character(x)) %in% c("true", "t", "1")

axyd_bind <- function(x) {
  x <- x[!vapply(x, function(z) is.null(z) || !nrow(z), logical(1))]
  if (!length(x)) {
    return(data.frame())
  }
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(z) {
    missing <- setdiff(cols, names(z))
    if (length(missing)) {
      z[missing] <- NA
    }
    z[cols]
  })
  do.call(rbind, x)
}

axyd_csv <- function(x, p) {
  dir.create(dirname(p), FALSE, TRUE)
  q <- paste0(p, ".tmp.", Sys.getpid())
  write.csv(x, q, row.names = FALSE, na = "")
  if (!file.rename(q, p)) {
    stop("Atomic CSV rename failed: ", p)
  }
}

axyd_rds <- function(x, p) {
  dir.create(dirname(p), FALSE, TRUE)
  q <- paste0(p, ".tmp.", Sys.getpid())
  saveRDS(x, q)
  if (!file.rename(q, p)) {
    stop("Atomic RDS rename failed: ", p)
  }
}

axyd_json <- function(x, p) {
  dir.create(dirname(p), FALSE, TRUE)
  q <- paste0(p, ".tmp.", Sys.getpid())
  writeLines(toJSON(x, pretty = TRUE, auto_unbox = TRUE, null = "null", na = "null"), q)
  if (!file.rename(q, p)) {
    stop("Atomic JSON rename failed: ", p)
  }
}

axyd_validate_config <- function(cfg) {
  stopifnot(cfg$donor == "Br6660", as.numeric(cfg$grid$width_um) == 2000, as.numeric(cfg$grid$height_um) ==
    2000, identical(as.numeric(unlist(cfg$grid$primary_shift_um)), c(0, 0)), identical(
    as.numeric(unlist(cfg$grid$shifted_sensitivity_um)),
    c(1000, 1000)
  ))
  if (cfg$model$covariance != "CR2_clustered_by_physical_slice" || cfg$model$df_method != "Satterthwaite" ||
    isTRUE(cfg$model$AP_density_interaction)) {
    stop("Frozen inference contract changed")
  }
  support <- cfg$support
  if (
    support$design != "unbalanced_repeated_absolute_tiles" ||
      as.integer(support$minimum_slices_per_absolute_tile) != 4L ||
      as.integer(support$minimum_contiguous_slices) != 8L ||
      as.numeric(support$minimum_ap_span_fraction) != 0.60
  ) {
    stop("Frozen unbalanced repeated-absolute-tile support contract changed")
  }
  ta <- cfg$tissue_area
  if (ta$source != "cell-coordinate-derived_analyzed-tissue_support_mask" || !isTRUE(ta$use_all_qc_passed_cells_regardless_of_CellType) ||
    ta$mask_method != "rasterized_buffered_cell_union" || as.numeric(ta$raster_resolution_um) !=
    25 || as.numeric(ta$primary_cell_buffer_um) != 100 || !identical(
    as.numeric(unlist(ta$sensitivity_cell_buffer_um)),
    c(75, 125)
  )) {
    stop("Frozen cell-coordinate support-mask contract changed")
  }
  ap <- as.numeric(unlist(cfg$samples))
  if (length(ap) != 11 || any(diff(ap) <= 0)) {
    stop("Frozen biological AP order invalid")
  }
}

axyd_identity <- function(cfg_path, input_dir, marker_path, sources) {
  man <- file.path(input_dir, "INPUT_CACHE_MANIFEST.json")
  if (!file.exists(man)) {
    stop("Input manifest missing")
  }
  maskman <- file.path(input_dir, "CELL_SUPPORT_MASK_MANIFEST.json")
  if (!file.exists(maskman)) {
    stop("Cell-support-mask manifest missing")
  }
  parts <- c(
    axyd_hash_file(cfg_path), axyd_hash_file(man), axyd_hash_file(maskman), axyd_hash_file(marker_path),
    vapply(sources, axyd_hash_file, character(1))
  )
  list(identity = axyd_hash_text(parts), components = list(
    config_md5 = parts[1], input_manifest_md5 = parts[2],
    mask_manifest_md5 = parts[3], marker_panel_md5 = parts[4], source_md5 = unname(parts[-(1:4)])
  ))
}

axyd_candidate_intervals <- function(ok, min_n) {
  idx <- which(ok)
  if (!length(idx)) {
    return(list())
  }
  starts <- idx[c(TRUE, diff(idx) > 1)]
  ends <- idx[c(diff(idx) > 1, TRUE)]
  out <- list()
  k <- 0L
  for (j in seq_along(starts)) {
    for (a in starts[j]:ends[j]) {
      for (b in a:ends[j]) {
        if (b - a + 1 >= min_n) {
          k <- k + 1L
          out[[k]] <- a:b
        }
      }
    }
  }
  out
}

axyd_select_interval <- function(d, cfg) {
  d <- d[order(d$AP_um), ]
  ok <- d$n_cells >= as.integer(cfg$support$minimum_cells_per_slice_celltype) & d$exposure >= as.numeric(cfg$support$minimum_exposure_per_slice_celltype)
  cand <- axyd_candidate_intervals(ok, as.integer(cfg$support$minimum_contiguous_slices))
  full_span <- diff(range(as.numeric(unlist(cfg$samples))))
  rows <- lapply(seq_along(cand), function(i) {
    ix <- cand[[i]]
    data.frame(
      candidate_id = i, start_index = min(ix), end_index = max(ix), n_slices = length(ix),
      AP_min_um = min(d$AP_um[ix]), AP_max_um = max(d$AP_um[ix]), AP_span_um = diff(range(d$AP_um[ix])),
      AP_span_fraction = diff(range(d$AP_um[ix])) / full_span, total_exposure = sum(d$exposure[ix]),
      total_cells = sum(d$n_cells[ix]), passes_span = diff(range(d$AP_um[ix])) / full_span >= as.numeric(cfg$support$minimum_ap_span_fraction)
    )
  })
  tab <- axyd_bind(rows)
  if (!nrow(tab) || !any(tab$passes_span)) {
    return(list(
      status = "not_tested", reason = "no_contiguous_interval_meets_abundance_exposure_and_AP_span",
      samples = character(), candidates = tab
    ))
  }
  q <- tab[tab$passes_span, ]
  q <- q[order(-q$n_slices, -q$total_exposure, -q$total_cells, q$start_index), ]
  win <- cand[[q$candidate_id[1]]]
  list(status = "tested", reason = "", samples = as.character(d$Sample[win]), candidates = transform(tab,
    selected = candidate_id == q$candidate_id[1]
  ))
}

axyd_build_support <- function(meta, whole, cfg, variant) {
  cts <- sort(unique(as.character(whole$CellType)))
  dec <- list()
  cand <- list()
  recurrence_out <- list()
  ret <- list()
  rows <- list()
  di <- ci <- si <- ri <- oi <- 0L

  for (ct in cts) {
    w <- whole[whole$CellType == ct, ]
    sel <- axyd_select_interval(w, cfg)

    if (nrow(sel$candidates)) {
      sel$candidates$CellType <- ct
      ci <- ci + 1L
      cand[[ci]] <- sel$candidates
    }

    if (sel$status != "tested") {
      di <- di + 1L
      dec[[di]] <- data.frame(
        CellType = ct, grid_variant = variant, status = "not_tested", failure_reason = sel$reason,
        n_slices = 0, n_repeated_tiles = 0, repeated_tile_footprint_nominal_mm2 = 0
      )
      next
    }

    z <- meta[meta$CellType == ct & meta$Sample %in% sel$samples, ]
    z$qualifying <- z$n_cells >= as.integer(cfg$support$minimum_cells_per_slice_celltype_tile) &
      z$exposure > as.numeric(cfg$support$minimum_exposure_per_slice_celltype_tile)
    sets <- lapply(sel$samples, function(s) unique(z$absolute_tile[z$Sample == s & z$qualifying]))
    names(sets) <- sel$samples

    qualifying_rows <- z[z$qualifying, ]
    recurrence <- if (nrow(qualifying_rows)) {
      axyd_bind(lapply(split(qualifying_rows, qualifying_rows$absolute_tile), function(tile_rows) {
        first_row <- tile_rows[1, ]
        data.frame(
          CellType = ct,
          grid_variant = variant,
          absolute_tile = first_row$absolute_tile,
          bx = first_row$bx,
          by = first_row$by,
          ML_left_um = first_row$ML_left_um,
          ML_right_um = first_row$ML_right_um,
          DV_bottom_um = first_row$DV_bottom_um,
          DV_top_um = first_row$DV_top_um,
          n_selected_slices = length(sel$samples),
          n_slices_present = length(unique(tile_rows$Sample)),
          AP_min_um_present = min(tile_rows$AP_um),
          AP_max_um_present = max(tile_rows$AP_um),
          AP_span_um_present = diff(range(tile_rows$AP_um)),
          stringsAsFactors = FALSE
        )
      }))
    } else {
      data.frame()
    }

    if (nrow(recurrence)) {
      recurrence$minimum_slices_required <- as.integer(
        cfg$support$minimum_slices_per_absolute_tile
      )
      recurrence$retained_repeated_tile <- recurrence$n_slices_present >=
        recurrence$minimum_slices_required
      recurrence$recurrence_fraction <- recurrence$n_slices_present /
        recurrence$n_selected_slices
      retained_tiles <- recurrence$absolute_tile[recurrence$retained_repeated_tile]
    } else {
      retained_tiles <- character()
    }

    q <- qualifying_rows[qualifying_rows$absolute_tile %in% retained_tiles, ]
    before <- w[match(sel$samples, w$Sample), ]
    rc <- vapply(sel$samples, function(s) sum(q$n_cells[q$Sample == s]), numeric(1))
    re <- vapply(sel$samples, function(s) sum(q$exposure[q$Sample == s]), numeric(1))
    cf <- rc / before$n_cells
    ef <- re / before$exposure

    checks <- c(
      repeated_tiles = length(retained_tiles) >= as.integer(cfg$support$minimum_repeated_tiles),
      repeated_tile_footprint = length(retained_tiles) * 4 >=
        as.numeric(cfg$support$minimum_repeated_tile_footprint_mm2),
      cell_retention = all(cf >= as.numeric(
        cfg$support$minimum_retained_cell_fraction_per_slice
      )),
      exposure_retention = all(ef >= as.numeric(
        cfg$support$minimum_retained_exposure_fraction_per_slice
      ))
    )
    pass <- all(checks)
    reason <- if (pass) {
      ""
    } else {
      paste(names(checks)[!checks], collapse = ";")
    }
    di <- di + 1L
    dec[[di]] <- data.frame(
      CellType = ct, grid_variant = variant, status = if (pass) {
        "tested"
      } else {
        "not_tested"
      }, failure_reason = reason, n_slices = length(sel$samples), selected_samples = paste(sel$samples,
        collapse = ";"
      ), AP_min_um = min(before$AP_um), AP_max_um = max(before$AP_um),
      minimum_slices_per_absolute_tile = as.integer(cfg$support$minimum_slices_per_absolute_tile),
      n_repeated_tiles = length(retained_tiles),
      repeated_tile_footprint_nominal_mm2 = length(retained_tiles) * 4,
      min_retained_cell_fraction = min(cf), min_retained_exposure_fraction = min(ef)
    )

    if (nrow(recurrence)) {
      recurrence$support_status <- if (pass) "tested" else "not_tested"
      recurrence$support_failure_reason <- reason
      si <- si + 1L
      recurrence_out[[si]] <- recurrence
    }

    for (j in seq_along(sel$samples)) {
      ri <- ri + 1L
      ret[[ri]] <- data.frame(
        CellType = ct, grid_variant = variant, Sample = sel$samples[j], AP_um = before$AP_um[j],
        n_qualifying_tiles = length(sets[[j]]),
        n_retained_repeated_tiles = length(unique(q$absolute_tile[q$Sample == sel$samples[j]])),
        n_repeated_tiles_overall = length(retained_tiles), cells_before = before$n_cells[j],
        cells_retained = rc[j], retained_cell_fraction = cf[j], exposure_before = before$exposure[j],
        exposure_retained = re[j], retained_exposure_fraction = ef[j], support_status = if (pass) {
          "tested"
        } else {
          "not_tested"
        }, failure_reason = reason
      )
    }

    if (pass) {
      qq <- q[order(match(q$Sample, sel$samples), q$absolute_tile), ]
      for (j in seq_len(nrow(qq))) {
        oi <- oi + 1L
        rows[[oi]] <- qq[j, ]
      }
    }
  }

  list(
    decisions = axyd_bind(dec), candidates = axyd_bind(cand),
    tile_recurrence = axyd_bind(recurrence_out), retention = axyd_bind(ret),
    rows = axyd_bind(rows)
  )
}

axyd_validate_area_and_density <- function(rows, area_path, cfg, mask_parameter_set = as.character(cfg$tissue_area$primary_parameter_set)) {
  if (!nzchar(area_path) || !file.exists(area_path)) {
    return(list(available = FALSE, reason = "cell_coordinate_support_area_file_missing", area = data.frame(
      status = "unavailable",
      failure_reason = "cell_coordinate_support_area_file_missing"
    ), density = data.frame(
      status = "unavailable",
      failure_reason = "density_denominator_unavailable"
    ), decomposition = data.frame(
      status = "unavailable",
      failure_reason = "density_denominator_unavailable"
    ), rows = rows))
  }
  a <- read.csv(area_path, stringsAsFactors = FALSE)
  req <- c(
    "Sample", "grid_variant", "bx", "by", "valid_tissue_area_mm2", "area_source", "area_source_sha256",
    "independent_of_gene_and_celltype"
  )
  if (!all(req %in% names(a))) {
    stop("Tissue-area CSV lacks: ", paste(setdiff(req, names(a)), collapse = ","))
  }
  if (!"mask_parameter_set" %in% names(a)) {
    stop("Cell-coordinate support area lacks mask_parameter_set")
  }
  a <- a[a$mask_parameter_set == mask_parameter_set, ]
  if (!nrow(a)) {
    stop("Mask parameter set absent from area table: ", mask_parameter_set)
  }
  if (!all(axyd_bool(a$independent_of_gene_and_celltype)) || any(!is.finite(a$valid_tissue_area_mm2) |
    a$valid_tissue_area_mm2 <= as.numeric(cfg$tissue_area$minimum_valid_tissue_area_mm2) | a$valid_tissue_area_mm2 >
    4)) {
    stop("Invalid or non-independent valid tissue area")
  }
  if (any(a$area_source != "cell-coordinate-derived analyzed-tissue support mask")) {
    stop("Area source is not the frozen cell-coordinate-derived support mask")
  }
  if (anyDuplicated(a[c("Sample", "grid_variant", "bx", "by")])) {
    stop("Duplicated tissue-area key")
  }
  z <- merge(rows, a, by = c("Sample", "grid_variant", "bx", "by"), all.x = TRUE, sort = FALSE)
  bad <- !is.finite(z$valid_tissue_area_mm2) | z$valid_tissue_area_mm2 <= 0
  if (any(bad)) {
    stop("Missing/nonpositive valid tissue area for retained row: ", paste(z$Sample[which(bad)[1]],
      z$CellType[which(bad)[1]], z$absolute_tile[which(bad)[1]],
      sep = " / "
    ))
  }
  z$density_cells_per_mm2 <- z$n_cells / z$valid_tissue_area_mm2
  if (any(!is.finite(z$density_cells_per_mm2) | z$density_cells_per_mm2 <= 0)) {
    stop("Invalid retained density")
  }
  z$log_density <- log(z$density_cells_per_mm2)
  keys <- interaction(z$CellType, z$Sample, z$grid_variant, drop = TRUE)
  means <- tapply(seq_len(nrow(z)), keys, function(i) weighted.mean(z$log_density[i], z$valid_tissue_area_mm2[i]))
  z$density_slice_mean <- as.numeric(means[as.character(keys)])
  z$density_within <- z$log_density - z$density_slice_mean
  ctkey <- interaction(z$CellType, z$grid_variant, drop = TRUE)
  between_mean <- tapply(z$density_slice_mean, ctkey, mean)
  z$density_between <- z$density_slice_mean - as.numeric(between_mean[as.character(ctkey)])
  spatial <- c(
    "bx", "by", "ML_left_um", "ML_right_um", "DV_bottom_um", "DV_top_um", "ML_center_um",
    "DV_center_um"
  )
  maskmeta <- intersect(c(
    "mask_parameter_set", "mask_method", "raster_resolution_um", "cell_buffer_um",
    "fill_holes_max_mm2", "remove_components_min_mm2", "connectivity", "n_qc_passed_cells"
  ), names(z))
  area <- unique(z[c(
    "Sample", "grid_variant", "absolute_tile", spatial, "valid_tissue_area_mm2", maskmeta,
    "area_source", "area_source_sha256", "independent_of_gene_and_celltype"
  )])
  density <- z[c(
    "Sample", "AP_um", "CellType", "grid_variant", "absolute_tile", spatial, "n_cells",
    "valid_tissue_area_mm2", maskmeta, "density_cells_per_mm2", "log_density"
  )]
  decomp <- z[c(
    "Sample", "AP_um", "CellType", "grid_variant", "absolute_tile", spatial, "valid_tissue_area_mm2",
    maskmeta, "log_density", "density_slice_mean", "density_within", "density_between"
  )]
  list(available = TRUE, reason = "", area = area, density = density, decomposition = decomp, rows = z)
}

axyd_require_cr2 <- function() {
  if (!requireNamespace("clubSandwich", quietly = TRUE)) {
    stop("CR2_INFERENCE_BLOCKED: R package clubSandwich is not installed; no naive covariance fallback is permitted")
  }
  set.seed(1)
  d <- data.frame(y = rpois(40, 5), tile = factor(rep(1:5, 8)), AP_centered_mm = rep(seq(-1, 1, length.out = 8),
    each = 5
  ), Sample = factor(rep(1:8, each = 5)), exposure = 1)
  f <- glm(y ~ tile + AP_centered_mm,
    offset = log(exposure), family = quasipoisson(link = "log"),
    data = d
  )
  v <- try(clubSandwich::vcovCR(f, cluster = d$Sample, type = "CR2"), silent = TRUE)
  if (inherits(v, "try-error")) {
    stop("CR2_INFERENCE_BLOCKED: clubSandwich cannot construct CR2 for quasipoisson glm")
  }
  q <- try(clubSandwich::coef_test(f, vcov = v, test = "Satterthwaite"), silent = TRUE)
  if (inherits(q, "try-error") || !"AP_centered_mm" %in% rownames(q)) {
    stop("CR2_INFERENCE_BLOCKED: Satterthwaite coefficient test validation failed")
  }
  TRUE
}

axyd_cr2_fit <- function(d, model_id) {
  form <- if (model_id == "M0") {
    raw_count ~ absolute_tile + AP_centered_mm
  } else {
    raw_count ~ absolute_tile + AP_centered_mm + density_within + density_between
  }
  fit <- try(glm(form, offset = log(exposure), family = quasipoisson(link = "log"), data = d), silent = TRUE)
  if (inherits(fit, "try-error") || !isTRUE(fit$converged)) {
    return(list(status = "model_failed", reason = "glm_failed_or_nonconverged"))
  }
  V <- try(clubSandwich::vcovCR(fit, cluster = d$Sample, type = "CR2"), silent = TRUE)
  if (inherits(V, "try-error") || any(!is.finite(V))) {
    return(list(status = "model_failed", reason = "CR2_covariance_failed"))
  }
  ct <- try(clubSandwich::coef_test(fit, vcov = V, test = "Satterthwaite"), silent = TRUE)
  if (inherits(ct, "try-error") || !"AP_centered_mm" %in% rownames(ct)) {
    return(list(status = "model_failed", reason = "Satterthwaite_test_failed"))
  }
  r <- ct["AP_centered_mm", ]
  pick <- function(pattern) {
    n <- grep(pattern, names(r), value = TRUE, ignore.case = TRUE)
    if (length(n)) {
      as.numeric(r[[n[1]]])
    } else {
      NA_real_
    }
  }
  beta <- unname(coef(fit)["AP_centered_mm"])
  se <- pick("^SE$")
  stat <- pick("tstat|t_stat")
  df <- pick("df_Satt|df")
  p <- pick("p_Satt|p_val|p$")
  crit <- qt(0.975, df)
  getterm <- function(term) {
    if (!term %in% rownames(ct)) {
      return(c(beta = NA, se = NA))
    }
    c(beta = unname(coef(fit)[term]), se = as.numeric(ct[term, grep("^SE$", names(ct), value = TRUE)[1]]))
  }
  dw <- getterm("density_within")
  db <- getterm("density_between")
  list(
    status = "tested", reason = "", fit = fit, V = V, beta = beta, se = se, stat = stat, df = df,
    p = p, lo = beta - crit * se, hi = beta + crit * se, dw = dw, db = db, coef_cor = if (model_id ==
      "M1") {
      V["AP_centered_mm", "density_between"] / sqrt(V["AP_centered_mm", "AP_centered_mm"] *
        V["density_between", "density_between"])
    } else {
      NA_real_
    }
  )
}

axyd_design_condition <- function(d) {
  w <- fitted(glm(raw_count ~ absolute_tile + AP_centered_mm + density_within + density_between,
    offset = log(exposure),
    family = quasipoisson(link = "log"), data = d
  ))
  Z <- model.matrix(~absolute_tile, d)
  X <- as.matrix(d[c("AP_centered_mm", "density_within", "density_between")])
  sw <- sqrt(pmax(w, .Machine$double.eps))
  Zw <- Z * sw
  Xw <- X * sw
  R <- Xw - Zw %*% qr.coef(qr(Zw), Xw)
  R <- scale(R, center = FALSE, scale = sqrt(colSums(R^2)))
  if (any(!is.finite(R))) {
    Inf
  } else {
    kappa(R)
  }
}
