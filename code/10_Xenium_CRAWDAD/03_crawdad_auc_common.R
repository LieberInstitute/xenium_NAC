#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(crawdad)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(pracma)
  library(readr)
  library(tidyr)
})

BR6660_DEPTHS <- c(
  Br6660_NAc1_580 = 580,
  Br6660_NAc2_1090 = 1090,
  Br6660_NAc3_1580 = 1580,
  Br6660_NAc4_2080 = 2080,
  Br6660_NAc5_2580 = 2580,
  Br6660_NAc6_3080 = 3080,
  Br6660_NAc7_3580 = 3580,
  Br6660_Nac10_4080 = 4080,
  Br6660_NAc8_4580 = 4580,
  Br6660_NAc9_5080 = 5080,
  Br6660_Nac11_5580 = 5580
)

ROI_COMPARISON_SETS <- list(
  primary_anatomical = c("lateral", "dorsomedial", "ventromedial"),
  full_exclusive_partition = c(
    "lateral", "dorsomedial", "ventromedial", "outside_global_roi"
  ),
  coarse_exclusive_partition = c("lateral", "medial", "outside_global_roi")
)

ROI_COLORS <- c(
  lateral = "#D55E00",
  dorsomedial = "#0072B2",
  ventromedial = "#009E73",
  medial = "#6A3D9A",
  outside_global_roi = "#7F7F7F",
  whole_tissue = "#111111"
)

parse_cli_tokens <- function(args, scalar_flags, multi_flags,
                             boolean_flags = c("overwrite", "validate-only")) {
  out <- list()
  i <- 1L
  while (i <= length(args)) {
    token <- args[[i]]
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token)
    name <- substring(token, 3L)
    if (name %in% boolean_flags) {
      out[[name]] <- TRUE
      i <- i + 1L
      next
    }
    if (!name %in% c(scalar_flags, multi_flags)) stop("Unknown option: ", token)
    j <- i + 1L
    while (j <= length(args) && !startsWith(args[[j]], "--")) j <- j + 1L
    values <- if (j == i + 1L) character() else args[seq.int(i + 1L, j - 1L)]
    if (!length(values)) stop("Missing value for ", token)
    if (name %in% scalar_flags && length(values) != 1L) {
      stop(token, " expects exactly one value")
    }
    out[[name]] <- if (name %in% scalar_flags) values[[1L]] else values
    i <- j
  }
  out
}

resolve_requested_samples <- function(requested, available, donor = "Br6660",
                                      minimum_samples = 1L) {
  available <- unique(as.character(available))
  if (identical(requested, "all")) {
    selected <- available
  } else {
    selected <- as.character(requested)
    missing <- setdiff(selected, available)
    if (length(missing)) stop("Requested samples are unavailable: ", paste(missing, collapse = ", "))
  }
  if (donor == "Br6660") {
    selected <- names(BR6660_DEPTHS)[names(BR6660_DEPTHS) %in% selected]
  } else {
    selected <- sort(selected)
  }
  if (length(selected) < minimum_samples) {
    stop("At least ", minimum_samples, " selected sample(s) are required")
  }
  selected
}

pair_key <- function(reference, neighbor) {
  paste(reference, neighbor, sep = " -> ")
}

slugify <- function(x) {
  ans <- gsub("[^A-Za-z0-9._-]+", "_", x)
  ans <- gsub("^_+|_+$", "", ans)
  ifelse(nchar(ans) == 0L, "unnamed", ans)
}

atomic_finalize <- function(tmp, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop("Could not finalize output: ", path)
    }
    unlink(tmp)
  }
  invisible(path)
}

safe_write_table <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ext <- if (grepl("\\.csv\\.gz$", path, ignore.case = TRUE)) {
    ".csv.gz"
  } else if (grepl("\\.parquet$", path, ignore.case = TRUE)) {
    ".parquet"
  } else if (grepl("\\.csv$", path, ignore.case = TRUE)) {
    ".csv"
  } else {
    stop("Unsupported table output extension: ", path)
  }
  tmp <- file.path(
    dirname(path),
    paste0(".", basename(path), ".tmp.", Sys.getpid(), ext)
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  if (ext == ".parquet") {
    arrow::write_parquet(as.data.frame(x), tmp)
  } else {
    readr::write_csv(as.data.frame(x), tmp, na = "")
  }
  if (!file.exists(tmp) || file.info(tmp)$size <= 0L) {
    stop("Table write produced an empty output: ", path)
  }
  atomic_finalize(tmp, path)
}

safe_write_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp.", Sys.getpid()))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  saveRDS(x, tmp, compress = "xz")
  atomic_finalize(tmp, path)
}

safe_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp.", Sys.getpid()))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  jsonlite::write_json(
    x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null",
    digits = NA
  )
  atomic_finalize(tmp, path)
}

safe_write_lines <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- file.path(dirname(path), paste0(".", basename(path), ".tmp.", Sys.getpid()))
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  readr::write_lines(x, tmp)
  atomic_finalize(tmp, path)
}

safe_save_plot <- function(plot, base_path, formats, width, height,
                           plot_name, required = TRUE,
                           status = "completed", message = NA_character_) {
  formats <- intersect(formats, c("png", "pdf"))
  if (!length(formats)) stop("No supported plot formats requested")
  rows <- list()
  for (format in formats) {
    path <- paste0(base_path, ".", format)
    dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
    tmp <- file.path(
      dirname(path),
      paste0(".", tools::file_path_sans_ext(basename(path)), ".tmp.",
             Sys.getpid(), ".", format)
    )
    plot_status <- status
    error_message <- message
    if (!identical(status, "excluded")) {
      tryCatch(
        {
          if (format == "png") {
            ggplot2::ggsave(
              tmp, plot, width = width, height = height, dpi = 300, bg = "white"
            )
          } else {
            ggplot2::ggsave(
              tmp, plot, width = width, height = height,
              device = grDevices::cairo_pdf, bg = "white"
            )
          }
          if (!file.exists(tmp) || file.info(tmp)$size <= 0L) {
            stop("Plot device produced an empty file")
          }
          atomic_finalize(tmp, path)
        },
        error = function(e) {
          plot_status <<- "failed"
          error_message <<- conditionMessage(e)
          if (file.exists(tmp)) unlink(tmp)
        }
      )
    } else {
      path <- NA_character_
    }
    rows[[length(rows) + 1L]] <- data.frame(
      plot_name = plot_name,
      status = plot_status,
      path = path,
      width = width,
      height = height,
      format = format,
      required = required,
      error_message = error_message,
      stringsAsFactors = FALSE
    )
  }
  manifest <- dplyr::bind_rows(rows)
  if (required && any(manifest$status == "failed")) {
    stop(
      "Required plot failed: ", plot_name, "; ",
      paste(na.omit(manifest$error_message), collapse = "; ")
    )
  }
  manifest
}

discover_module02_auc_inputs <- function(input_root, regions,
                                         neighborhood_distance_um) {
  input_root <- normalizePath(input_root, mustWork = TRUE)
  rows <- list()
  for (region in regions) {
    distance_root <- file.path(
      input_root, paste0("neighdist_", neighborhood_distance_um)
    )
    region_root <- file.path(distance_root, region)
    search_root <- file.path(region_root, "combined")
    canonical_input <- file.path(search_root, "calculate_auc_input.parquet")
    candidates <- canonical_input[file.exists(canonical_input)]
    if (length(candidates) != 1L) {
      all_candidates <- if (dir.exists(search_root)) {
        list.files(search_root, pattern = "\\.parquet$", recursive = TRUE,
                   full.names = TRUE)
      } else {
        character()
      }
      stop(
        "Expected exactly one Module 02 AUC input for region=", region,
        ", neighdist=", neighborhood_distance_um, "; found ",
        length(candidates), ". Candidates under the region root: ",
        if (length(all_candidates)) paste(all_candidates, collapse = "; ") else "<none>"
      )
    }
    eligibility_file <- file.path(region_root, "qc", "reference_eligibility.csv")
    if (!file.exists(eligibility_file)) {
      stop("Missing Module 02 eligibility table: ", eligibility_file)
    }
    rows[[length(rows) + 1L]] <- data.frame(
      region = region,
      neighborhood_distance_um = neighborhood_distance_um,
      auc_input_file = normalizePath(candidates[[1]], mustWork = TRUE),
      eligibility_file = normalizePath(eligibility_file, mustWork = TRUE),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

resolve_alias <- function(names_in, canonical, aliases, required = TRUE) {
  candidates <- intersect(c(canonical, aliases), names_in)
  if (canonical %in% candidates) return(canonical)
  if (length(candidates) > 1L) {
    stop("Ambiguous aliases for ", canonical, ": ", paste(candidates, collapse = ", "))
  }
  if (!length(candidates)) {
    if (required) stop("Could not map required column: ", canonical)
    return(NA_character_)
  }
  candidates[[1]]
}

read_and_normalize_auc_input <- function(path) {
  raw <- as.data.frame(arrow::read_parquet(path))
  aliases <- list(
    sample = c("id"),
    region = c("analysis_region"),
    reference = character(),
    neighbor = character(),
    scale = c("scale_um"),
    Z = c("z", "z_score"),
    reference_eligible = character(),
    donor = character(),
    z_height = c("depth"),
    neighborhood_distance_um = character(),
    permutation = c("perm"),
    threshold_used = character(),
    task_id = character()
  )
  required <- c(
    "sample", "region", "reference", "neighbor", "scale", "Z",
    "reference_eligible", "z_height", "permutation"
  )
  map <- vapply(
    names(aliases),
    function(x) resolve_alias(names(raw), x, aliases[[x]], x %in% required),
    character(1)
  )
  out <- data.frame(row_id___ = seq_len(nrow(raw)))
  for (canonical in names(map)) {
    source <- map[[canonical]]
    if (!is.na(source)) out[[canonical]] <- raw[[source]]
  }
  out$row_id___ <- NULL
  out$source_file <- normalizePath(path, mustWork = TRUE)
  for (column in c("sample", "region", "reference", "neighbor")) {
    out[[column]] <- as.character(out[[column]])
    if (anyNA(out[[column]]) || any(!nzchar(out[[column]]))) {
      stop("Missing/empty values in ", column, " for ", path)
    }
  }
  for (column in c("scale", "Z", "z_height", "permutation")) {
    out[[column]] <- suppressWarnings(as.numeric(out[[column]]))
    if (column != "Z" && any(!is.finite(out[[column]]))) {
      stop("Non-finite values in ", column, " for ", path)
    }
  }
  if (!is.logical(out$reference_eligible)) {
    value <- tolower(as.character(out$reference_eligible))
    if (any(!value %in% c("true", "false"))) {
      stop("reference_eligible cannot be converted to logical in ", path)
    }
    out$reference_eligible <- value == "true"
  }
  list(
    data = out,
    column_map = data.frame(
      canonical = names(map), original = unname(map), stringsAsFactors = FALSE
    ),
    original_columns = names(raw)
  )
}

select_auc_scales <- function(dat, expected_scales) {
  if (!"scale" %in% names(dat)) stop("AUC input is missing scale")
  expected_scales <- sort(unique(as.numeric(expected_scales)))
  if (!length(expected_scales) || any(!is.finite(expected_scales))) {
    stop("Expected scales must be finite numeric values")
  }
  available_scales <- sort(unique(dat$scale[is.finite(dat$scale)]))
  missing_scales <- setdiff(expected_scales, available_scales)
  if (length(missing_scales)) {
    stop(
      "Requested AUC scales are absent from the Module 02 input: ",
      paste(missing_scales, collapse = ", ")
    )
  }
  selected <- dat |>
    dplyr::filter(scale %in% expected_scales)
  if (!nrow(selected)) stop("Scale filtering removed every AUC input row")
  list(
    data = selected,
    qc = data.frame(
      available_scales_um = paste(available_scales, collapse = ";"),
      selected_scales_um = paste(expected_scales, collapse = ";"),
      excluded_scales_um = paste(
        setdiff(available_scales, expected_scales), collapse = ";"
      ),
      n_rows_before_scale_filter = nrow(dat),
      n_rows_after_scale_filter = nrow(selected),
      n_rows_excluded_by_scale_filter = nrow(dat) - nrow(selected),
      stringsAsFactors = FALSE
    )
  )
}

read_reference_eligibility <- function(path) {
  raw <- suppressMessages(readr::read_csv(path, show_col_types = FALSE))
  required <- c(
    "sample", "depth", "analysis_region", "cell_type", "n_cells",
    "present", "reference_eligible"
  )
  missing <- setdiff(required, names(raw))
  if (length(missing)) {
    stop("Eligibility table is missing columns: ", paste(missing, collapse = ", "))
  }
  out <- raw |>
    dplyr::transmute(
      sample = as.character(sample),
      z_height = as.numeric(depth),
      region = as.character(analysis_region),
      reference = as.character(cell_type),
      n_reference_cells = as.numeric(n_cells),
      present = as.logical(present),
      reference_eligible = as.logical(reference_eligible),
      exclusion_reason = dplyr::case_when(
        !present ~ "reference_absent",
        !reference_eligible ~ "reference_cell_count_below_minimum",
        TRUE ~ NA_character_
      ),
      eligibility_source_file = normalizePath(path, mustWork = TRUE)
    )
  key <- c("sample", "region", "reference")
  if (anyDuplicated(out[key])) stop("Duplicate eligibility keys in ", path)
  out
}

validate_depth_metadata <- function(dat, donor = "Br6660", require_all = FALSE) {
  mapping <- dat |>
    dplyr::distinct(sample, z_height) |>
    dplyr::arrange(z_height, sample)
  duplicate_samples <- mapping |>
    dplyr::count(sample) |>
    dplyr::filter(n != 1L)
  if (nrow(duplicate_samples)) stop("Sample-to-depth mapping is not one-to-one")
  if (donor == "Br6660") {
    unknown <- setdiff(mapping$sample, names(BR6660_DEPTHS))
    if (length(unknown)) stop("Unexpected Br6660 sample IDs: ", paste(unknown, collapse = ", "))
    expected <- unname(BR6660_DEPTHS[mapping$sample])
    if (!isTRUE(all.equal(mapping$z_height, expected))) {
      stop("Observed Br6660 sample-to-depth mapping conflicts with the fixed expectation")
    }
    if (require_all && !setequal(mapping$sample, names(BR6660_DEPTHS))) {
      stop("The full Br6660 11-slice set is required but is incomplete")
    }
  }
  if (is.unsorted(mapping$z_height, strictly = TRUE)) {
    stop("Depth metadata are duplicated or non-monotonic")
  }
  mapping
}

validate_comparison_set <- function(regions) {
  regions <- unique(regions)
  if ("medial" %in% regions && any(c("dorsomedial", "ventromedial") %in% regions)) {
    stop("Overlapping ROI set rejected: medial cannot be combined with dorsomedial or ventromedial")
  }
  if ("whole_tissue" %in% regions && length(regions) > 1L) {
    stop("Overlapping ROI set rejected: whole_tissue cannot be combined with constituent ROIs")
  }
  invisible(TRUE)
}

validate_auc_curves <- function(dat, expected_scales, expected_permutations) {
  required <- c(
    "id", "reference", "neighbor", "scale", "Z", "permutation",
    "reference_eligible"
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) stop("Curve input is missing columns: ", paste(missing, collapse = ", "))
  duplicate_key <- c("id", "reference", "neighbor", "scale", "permutation")
  duplicated_rows <- duplicated(dat[duplicate_key]) | duplicated(dat[duplicate_key], fromLast = TRUE)
  duplicate_counts <- dat[duplicated_rows, , drop = FALSE] |>
    dplyr::count(id, reference, neighbor, name = "n_duplicate_rows")

  per_scale <- dat |>
    dplyr::group_by(id, reference, neighbor, scale) |>
    dplyr::summarise(
      permutations_present = paste(sort(unique(permutation)), collapse = ","),
      n_permutations_present_scale = dplyr::n_distinct(permutation),
      missing_permutations = paste(
        setdiff(expected_permutations, sort(unique(permutation))), collapse = ","
      ),
      scale_permutations_complete =
        setequal(sort(unique(permutation)), expected_permutations),
      .groups = "drop"
    )

  qc <- dat |>
    dplyr::group_by(id, reference, neighbor) |>
    dplyr::summarise(
      n_scales_expected = length(expected_scales),
      n_scales_present = dplyr::n_distinct(scale),
      scales_present = paste(sort(unique(scale)), collapse = ","),
      missing_scales = paste(setdiff(expected_scales, sort(unique(scale))), collapse = ","),
      unexpected_scales = paste(setdiff(sort(unique(scale)), expected_scales), collapse = ","),
      n_permutations_expected = length(expected_permutations),
      n_permutations_present = dplyr::n_distinct(permutation),
      n_nonfinite_z = sum(!is.finite(Z)),
      reference_eligible = all(reference_eligible %in% TRUE),
      .groups = "drop"
    ) |>
    dplyr::left_join(
      per_scale |>
        dplyr::group_by(id, reference, neighbor) |>
        dplyr::summarise(
          missing_permutations_by_scale = paste(
            paste0(scale[missing_permutations != ""], ":", missing_permutations[missing_permutations != ""]),
            collapse = ";"
          ),
          all_scales_have_expected_permutations = all(scale_permutations_complete),
          .groups = "drop"
        ),
      by = c("id", "reference", "neighbor")
    ) |>
    dplyr::left_join(
      duplicate_counts, by = c("id", "reference", "neighbor")
    ) |>
    dplyr::mutate(
      n_duplicate_rows = dplyr::coalesce(n_duplicate_rows, 0L),
      curve_complete =
        n_scales_present == n_scales_expected &
        missing_scales == "" &
        unexpected_scales == "" &
        all_scales_have_expected_permutations &
        n_nonfinite_z == 0L &
        n_duplicate_rows == 0L &
        reference_eligible,
      curve_incomplete_reason = dplyr::case_when(
        n_duplicate_rows > 0L ~ "duplicate_rows",
        !reference_eligible ~ "reference_ineligible",
        unexpected_scales != "" ~ "unexpected_scales",
        missing_scales != "" ~ "missing_scales",
        !all_scales_have_expected_permutations ~ "missing_permutations",
        n_nonfinite_z > 0L ~ "nonfinite_z",
        TRUE ~ NA_character_
      ),
      pair_key = pair_key(reference, neighbor)
    ) |>
    dplyr::arrange(id, reference, neighbor)
  qc
}

filter_complete_shared_pairs <- function(dat, curve_qc, ids, eligibility) {
  pairs <- dat |>
    dplyr::distinct(reference, neighbor) |>
    dplyr::arrange(reference, neighbor)
  coverage <- tidyr::crossing(
    id = as.character(ids), pairs
  ) |>
    dplyr::left_join(
      curve_qc |>
        dplyr::select(
          id, reference, neighbor, curve_complete, curve_incomplete_reason
        ),
      by = c("id", "reference", "neighbor")
    ) |>
    dplyr::left_join(
      eligibility |>
        dplyr::select(id, reference, reference_eligible, exclusion_reason),
      by = c("id", "reference")
    ) |>
    dplyr::mutate(
      reference_eligible = dplyr::coalesce(reference_eligible, FALSE),
      id_status = dplyr::case_when(
        !reference_eligible ~ dplyr::coalesce(
          exclusion_reason, "reference_eligibility_missing"
        ),
        is.na(curve_complete) ~ "curve_absent",
        !curve_complete ~ dplyr::coalesce(curve_incomplete_reason, "curve_incomplete"),
        TRUE ~ "complete"
      )
    )
  pair_qc <- coverage |>
    dplyr::group_by(reference, neighbor) |>
    dplyr::summarise(
      pair_key = paste(dplyr::first(reference), dplyr::first(neighbor), sep = " -> "),
      n_ids_expected = length(ids),
      n_ids_complete = sum(id_status == "complete"),
      pair_included = all(id_status == "complete"),
      exclusion_reasons = paste(
        paste0(id[id_status != "complete"], ":", id_status[id_status != "complete"]),
        collapse = ";"
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(reference, neighbor)
  retained <- pair_qc |>
    dplyr::filter(pair_included) |>
    dplyr::select(reference, neighbor)
  if (!nrow(retained)) stop("No complete shared directional pairs remain")
  filtered <- dat |>
    dplyr::inner_join(retained, by = c("reference", "neighbor")) |>
    dplyr::arrange(id, reference, neighbor, scale, permutation)
  if (any(!is.finite(filtered$Z))) {
    stop("Non-finite Z remained after complete-pair filtering")
  }
  list(data = filtered, pair_qc = pair_qc, coverage = coverage,
       shared_pairs = pair_qc |> dplyr::filter(pair_included))
}

run_official_crawdad_auc <- function(dat, id_order) {
  dat$id <- factor(as.character(dat$id), levels = id_order, ordered = TRUE)
  dat <- dat |>
    dplyr::arrange(id, reference, neighbor, scale, permutation) |>
    dplyr::mutate(id = as.character(id))
  dat_list <- split(dat, factor(dat$id, levels = id_order), drop = TRUE)
  auc <- crawdad::calculateAUC(dat_list, sharedPairs = TRUE) |>
    as.data.frame() |>
    dplyr::mutate(
      id = as.character(id),
      reference = as.character(reference),
      neighbor = as.character(neighbor),
      pair_key = pair_key(reference, neighbor)
    ) |>
    dplyr::arrange(match(id, id_order), reference, neighbor)
  counts <- auc |>
    dplyr::count(id, name = "n_pairs")
  if (nrow(counts) != length(id_order) || length(unique(counts$n_pairs)) != 1L) {
    stop("Official AUC output is unbalanced across IDs")
  }
  identities <- split(auc$pair_key, auc$id)
  reference_pairs <- sort(unique(identities[[1]]))
  if (!all(vapply(identities, function(x) identical(sort(unique(x)), reference_pairs), logical(1)))) {
    stop("Official AUC output does not contain identical pairs across IDs")
  }
  if (anyDuplicated(auc[c("id", "reference", "neighbor")])) {
    stop("Official AUC output has duplicate id-pair rows")
  }
  if (any(!is.finite(auc$auc))) stop("Official AUC output contains non-finite values")

  spot <- auc[1, , drop = FALSE]
  curve <- dat |>
    dplyr::filter(
      id == spot$id[[1]], reference == spot$reference[[1]],
      neighbor == spot$neighbor[[1]]
    ) |>
    dplyr::group_by(scale) |>
    dplyr::summarise(mean_z = mean(Z, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(scale)
  manual <- pracma::trapz(curve$scale, curve$mean_z)
  if (!isTRUE(all.equal(manual, spot$auc[[1]], tolerance = 1e-10))) {
    stop("Official AUC/manual trapz spot-check disagreement")
  }
  list(auc = auc, spot_check = data.frame(
    id = spot$id[[1]], reference = spot$reference[[1]],
    neighbor = spot$neighbor[[1]], official_auc = spot$auc[[1]],
    manual_auc = manual, absolute_difference = abs(spot$auc[[1]] - manual)
  ))
}

build_auc_matrix <- function(auc, id_order) {
  wide <- auc |>
    dplyr::select(id, pair_key, auc) |>
    tidyr::pivot_wider(names_from = pair_key, values_from = auc)
  wide <- wide[match(id_order, wide$id), , drop = FALSE]
  if (anyNA(wide)) stop("AUC matrix contains missing values")
  mat <- as.matrix(wide[, setdiff(names(wide), "id"), drop = FALSE])
  rownames(mat) <- wide$id
  list(matrix = mat, table = data.frame(id = wide$id, wide[, -1, drop = FALSE],
                                        check.names = FALSE))
}

compute_auc_variance_table <- function(auc) {
  auc |>
    dplyr::group_by(reference, neighbor, pair_key) |>
    dplyr::summarise(
      variance = stats::var(auc),
      mean_auc = mean(auc),
      median_auc = stats::median(auc),
      min_auc = min(auc),
      max_auc = max(auc),
      n_ids = dplyr::n_distinct(id),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(variance), reference, neighbor)
}

compute_pca_components <- function(auc_matrix) {
  empty_scores <- data.frame(id = character(), PC1 = numeric(), PC2 = numeric())
  empty_loadings <- data.frame(pair_key = character(), PC1 = numeric(), PC2 = numeric())
  empty_variance <- data.frame(
    component = character(), variance_explained = numeric(),
    cumulative_variance_explained = numeric()
  )
  variances <- apply(auc_matrix, 2, stats::var)
  if (nrow(auc_matrix) < 3L || sum(is.finite(variances) & variances > 0) < 2L) {
    return(list(computable = FALSE, reason =
      "PCA requires at least three IDs and two nonzero-variance AUC features",
      scores = empty_scores, loadings = empty_loadings, variance = empty_variance))
  }
  scaled <- scale(auc_matrix)
  if (any(!is.finite(scaled))) {
    return(list(computable = FALSE, reason =
      "Official scale() standardization produced non-finite values",
      scores = empty_scores, loadings = empty_loadings, variance = empty_variance))
  }
  pca <- stats::prcomp(scaled)
  scores <- as.data.frame(pca$x) |>
    tibble::rownames_to_column("id")
  loadings <- as.data.frame(pca$rotation) |>
    tibble::rownames_to_column("pair_key")
  explained <- pca$sdev^2 / sum(pca$sdev^2)
  variance <- data.frame(
    component = paste0("PC", seq_along(explained)),
    variance_explained = explained,
    cumulative_variance_explained = cumsum(explained)
  )
  list(computable = TRUE, reason = NA_character_, scores = scores,
       loadings = loadings, variance = variance, pca = pca)
}

placeholder_plot <- function(title, message) {
  ggplot2::ggplot() +
    ggplot2::annotate("text", x = 0, y = 0, label = message, size = 5) +
    ggplot2::xlim(-1, 1) + ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) + ggplot2::theme_void()
}

make_auc_heatmap <- function(auc, id_order, variance_table, title) {
  pair_order <- variance_table$pair_key
  auc |>
    dplyr::mutate(
      id = factor(id, levels = rev(id_order)),
      pair_key = factor(pair_key, levels = pair_order)
    ) |>
    ggplot2::ggplot(ggplot2::aes(x = pair_key, y = id, fill = auc)) +
    ggplot2::geom_tile() +
    ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B") +
    ggplot2::labs(x = "Directional pair", y = NULL, fill = "Signed AUC", title = title) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, hjust = 1, vjust = 0.5, size = 5)
    )
}

make_top_auc_plot <- function(auc, top_pairs, analysis_axis, id_metadata, title) {
  d <- auc |>
    dplyr::filter(pair_key %in% top_pairs)
  metadata_to_add <- setdiff(names(id_metadata), names(d))
  if (length(metadata_to_add)) {
    d <- dplyr::left_join(
      d, id_metadata |> dplyr::select(id, dplyr::all_of(metadata_to_add)), by = "id"
    )
  }
  if (analysis_axis == "depth") {
    ggplot2::ggplot(
      d, ggplot2::aes(x = z_height, y = auc, group = pair_key, color = pair_key)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey55") +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::facet_wrap(~pair_key, scales = "free_y") +
      ggplot2::labs(x = "AP depth (µm)", y = "Signed AUC", color = NULL, title = title) +
      ggplot2::theme_bw() + ggplot2::theme(legend.position = "none")
  } else {
    ggplot2::ggplot(
      d, ggplot2::aes(x = id, y = auc, color = id, group = pair_key)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "grey55") +
      ggplot2::geom_point(size = 2.5) +
      ggplot2::facet_wrap(~pair_key, scales = "free_y") +
      ggplot2::scale_color_manual(values = ROI_COLORS) +
      ggplot2::labs(x = "Region", y = "Signed AUC", color = NULL, title = title) +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1))
  }
}

make_multiscale_plot <- function(dat, selected_pairs, analysis_axis,
                                 id_metadata, title) {
  d <- dat |>
    dplyr::mutate(pair_key = pair_key(reference, neighbor)) |>
    dplyr::filter(pair_key %in% selected_pairs) |>
    dplyr::group_by(id, reference, neighbor, pair_key, scale) |>
    dplyr::summarise(mean_z = mean(Z), .groups = "drop")
  metadata_to_add <- setdiff(names(id_metadata), names(d))
  if (length(metadata_to_add)) {
    d <- dplyr::left_join(
      d, id_metadata |> dplyr::select(id, dplyr::all_of(metadata_to_add)), by = "id"
    )
  }
  if (analysis_axis == "depth") {
    ggplot2::ggplot(
      d, ggplot2::aes(x = scale, y = mean_z, group = id, color = z_height)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "black") +
      ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
      ggplot2::facet_wrap(~pair_key, scales = "free_y") +
      ggplot2::scale_color_viridis_c() +
      ggplot2::labs(x = "Shuffle scale (µm)", y = "Mean Z", color = "AP depth (µm)", title = title) +
      ggplot2::theme_bw()
  } else {
    ggplot2::ggplot(
      d, ggplot2::aes(x = scale, y = mean_z, group = id, color = id)
    ) +
      ggplot2::geom_hline(yintercept = 0, color = "black") +
      ggplot2::geom_line() + ggplot2::geom_point(size = 0.8) +
      ggplot2::facet_wrap(~pair_key, scales = "free_y") +
      ggplot2::scale_color_manual(values = ROI_COLORS) +
      ggplot2::labs(x = "Shuffle scale (µm)", y = "Mean Z", color = "Region", title = title) +
      ggplot2::theme_bw()
  }
}

write_task_metadata <- function(metadata, path) {
  safe_write_json(metadata, path)
}

run_auc_comparison_task <- function(dat, eligibility, id_order, id_metadata,
                                    analysis_axis, task_labels,
                                    expected_scales, expected_permutations,
                                    top_n_pairs, highlight_pairs,
                                    output_root, plot_root, plot_formats,
                                    input_manifest, column_maps,
                                    overwrite = FALSE) {
  started <- Sys.time()
  required_dirs <- c("inputs", "qc", "results", "metadata")
  if (dir.exists(output_root) && !overwrite) {
    stop("Output task directory exists; use --overwrite: ", output_root)
  }
  for (d in required_dirs) dir.create(file.path(output_root, d), recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_root, recursive = TRUE, showWarnings = FALSE)
  safe_write_table(input_manifest, file.path(output_root, "inputs", "input_manifest.csv"))

  curve_qc <- validate_auc_curves(dat, expected_scales, expected_permutations)
  safe_write_table(curve_qc, file.path(output_root, "qc", "curve_completeness.parquet"))
  if (any(curve_qc$n_duplicate_rows > 0L)) {
    stop("Duplicate id × reference × neighbor × scale × permutation rows detected")
  }
  filtered <- filter_complete_shared_pairs(dat, curve_qc, id_order, eligibility)
  official <- run_official_crawdad_auc(filtered$data, id_order)
  auc <- official$auc |>
    dplyr::left_join(id_metadata, by = "id")
  for (name in names(task_labels)) auc[[name]] <- task_labels[[name]]
  if (!"donor" %in% names(auc)) auc$donor <- unique(dat$donor)[[1]]
  auc <- auc |>
    dplyr::relocate(
      dplyr::any_of(c("donor", names(task_labels), "id", "sample", "z_height", "region")),
      reference, neighbor, pair, pair_key, auc
    )

  matrix_result <- build_auc_matrix(auc, id_order)
  variance <- compute_auc_variance_table(auc)
  top <- head(variance, top_n_pairs)
  pca <- compute_pca_components(matrix_result$matrix)

  safe_write_table(filtered$pair_qc, file.path(output_root, "qc", "pair_inclusion_exclusion.csv"))
  safe_write_table(eligibility, file.path(output_root, "qc", "reference_eligibility_summary.csv"))
  safe_write_table(filtered$shared_pairs, file.path(output_root, "qc", "shared_pairs.csv"))
  safe_write_table(official$spot_check, file.path(output_root, "qc", "official_auc_spot_check.csv"))
  safe_write_table(auc, file.path(output_root, "results", "auc_samples.parquet"))
  safe_write_table(auc, file.path(output_root, "results", "auc_samples.csv.gz"))
  safe_write_rds(auc, file.path(output_root, "results", "auc_samples.rds"))
  safe_write_table(matrix_result$table, file.path(output_root, "results", "auc_matrix_id_by_pair.csv.gz"))
  safe_write_table(variance, file.path(output_root, "results", "auc_variance_by_pair.csv"))
  safe_write_table(pca$scores, file.path(output_root, "results", "pca_scores.csv"))
  safe_write_table(pca$loadings, file.path(output_root, "results", "pca_loadings.csv"))
  safe_write_table(pca$variance, file.path(output_root, "results", "pca_variance_explained.csv"))
  safe_write_table(top, file.path(output_root, "results", "top_variable_pairs.csv"))
  capture.output(sessionInfo(), file = file.path(output_root, "metadata", "session_info.txt"))

  plot_manifest <- list()
  add_plot <- function(...) {
    plot_manifest[[length(plot_manifest) + 1L]] <<- safe_save_plot(...)
  }
  add_plot(
    crawdad::vizVarianceSamples(official$auc),
    file.path(plot_root, "official_variance_samples"), plot_formats, 8, 7,
    "official_variance_samples", required = TRUE
  )
  if (pca$computable) {
    official_pca <- tryCatch(crawdad::vizPCASamples(official$auc), error = identity)
    if (inherits(official_pca, "error")) {
      add_plot(
        placeholder_plot("Official CRAWDAD PCA", conditionMessage(official_pca)),
        file.path(plot_root, "official_pca_samples"), plot_formats, 7, 6,
        "official_pca_samples", required = FALSE,
        status = "not_computable", message = conditionMessage(official_pca)
      )
    } else {
      add_plot(
        official_pca, file.path(plot_root, "official_pca_samples"),
        plot_formats, 7, 6, "official_pca_samples", required = FALSE
      )
    }
    score_plot_data <- pca$scores |>
      dplyr::left_join(id_metadata, by = "id")
    if (analysis_axis == "depth") {
      project_pca <- ggplot2::ggplot(
        score_plot_data, ggplot2::aes(PC1, PC2, color = z_height, label = id)
      ) + ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.7, size = 3) +
        ggplot2::scale_color_viridis_c() + ggplot2::coord_equal() +
        ggplot2::labs(color = "AP depth (µm)", title = "PCA of signed CRAWDAD AUC") +
        ggplot2::theme_bw()
    } else {
      project_pca <- ggplot2::ggplot(
        score_plot_data, ggplot2::aes(PC1, PC2, color = id, label = id)
      ) + ggplot2::geom_point(size = 3) + ggplot2::geom_text(vjust = -0.7, size = 3) +
        ggplot2::scale_color_manual(values = ROI_COLORS) + ggplot2::coord_equal() +
        ggplot2::labs(color = "Region", title = "PCA of signed CRAWDAD AUC") +
        ggplot2::theme_bw()
    }
    add_plot(
      project_pca, file.path(plot_root, "project_pca_samples"),
      plot_formats, 8, 6, "project_pca_samples", required = TRUE
    )
  } else {
    for (name in c("official_pca_samples", "project_pca_samples")) {
      add_plot(
        placeholder_plot("PCA not computable", pca$reason),
        file.path(plot_root, name), plot_formats, 7, 6, name,
        required = FALSE, status = "not_computable", message = pca$reason
      )
    }
  }
  add_plot(
    make_auc_heatmap(auc, id_order, variance, "Signed AUC by directional pair"),
    file.path(plot_root, "auc_heatmap"), plot_formats, 18, 6,
    "auc_heatmap", required = TRUE
  )
  add_plot(
    make_top_auc_plot(auc, top$pair_key, analysis_axis, id_metadata,
                      "Top variable directional-pair AUC"),
    file.path(plot_root, "top_variable_pair_auc"), plot_formats, 12, 8,
    "top_variable_pair_auc", required = TRUE
  )
  add_plot(
    make_multiscale_plot(filtered$data, top$pair_key, analysis_axis,
                         id_metadata, "Multi-scale Z trends for top AUC-variable pairs"),
    file.path(plot_root, "top_variable_pair_multiscale_z"), plot_formats, 12, 8,
    "top_variable_pair_multiscale_z", required = TRUE
  )

  for (highlight in highlight_pairs) {
    parts <- strsplit(highlight, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L) stop("Highlight pair must be reference:neighbor: ", highlight)
    key <- pair_key(parts[[1]], parts[[2]])
    plot_name <- paste0("highlight_", slugify(parts[[1]]), "__to__", slugify(parts[[2]]))
    if (key %in% filtered$shared_pairs$pair_key) {
      add_plot(
        make_multiscale_plot(filtered$data, key, analysis_axis, id_metadata,
                             paste0("Directional relationship: ", key)),
        file.path(plot_root, plot_name), plot_formats, 9, 6,
        plot_name, required = TRUE
      )
    } else {
      reason <- filtered$pair_qc |>
        dplyr::filter(reference == parts[[1]], neighbor == parts[[2]]) |>
        dplyr::pull(exclusion_reasons)
      add_plot(
        placeholder_plot(paste0("Excluded: ", key), paste(reason, collapse = "; ")),
        file.path(plot_root, plot_name), plot_formats, 9, 6,
        plot_name, required = FALSE, status = "excluded",
        message = if (length(reason)) paste(reason, collapse = "; ") else "pair_absent"
      )
    }
  }
  plot_manifest_df <- dplyr::bind_rows(plot_manifest)
  safe_write_table(plot_manifest_df, file.path(output_root, "qc", "plot_manifest.csv"))
  errors <- plot_manifest_df |>
    dplyr::filter(status == "failed")
  safe_write_lines(
    if (nrow(errors)) paste(errors$plot_name, errors$error_message, sep = "\t") else character(),
    file.path(output_root, "qc", "plot_errors.log")
  )

  ended <- Sys.time()
  reason_counts <- filtered$pair_qc |>
    dplyr::filter(!pair_included) |>
    dplyr::count(exclusion_reasons, name = "n_pairs")
  metadata <- c(
    list(
      module_name = "module03_crawdad_auc_across_samples",
      analysis_axis = analysis_axis,
      donor = unique(dat$donor),
      ids = id_order,
      id_metadata = id_metadata,
      neighborhood_distance_um = unique(dat$neighborhood_distance_um),
      expected_scales = expected_scales,
      observed_scales = sort(unique(dat$scale)),
      expected_permutations = expected_permutations,
      observed_permutations = sort(unique(dat$permutation)),
      n_input_rows = nrow(dat),
      n_eligible_rows = sum(dat$reference_eligible),
      n_candidate_directional_pairs = nrow(filtered$pair_qc),
      n_complete_shared_directional_pairs = nrow(filtered$shared_pairs),
      excluded_pair_reasons = reason_counts,
      crawdad_calls = c(
        "crawdad::calculateAUC(dat_list, sharedPairs = TRUE)",
        "crawdad::vizVarianceSamples(auc_samples)",
        "crawdad::vizPCASamples(auc_samples)"
      ),
      sharedPairs = TRUE,
      crawdad_version = as.character(utils::packageVersion("crawdad")),
      package_versions = list(
        arrow = as.character(utils::packageVersion("arrow")),
        crawdad = as.character(utils::packageVersion("crawdad")),
        dplyr = as.character(utils::packageVersion("dplyr")),
        ggplot2 = as.character(utils::packageVersion("ggplot2")),
        pracma = as.character(utils::packageVersion("pracma")),
        tidyr = as.character(utils::packageVersion("tidyr"))
      ),
      R_version = R.version.string,
      source_files = unique(dat$source_file),
      column_maps = column_maps,
      output_root = normalizePath(output_root, mustWork = TRUE),
      plot_root = normalizePath(plot_root, mustWork = TRUE),
      output_files = sort(list.files(output_root, recursive = TRUE, full.names = TRUE)),
      plot_files = sort(list.files(plot_root, recursive = TRUE, full.names = TRUE)),
      start_time = format(started, tz = "America/New_York"),
      end_time = format(ended, tz = "America/New_York"),
      elapsed_seconds = as.numeric(difftime(ended, started, units = "secs")),
      final_status = "completed"
    ),
    task_labels
  )
  write_task_metadata(metadata, file.path(output_root, "metadata", "task_metadata.json"))
  list(
    status = "completed", auc = auc, variance = variance, top = top,
    pair_qc = filtered$pair_qc, plot_manifest = plot_manifest_df,
    metadata = metadata
  )
}

run_internal_validation_tests <- function() {
  expect_error <- function(expr) {
    inherits(try(force(expr), silent = TRUE), "try-error")
  }
  stopifnot(pair_key("A", "B") != pair_key("B", "A"))
  stopifnot(expect_error(validate_comparison_set(c("medial", "dorsomedial"))))
  stopifnot(expect_error(validate_comparison_set(c("whole_tissue", "lateral"))))
  validate_comparison_set(ROI_COMPARISON_SETS$primary_anatomical)
  depth <- data.frame(sample = names(BR6660_DEPTHS), z_height = unname(BR6660_DEPTHS))
  stopifnot(nrow(validate_depth_metadata(depth, require_all = TRUE)) == 11L)

  synth <- tidyr::crossing(
    id = c("id1", "id2"), reference = c("A", "B"), neighbor = c("A", "B"),
    scale = c(100, 200), permutation = c(1, 2)
  ) |>
    dplyr::mutate(Z = seq_len(dplyr::n()) / 10, reference_eligible = TRUE)
  qc <- validate_auc_curves(synth, c(100, 200), c(1, 2))
  stopifnot(all(qc$curve_complete))
  synth_with_extra_scale <- dplyr::bind_rows(
    synth,
    synth |>
      dplyr::filter(scale == 100) |>
      dplyr::mutate(scale = 50)
  )
  scale_selected <- select_auc_scales(synth_with_extra_scale, c(100, 200))
  stopifnot(
    setequal(unique(scale_selected$data$scale), c(100, 200)),
    identical(scale_selected$qc$excluded_scales_um, "50"),
    all(validate_auc_curves(
      scale_selected$data, c(100, 200), c(1, 2)
    )$curve_complete)
  )
  missing_scale <- synth |>
    dplyr::filter(!(id == "id2" & reference == "A" & neighbor == "B" & scale == 200))
  qc_scale <- validate_auc_curves(missing_scale, c(100, 200), c(1, 2))
  stopifnot(any(qc_scale$curve_incomplete_reason == "missing_scales"))
  missing_perm <- synth |>
    dplyr::filter(!(id == "id2" & reference == "A" & neighbor == "B" &
                      scale == 200 & permutation == 2))
  qc_perm <- validate_auc_curves(missing_perm, c(100, 200), c(1, 2))
  stopifnot(any(qc_perm$curve_incomplete_reason == "missing_permutations"))

  eligibility <- tidyr::crossing(id = c("id1", "id2"), reference = c("A", "B")) |>
    dplyr::mutate(
      reference_eligible = !(id == "id2" & reference == "A"),
      exclusion_reason = dplyr::if_else(
        reference_eligible, NA_character_, "reference_cell_count_below_minimum"
      )
    )
  filtered <- filter_complete_shared_pairs(
    synth, qc, c("id1", "id2"), eligibility
  )
  stopifnot(!any(filtered$shared_pairs$reference == "A"))
  stopifnot(any(filtered$shared_pairs$neighbor == "A"))
  filtered_missing <- filter_complete_shared_pairs(
    missing_scale, qc_scale, c("id1", "id2"),
    eligibility |> dplyr::mutate(reference_eligible = TRUE, exclusion_reason = NA_character_)
  )
  stopifnot(!any(filtered_missing$shared_pairs$reference == "A" &
                   filtered_missing$shared_pairs$neighbor == "B"))

  official <- run_official_crawdad_auc(synth, c("id1", "id2"))
  stopifnot(official$spot_check$absolute_difference < 1e-10)
  depth_ids <- transform(synth, id = as.character(id))
  region_ids <- transform(synth, id = ifelse(id == "id1", "lateral", "dorsomedial"))
  stopifnot(all(unique(depth_ids$id) == c("id1", "id2")))
  stopifnot(setequal(unique(region_ids$id), c("lateral", "dorsomedial")))
  data.frame(
    test = c(
      "directionality", "reference_only_eligibility", "missing_scale",
      "missing_permutation", "shared_pair_exclusion", "official_auc_spot_check",
      "axis_specific_ids", "reject_medial_overlap", "reject_whole_tissue_overlap",
      "Br6660_depth_order", "explicit_scale_filtering"
    ),
    status = "passed",
    stringsAsFactors = FALSE
  )
}
