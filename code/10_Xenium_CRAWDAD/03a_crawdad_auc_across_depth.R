#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "10_Xenium_CRAWDAD/03a_crawdad_auc_across_depth.R"
source(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), "03_crawdad_auc_common.R"))

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
defaults <- list(
  `input-root` = file.path(project_root, "..", "processed-data", "10_Xenium_CRAWDAD", "module02_crawdad_within_depth"),
  donor = "Br6660",
  regions = c("lateral", "dorsomedial", "ventromedial", "outside_global_roi"),
  samples = "all",
  `neighborhood-distance-um` = "50",
  `expected-scales-um` = as.character(seq(200, 1000, 100)),
  `expected-permutations` = "3",
  `top-n-pairs` = "20",
  `highlight-pairs` = c("D1_Island_A:D1_Island_B", "D1_Island_B:D1_Island_A"),
  `output-dir` = file.path(project_root, "..", "processed-data", "10_Xenium_CRAWDAD", "module03_crawdad_auc_across_depth"),
  `plot-dir` = file.path(project_root, "..", "plots", "10_Xenium_CRAWDAD", "module03_crawdad_auc_across_depth"),
  `plot-formats` = "png", overwrite = FALSE, `validate-only` = FALSE
)
parsed <- parse_cli_tokens(
  commandArgs(trailingOnly = TRUE),
  scalar_flags = c("input-root", "donor", "neighborhood-distance-um", "top-n-pairs", "output-dir", "plot-dir"),
  multi_flags = c("regions", "samples", "expected-scales-um", "expected-permutations", "highlight-pairs", "plot-formats")
)
opt <- utils::modifyList(defaults, parsed)
opt$`neighborhood-distance-um` <- as.numeric(opt$`neighborhood-distance-um`)
opt$`expected-scales-um` <- as.numeric(opt$`expected-scales-um`)
perm_arg <- as.integer(opt$`expected-permutations`)
opt$`expected-permutations` <- if (length(perm_arg) == 1L) seq_len(perm_arg) else perm_arg
opt$`top-n-pairs` <- as.integer(opt$`top-n-pairs`)
opt$`plot-formats` <- unique(tolower(opt$`plot-formats`))

allowed_regions <- c("lateral", "dorsomedial", "ventromedial", "medial", "outside_global_roi", "whole_tissue")
if (!length(opt$regions) || any(!opt$regions %in% allowed_regions)) stop("Invalid --regions value")
if (!is.finite(opt$`neighborhood-distance-um`) || opt$`neighborhood-distance-um` <= 0) stop("Invalid neighborhood distance")
if (!length(opt$`expected-scales-um`) || any(!is.finite(opt$`expected-scales-um`)) || anyDuplicated(opt$`expected-scales-um`)) stop("Invalid expected scales")
if (!length(opt$`expected-permutations`) || any(opt$`expected-permutations` <= 0) || anyDuplicated(opt$`expected-permutations`)) stop("Invalid expected permutations")
if (!is.finite(opt$`top-n-pairs`) || opt$`top-n-pairs` < 1L) stop("Invalid --top-n-pairs")
if (!identical(opt$`plot-formats`, "png")) stop("Module 03 plots must use --plot-formats png")
invisible(lapply(opt$`highlight-pairs`, function(x) if (length(strsplit(x, ":", fixed = TRUE)[[1]]) != 2L) stop("Invalid highlight pair: ", x)))
distance_token <- paste0("neighdist_", opt$`neighborhood-distance-um`)
opt$`output-dir` <- file.path(opt$`output-dir`, distance_token)
opt$`plot-dir` <- file.path(opt$`plot-dir`, distance_token)

if (isTRUE(opt$`validate-only`)) {
  print(run_internal_validation_tests())
  cat("Module 03 shared validation passed.\n")
  quit(status = 0L)
}

loaded <- list()
for (region in opt$regions) {
  manifest <- discover_module02_auc_inputs(
    opt$`input-root`, region, opt$`neighborhood-distance-um`
  )
  normalized <- read_and_normalize_auc_input(manifest$auc_input_file[[1]])
  scale_selection <- select_auc_scales(
    normalized$data, opt$`expected-scales-um`
  )
  dat <- scale_selection$data
  if ("donor" %in% names(dat) && any(dat$donor != opt$donor)) stop("Donor mismatch in ", region)
  if (any(dat$region != region)) stop("Region mismatch in ", region)
  if ("neighborhood_distance_um" %in% names(dat) &&
      any(dat$neighborhood_distance_um != opt$`neighborhood-distance-um`)) stop("Neighborhood-distance mismatch in ", region)
  eligibility <- read_reference_eligibility(manifest$eligibility_file[[1]])
  if (any(eligibility$region != region)) stop("Eligibility-region mismatch in ", region)
  manifest <- dplyr::bind_cols(manifest, scale_selection$qc)
  loaded[[region]] <- list(
    data = dat, eligibility = eligibility, manifest = manifest,
    column_map = normalized$column_map
  )
}

all_depths <- dplyr::bind_rows(lapply(loaded, function(x) x$data[c("sample", "z_height")])) |>
  dplyr::distinct()
available_samples <- Reduce(intersect, lapply(loaded, function(x) unique(x$data$sample)))
selected_samples <- resolve_requested_samples(opt$samples, available_samples, opt$donor, 2L)
validate_depth_metadata(
  all_depths |> dplyr::filter(sample %in% selected_samples), opt$donor,
  require_all = identical(opt$samples, "all")
)
for (region in names(loaded)) {
  observed <- loaded[[region]]$data |>
    dplyr::filter(sample %in% selected_samples) |>
    dplyr::distinct(sample, z_height)
  if (nrow(observed) != length(selected_samples)) stop("Incomplete sample coverage for ", region)
}

cat("Resolved Module 03 across-depth plan:\n")
cat("  tasks:", length(opt$regions), "\n")
cat("  regions:", paste(opt$regions, collapse = ", "), "\n")
cat("  IDs per task:", length(selected_samples), paste(selected_samples, collapse = ", "), "\n")
cat("  expected scales:", paste(opt$`expected-scales-um`, collapse = ", "), "\n")
cat("  output root:", opt$`output-dir`, "\n")
cat("  plot root:", opt$`plot-dir`, "\n")

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)
results <- list()
task_rows <- list()
for (region in opt$regions) {
  x <- loaded[[region]]
  dat <- x$data |>
    dplyr::filter(sample %in% selected_samples) |>
    dplyr::mutate(id = sample, donor = opt$donor) |>
    dplyr::arrange(match(sample, selected_samples), reference, neighbor, scale, permutation)
  eligibility <- x$eligibility |>
    dplyr::filter(sample %in% selected_samples) |>
    dplyr::mutate(id = sample)
  id_metadata <- dat |>
    dplyr::distinct(id, sample, z_height) |>
    dplyr::arrange(match(id, selected_samples))
  input_manifest <- x$manifest |>
    dplyr::mutate(
      donor = opt$donor,
      selected_samples = paste(.env$selected_samples, collapse = ";"),
      n_selected_samples = length(.env$selected_samples),
      n_selected_rows = nrow(dat)
    )
  started <- Sys.time()
  ans <- tryCatch(
    run_auc_comparison_task(
      dat, eligibility, selected_samples, id_metadata, "depth",
      task_labels = list(region = region),
      expected_scales = opt$`expected-scales-um`,
      expected_permutations = opt$`expected-permutations`,
      top_n_pairs = opt$`top-n-pairs`, highlight_pairs = opt$`highlight-pairs`,
      output_root = file.path(opt$`output-dir`, region),
      plot_root = file.path(opt$`plot-dir`, region),
      plot_formats = opt$`plot-formats`, input_manifest = input_manifest,
      column_maps = list(x$column_map), overwrite = isTRUE(opt$overwrite)
    ),
    error = identity
  )
  if (inherits(ans, "error")) {
    status <- "failed"
    message <- conditionMessage(ans)
    warning("Across-depth task failed for ", region, ": ", message)
  } else {
    status <- "completed"
    message <- NA_character_
    results[[region]] <- ans
  }
  task_rows[[length(task_rows) + 1L]] <- data.frame(
    region = region, status = status, n_ids = length(selected_samples),
    started = format(started, tz = "America/New_York"),
    ended = format(Sys.time(), tz = "America/New_York"), error_message = message
  )
}

task_manifest <- dplyr::bind_rows(task_rows)
combined_dir <- file.path(opt$`output-dir`, "combined")
dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
safe_write_table(task_manifest, file.path(combined_dir, "task_manifest.csv"))
if (length(results)) {
  safe_write_table(
    dplyr::bind_rows(lapply(results, `[[`, "auc")) |>
      dplyr::arrange(region, z_height, reference, neighbor),
    file.path(combined_dir, "auc_samples_all_regions.parquet")
  )
  safe_write_table(
    dplyr::bind_rows(lapply(names(results), function(region) {
      results[[region]]$pair_qc |> dplyr::mutate(region = region, .before = 1)
    })),
    file.path(combined_dir, "pair_inclusion_exclusion_all_regions.csv")
  )
}
if (any(task_manifest$status == "failed")) quit(status = 1L)
cat("All requested across-depth tasks completed.\n")
