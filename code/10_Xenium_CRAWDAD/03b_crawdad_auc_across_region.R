#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg)) sub("^--file=", "", script_arg[[1]]) else "10_Xenium_CRAWDAD/03b_crawdad_auc_across_region.R"
source(file.path(dirname(normalizePath(script_path, mustWork = TRUE)), "03_crawdad_auc_common.R"))

project_root <- normalizePath(file.path(dirname(script_path), ".."), mustWork = TRUE)
defaults <- list(
  `input-root` = file.path(project_root, "..", "processed-data", "10_Xenium_CRAWDAD", "module02_crawdad_within_depth"),
  donor = "Br6660", `comparison-sets` = "full_exclusive_partition", samples = "all",
  `neighborhood-distance-um` = "50",
  `expected-scales-um` = as.character(seq(200, 1000, 100)),
  `expected-permutations` = "3", `top-n-pairs` = "20",
  `highlight-pairs` = c("D1_Island_A:D1_Island_B", "D1_Island_B:D1_Island_A"),
  `output-dir` = file.path(project_root, "..", "processed-data", "10_Xenium_CRAWDAD", "module03_crawdad_auc_across_region"),
  `plot-dir` = file.path(project_root, "..", "plots", "10_Xenium_CRAWDAD", "module03_crawdad_auc_across_region"),
  `plot-formats` = "png", overwrite = FALSE, `validate-only` = FALSE
)
parsed <- parse_cli_tokens(
  commandArgs(trailingOnly = TRUE),
  scalar_flags = c("input-root", "donor", "neighborhood-distance-um", "top-n-pairs", "output-dir", "plot-dir"),
  multi_flags = c("comparison-sets", "samples", "expected-scales-um", "expected-permutations", "highlight-pairs", "plot-formats")
)
opt <- utils::modifyList(defaults, parsed)
opt$`neighborhood-distance-um` <- as.numeric(opt$`neighborhood-distance-um`)
opt$`expected-scales-um` <- as.numeric(opt$`expected-scales-um`)
perm_arg <- as.integer(opt$`expected-permutations`)
opt$`expected-permutations` <- if (length(perm_arg) == 1L) seq_len(perm_arg) else perm_arg
opt$`top-n-pairs` <- as.integer(opt$`top-n-pairs`)
opt$`plot-formats` <- unique(tolower(opt$`plot-formats`))

if (!length(opt$`comparison-sets`) || any(!opt$`comparison-sets` %in% names(ROI_COMPARISON_SETS))) stop("Invalid --comparison-sets value")
invisible(lapply(opt$`comparison-sets`, function(x) validate_comparison_set(ROI_COMPARISON_SETS[[x]])))
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

required_regions <- unique(unlist(ROI_COMPARISON_SETS[opt$`comparison-sets`], use.names = FALSE))
loaded <- list()
for (region in required_regions) {
  manifest <- discover_module02_auc_inputs(opt$`input-root`, region, opt$`neighborhood-distance-um`)
  normalized <- read_and_normalize_auc_input(manifest$auc_input_file[[1]])
  scale_selection <- select_auc_scales(
    normalized$data, opt$`expected-scales-um`
  )
  dat <- scale_selection$data
  if ("donor" %in% names(dat) && any(dat$donor != opt$donor)) stop("Donor mismatch in ", region)
  if (any(dat$region != region)) stop("Region mismatch in ", region)
  if ("neighborhood_distance_um" %in% names(dat) && any(dat$neighborhood_distance_um != opt$`neighborhood-distance-um`)) stop("Neighborhood-distance mismatch in ", region)
  eligibility <- read_reference_eligibility(manifest$eligibility_file[[1]])
  if (any(eligibility$region != region)) stop("Eligibility-region mismatch in ", region)
  manifest <- dplyr::bind_cols(manifest, scale_selection$qc)
  loaded[[region]] <- list(data = dat, eligibility = eligibility, manifest = manifest,
                           column_map = normalized$column_map)
}

available_samples <- Reduce(intersect, lapply(loaded, function(x) unique(x$data$sample)))
selected_samples <- resolve_requested_samples(opt$samples, available_samples, opt$donor, 1L)
depth_maps <- lapply(loaded, function(x) x$data |>
  dplyr::filter(sample %in% selected_samples) |>
  dplyr::distinct(sample, z_height))
canonical_depths <- validate_depth_metadata(dplyr::bind_rows(depth_maps) |> dplyr::distinct(), opt$donor,
                                            require_all = identical(opt$samples, "all"))
for (region in names(depth_maps)) {
  if (!setequal(paste(depth_maps[[region]]$sample, depth_maps[[region]]$z_height),
                paste(canonical_depths$sample, canonical_depths$z_height))) {
    stop("Sample/depth metadata conflict or missing coverage for ", region)
  }
}

n_tasks <- length(selected_samples) * length(opt$`comparison-sets`)
cat("Resolved Module 03 across-region plan:\n")
cat("  tasks:", n_tasks, "\n")
cat("  comparison sets:", paste(opt$`comparison-sets`, collapse = ", "), "\n")
for (set_name in opt$`comparison-sets`) cat("   ", set_name, ":", paste(ROI_COMPARISON_SETS[[set_name]], collapse = ", "), "\n")
cat("  slice IDs:", length(selected_samples), paste(selected_samples, collapse = ", "), "\n")
cat("  expected scales:", paste(opt$`expected-scales-um`, collapse = ", "), "\n")
cat("  output root:", opt$`output-dir`, "\n")
cat("  plot root:", opt$`plot-dir`, "\n")

dir.create(opt$`output-dir`, recursive = TRUE, showWarnings = FALSE)
dir.create(opt$`plot-dir`, recursive = TRUE, showWarnings = FALSE)
results <- list()
task_rows <- list()
for (set_name in opt$`comparison-sets`) {
  regions <- ROI_COMPARISON_SETS[[set_name]]
  for (sample_id in selected_samples) {
    depth <- unname(BR6660_DEPTHS[[sample_id]])
    dat <- dplyr::bind_rows(lapply(regions, function(region) loaded[[region]]$data |>
      dplyr::filter(sample == sample_id))) |>
      dplyr::mutate(id = region, donor = opt$donor) |>
      dplyr::arrange(match(region, regions), reference, neighbor, scale, permutation)
    if (!setequal(unique(dat$region), regions)) stop("Missing region data for ", set_name, "/", sample_id)
    eligibility <- dplyr::bind_rows(lapply(regions, function(region) loaded[[region]]$eligibility |>
      dplyr::filter(sample == sample_id))) |>
      dplyr::mutate(id = region)
    id_metadata <- data.frame(id = regions, region = regions, stringsAsFactors = FALSE)
    input_manifest <- dplyr::bind_rows(lapply(regions, function(region) loaded[[region]]$manifest)) |>
      dplyr::mutate(donor = opt$donor, comparison_set = set_name, sample = sample_id,
                    z_height = depth, n_selected_rows = vapply(regions, function(r) sum(dat$region == r), integer(1)))
    key <- paste(set_name, sample_id, sep = "::")
    started <- Sys.time()
    ans <- tryCatch(
      run_auc_comparison_task(
        dat, eligibility, regions, id_metadata, "region",
        task_labels = list(comparison_set = set_name, sample = sample_id, z_height = depth),
        expected_scales = opt$`expected-scales-um`, expected_permutations = opt$`expected-permutations`,
        top_n_pairs = opt$`top-n-pairs`, highlight_pairs = opt$`highlight-pairs`,
        output_root = file.path(opt$`output-dir`, set_name, sample_id),
        plot_root = file.path(opt$`plot-dir`, set_name, sample_id),
        plot_formats = opt$`plot-formats`, input_manifest = input_manifest,
        column_maps = lapply(regions, function(r) loaded[[r]]$column_map),
        overwrite = isTRUE(opt$overwrite)
      ), error = identity
    )
    if (inherits(ans, "error")) {
      status <- "failed"; message <- conditionMessage(ans)
      warning("Across-region task failed for ", set_name, "/", sample_id, ": ", message)
    } else {
      status <- "completed"; message <- NA_character_; results[[key]] <- ans
    }
    task_rows[[length(task_rows) + 1L]] <- data.frame(
      comparison_set = set_name, sample = sample_id, z_height = depth,
      status = status, n_ids = length(regions),
      started = format(started, tz = "America/New_York"), ended = format(Sys.time(), tz = "America/New_York"),
      error_message = message
    )
  }
}

task_manifest <- dplyr::bind_rows(task_rows)
safe_write_table(task_manifest, file.path(opt$`output-dir`, "task_manifest.csv"))
for (set_name in opt$`comparison-sets`) {
  set_keys <- grep(paste0("^", set_name, "::"), names(results), value = TRUE)
  if (!length(set_keys)) next
  combined_dir <- file.path(opt$`output-dir`, set_name, "combined")
  dir.create(combined_dir, recursive = TRUE, showWarnings = FALSE)
  auc_all <- dplyr::bind_rows(lapply(results[set_keys], `[[`, "auc")) |>
    dplyr::arrange(z_height, region, reference, neighbor)
  pair_rows <- dplyr::bind_rows(lapply(set_keys, function(key) {
    sample_id <- sub("^[^:]+::", "", key)
    results[[key]]$pair_qc |>
      dplyr::mutate(comparison_set = set_name, sample = sample_id,
                    z_height = unname(BR6660_DEPTHS[[sample_id]]), .before = 1)
  }))
  universe <- pair_rows |> dplyr::distinct(reference, neighbor, pair_key)
  contribution <- tidyr::crossing(sample = selected_samples, universe) |>
    dplyr::left_join(pair_rows |> dplyr::select(sample, reference, neighbor, pair_included, exclusion_reasons),
                     by = c("sample", "reference", "neighbor")) |>
    dplyr::mutate(
      comparison_set = set_name,
      z_height = unname(BR6660_DEPTHS[sample]),
      pair_included = dplyr::coalesce(pair_included, FALSE),
      exclusion_reasons = dplyr::coalesce(exclusion_reasons, "pair_absent_from_task")
    ) |>
    dplyr::arrange(z_height, reference, neighbor)
  coverage <- contribution |>
    dplyr::group_by(comparison_set, reference, neighbor, pair_key) |>
    dplyr::summarise(
      n_depths_expected = length(selected_samples),
      n_depths_eligible_and_shared = sum(pair_included),
      fraction_depths_eligible_and_shared = mean(pair_included), .groups = "drop"
    )
  stats <- auc_all |>
    dplyr::group_by(comparison_set, region, reference, neighbor, pair_key) |>
    dplyr::summarise(
      n_depths = dplyr::n_distinct(sample), mean_auc = mean(auc), median_auc = median(auc),
      variance_auc = stats::var(auc), min_auc = min(auc), max_auc = max(auc), .groups = "drop"
    )
  top_frequency <- dplyr::bind_rows(lapply(set_keys, function(key) {
    sample_id <- sub("^[^:]+::", "", key)
    results[[key]]$top |> dplyr::mutate(sample = sample_id)
  })) |>
    dplyr::count(reference, neighbor, pair_key, name = "n_depths_top_pair") |>
    dplyr::mutate(fraction_depths_top_pair = n_depths_top_pair / length(selected_samples), comparison_set = set_name)
  safe_write_table(auc_all, file.path(combined_dir, "auc_samples_all_depths.parquet"))
  safe_write_table(auc_all, file.path(combined_dir, "auc_samples_all_depths.csv.gz"))
  safe_write_table(coverage, file.path(combined_dir, "pair_depth_coverage.csv"))
  safe_write_table(stats, file.path(combined_dir, "auc_summary_by_region_and_pair.csv"))
  safe_write_table(top_frequency, file.path(combined_dir, "top_variable_pair_frequency.csv"))
  safe_write_table(contribution, file.path(combined_dir, "depth_contribution_manifest.csv"))
  safe_write_table(task_manifest |> dplyr::filter(comparison_set == set_name),
                   file.path(combined_dir, "task_manifest.csv"))
}
if (any(task_manifest$status == "failed")) quit(status = 1L)
cat("All requested across-region tasks completed.\n")
