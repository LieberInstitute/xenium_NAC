#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(arrow)
  library(crawdad)
  library(data.table)
  library(digest)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
})

# CRAWDAD within-depth analysis for Module 01 ROI assignments.

CELLTYPE_UNIVERSE_20 <- c(
  "Astro_A", "Astro_B", "Astrocyte_Oligo", "CHAT", "D1_Island_A",
  "D1_Island_B", "DRD1_MSN", "DRD2_MSN", "Ependymal", "Excitatory",
  "Fibroblast_A", "Fibroblast_B", "Inh_PVALB", "Inh_SST", "Microglia_A",
  "Microglia_B", "Microglia_Oligo", "MSN_Oligo", "OPC", "WM"
)

THREAD_ENV <- c(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)
do.call(Sys.setenv, as.list(THREAD_ENV))

usage_text <- function() {
  cat(
"Usage:
  Rscript 02_crawdad_within_depth.R [options]

Required:
  --input-assignment-parquet PATH   Module 01 parquet, csv, or csv.gz table

Region selection (default: lateral, dorsomedial, ventromedial):
  --regions NAME [NAME ...]         Atomic ROIs analyzed separately
  --combine-rois NAME=ROI,ROI       Named ROI union; option may be repeated
  --whole-tissue                    Analyze every retained tissue cell

Common options:
  --donor NAME                      Donor filter [Br6660]
  --samples all|SAMPLE [...]        Samples to analyze [all]
  --scale-range MIN MAX             Shuffle-scale bounds [200 1000]
  --scale-interval NUM              Shuffle-scale interval [100]
  --scales-um NUM [...]             Explicit scales (backward-compatible)
  --neighborhood-distances-um NUM   One neighborhood radius per run [50]
  --permutations INT                Permutations [3]
  --seed INT                        Requested CRAWDAD seed [1]
  --threshold-mode global_fixed_universe|local
                                    Bonferroni policy [global_fixed_universe]
  --alpha NUM                       Family-wise alpha [0.05]
  --min-reference-cells INT         Reference eligibility cutoff [20]
  --min-total-cells INT             Skip tasks below this count [100]
  --include-celltypes TYPE [...]    Exploratory global inclusion filter
  --exclude-celltypes TYPE [...]    Exploratory global exclusion filter
  --highlight-neighbor-celltypes TYPE [...]
  --celltype-universe-file PATH     One ordered cell type per line
  --drop-celltypes-below-count INT  Exploratory task-local removal [off]
  --total-cpus INT                  Explicit allocated CPUs
  --n-jobs INT                      Concurrent sample-region tasks
  --output-dir PATH                 Processed-data output directory
  --plot-dir PATH                   Plot output directory
  --plot-formats png [pdf]          Plot formats [png]
  --z-score-limit NUM               Shared heatmap/dot-plot limit [10]
  --resume / --no-resume            Resume valid completed tasks [resume]
  --overwrite                       Replace outputs for selected tasks
  --fail-fast                       Stop after a task failure
  --pairwise-plots                  Also plot every reference-neighbor pair
  --help                            Show this message

Column overrides:
  --cell-id-column NAME             [cell_id]
  --sample-column NAME              [sample]
  --celltype-column NAME            [cell_type]
  --x-column NAME                   [x_aligned]
  --y-column NAME                   [y_aligned]
  --depth-column NAME               [z_height]
  --roi-column NAME                 [primary_roi]
  --donor-column NAME               [donor]

Example: all six requested comparisons
  Rscript 02_crawdad_within_depth.R \\
    --input-assignment-parquet ../processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois/cell_assignments/all_cells_roi_assignment.parquet \\
    --donor Br6660 --samples all \\
    --regions lateral dorsomedial ventromedial outside_global_roi \\
    --combine-rois medial=dorsomedial,ventromedial --whole-tissue \\
    --scale-range 200 1000 --scale-interval 100 \\
    --neighborhood-distances-um 50 \\
    --permutations 3 --seed 1 --total-cpus 10 --n-jobs 10 \\
    --output-dir ../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --plot-dir ../plots/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --resume
", sep = "")
}

default_options <- function() {
  list(
    input_assignment_parquet = NULL,
    donor = "Br6660",
    samples = "all",
    regions = character(),
    combine_rois = character(),
    whole_tissue = FALSE,
    scale_range = c(200, 1000),
    scale_interval = 100,
    scales_um = NULL,
    neighborhood_distances_um = 50,
    permutations = 3L,
    seed = 1L,
    grid_shape = "square",
    threshold_mode = "global_fixed_universe",
    alpha = 0.05,
    min_reference_cells = 20L,
    min_total_cells = 100L,
    include_celltypes = character(),
    exclude_celltypes = character(),
    highlight_neighbor_celltypes = character(),
    celltype_universe_file = NULL,
    drop_celltypes_below_count = NA_integer_,
    total_cpus = NA_integer_,
    n_jobs = NA_integer_,
    output_dir = NULL,
    plot_dir = NULL,
    plot_formats = "png",
    z_score_limit = 10,
    resume = TRUE,
    overwrite = FALSE,
    fail_fast = FALSE,
    pairwise_plots = FALSE,
    cell_id_column = "cell_id",
    sample_column = "sample",
    celltype_column = "cell_type",
    x_column = "x_aligned",
    y_column = "y_aligned",
    depth_column = "z_height",
    roi_column = "primary_roi",
    donor_column = "donor",
    self_test = FALSE
  )
}

parse_cli <- function(argv) {
  opts <- default_options()
  supplied <- list(
    scale_range = "--scale-range" %in% argv,
    scale_interval = "--scale-interval" %in% argv,
    scales_um = "--scales-um" %in% argv
  )
  boolean_keys <- c(
    "whole_tissue", "overwrite", "fail_fast", "pairwise_plots", "self_test"
  )
  multi_keys <- c(
    "samples", "regions", "combine_rois", "scale_range", "scales_um",
    "neighborhood_distances_um", "highlight_neighbor_celltypes",
    "include_celltypes", "exclude_celltypes", "plot_formats"
  )
  integer_keys <- c(
    "permutations", "seed", "min_reference_cells", "min_total_cells",
    "drop_celltypes_below_count", "total_cpus", "n_jobs"
  )
  numeric_keys <- c("alpha", "z_score_limit", "scale_interval")
  if ("--help" %in% argv || "-h" %in% argv) {
    usage_text()
    quit(save = "no", status = 0)
  }
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!startsWith(token, "--")) stop("Unexpected positional argument: ", token)
    if (token == "--resume") {
      opts$resume <- TRUE
      i <- i + 1L
      next
    }
    if (token == "--no-resume") {
      opts$resume <- FALSE
      i <- i + 1L
      next
    }
    key <- gsub("-", "_", substring(token, 3), fixed = TRUE)
    if (key == "min_task_cells") key <- "min_total_cells"
    if (!key %in% names(opts)) stop("Unknown option: ", token)
    if (key %in% boolean_keys) {
      opts[[key]] <- TRUE
      i <- i + 1L
      next
    }
    j <- i + 1L
    if (j > length(argv) || startsWith(argv[[j]], "--")) {
      stop("Option requires a value: ", token)
    }
    if (key %in% multi_keys) {
      values <- character()
      while (j <= length(argv) && !startsWith(argv[[j]], "--")) {
        values <- c(values, argv[[j]])
        j <- j + 1L
      }
      if (key == "combine_rois") {
        opts[[key]] <- c(opts[[key]], values)
      } else {
        opts[[key]] <- values
      }
      i <- j
    } else {
      opts[[key]] <- argv[[j]]
      i <- j + 1L
    }
  }
  for (key in integer_keys) {
    if (!is.null(opts[[key]]) && !is.na(opts[[key]])) {
      opts[[key]] <- suppressWarnings(as.integer(opts[[key]]))
      if (is.na(opts[[key]])) stop("--", gsub("_", "-", key), " must be an integer")
    }
  }
  for (key in numeric_keys) {
    opts[[key]] <- suppressWarnings(as.numeric(opts[[key]]))
    if (is.na(opts[[key]])) stop("--", gsub("_", "-", key), " must be numeric")
  }
  opts$scale_range <- as.numeric(opts$scale_range)
  if (length(opts$scale_range) != 2L || anyNA(opts$scale_range) ||
      opts$scale_range[[1]] <= 0 || opts$scale_range[[2]] < opts$scale_range[[1]]) {
    stop("--scale-range requires two positive values: MIN MAX")
  }
  if (!is.finite(opts$scale_interval) || opts$scale_interval <= 0) {
    stop("--scale-interval must be a positive number")
  }
  if (supplied$scales_um && (supplied$scale_range || supplied$scale_interval)) {
    stop("Do not combine --scales-um with --scale-range or --scale-interval")
  }
  if (supplied$scales_um) {
    opts$scale_input_mode <- "explicit_vector"
    opts$scales_um <- sort(unique(as.numeric(opts$scales_um)))
  } else {
    opts$scale_input_mode <- "range_interval"
    opts$scales_um <- seq(
      opts$scale_range[[1]], opts$scale_range[[2]], by = opts$scale_interval
    )
  }
  opts$scale_range_um <- opts$scale_range
  opts$scale_interval_um <- opts$scale_interval
  opts$scales_um_resolved <- opts$scales_um
  opts$neighborhood_distances_um <- as.numeric(opts$neighborhood_distances_um)
  if (!length(opts$scales_um) || any(!is.finite(opts$scales_um)) ||
      any(opts$scales_um <= 0)) {
    stop("Resolved scales must be finite, positive, and nonempty")
  }
  if (anyNA(opts$neighborhood_distances_um) ||
      any(opts$neighborhood_distances_um <= 0)) {
    stop("--neighborhood-distances-um requires positive numbers")
  }
  if (length(opts$neighborhood_distances_um) != 1L) {
    stop(
      "Run one --neighborhood-distances-um value per invocation so every ",
      "distance has an independent, identical neighdist_<D>/ output tree"
    )
  }
  if (identical(opts$threshold_mode, "global")) {
    opts$threshold_mode <- "global_fixed_universe"
  }
  opts
}

slugify <- function(x) {
  ans <- gsub("[^A-Za-z0-9._-]+", "_", x)
  ans <- gsub("^_+|_+$", "", ans)
  ifelse(nchar(ans) == 0L, "unnamed", ans)
}

atomic_write_csv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) stop("Could not finalize ", path)
}

atomic_write_csv_gz <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid(), ".gz")
  readr::write_csv(x, tmp, na = "")
  if (!file.rename(tmp, path)) stop("Could not finalize ", path)
}

atomic_write_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  jsonlite::write_json(
    x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", na = "null"
  )
  if (!file.rename(tmp, path)) stop("Could not finalize ", path)
}

atomic_write_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp.", Sys.getpid())
  saveRDS(x, tmp, compress = "xz")
  if (!file.rename(tmp, path)) stop("Could not finalize ", path)
}

read_assignment <- function(path) {
  lower <- tolower(path)
  if (grepl("\\.parquet$", lower)) {
    as.data.frame(arrow::read_parquet(path))
  } else if (grepl("\\.csv(\\.gz)?$", lower)) {
    as.data.frame(readr::read_csv(path, show_col_types = FALSE, progress = FALSE))
  } else {
    stop("Input must end in .parquet, .csv, or .csv.gz")
  }
}

canonicalize_assignment <- function(raw, opts) {
  mapping <- c(
    cell_id = opts$cell_id_column,
    sample = opts$sample_column,
    cell_type = opts$celltype_column,
    x_aligned = opts$x_column,
    y_aligned = opts$y_column,
    z_height = opts$depth_column,
    primary_roi = opts$roi_column
  )
  if (!is.null(opts$donor_column) && opts$donor_column %in% names(raw)) {
    mapping <- c(mapping, donor = opts$donor_column)
  }
  missing <- setdiff(unname(mapping), names(raw))
  if (length(missing)) stop("Missing input columns: ", paste(missing, collapse = ", "))
  dat <- raw[, unname(mapping), drop = FALSE]
  names(dat) <- names(mapping)
  dat
}

validate_assignment <- function(dat) {
  character_cols <- c("cell_id", "sample", "cell_type", "primary_roi")
  for (nm in character_cols) {
    dat[[nm]] <- as.character(dat[[nm]])
    if (anyNA(dat[[nm]]) || any(!nzchar(dat[[nm]]))) {
      stop("Missing or empty values in ", nm)
    }
  }
  for (nm in c("x_aligned", "y_aligned", "z_height")) {
    dat[[nm]] <- as.numeric(dat[[nm]])
    if (any(!is.finite(dat[[nm]]))) stop("Non-finite values in ", nm)
  }
  key <- paste(dat$sample, dat$cell_id, sep = "\r")
  if (anyDuplicated(key)) stop("Duplicate cell_id values within sample")
  depth_n <- tapply(dat$z_height, dat$sample, function(x) length(unique(x)))
  if (any(depth_n != 1L)) stop("Each sample must have exactly one z_height")
  dat
}

read_celltype_universe <- function(path) {
  if (is.null(path)) return(CELLTYPE_UNIVERSE_20)
  x <- trimws(readLines(path, warn = FALSE))
  x <- x[nzchar(x) & !startsWith(x, "#")]
  if (anyDuplicated(x)) stop("Cell-type universe contains duplicates")
  x
}

parse_region_definitions <- function(opts, available_rois) {
  defs <- list()
  for (region in unique(opts$regions)) {
    if (!region %in% available_rois) stop("Unknown ROI requested: ", region)
    defs[[region]] <- region
  }
  for (spec in opts$combine_rois) {
    bits <- strsplit(spec, "=", fixed = TRUE)[[1]]
    if (length(bits) != 2L || !nzchar(bits[[1]]) || !nzchar(bits[[2]])) {
      stop("Invalid --combine-rois value: ", spec)
    }
    components <- unique(trimws(strsplit(bits[[2]], ",", fixed = TRUE)[[1]]))
    if (any(!components %in% available_rois)) {
      stop("Unknown component ROI in ", spec)
    }
    name <- trimws(bits[[1]])
    if (!is.null(defs[[name]]) &&
        !identical(sort(defs[[name]]), sort(components))) {
      stop("Conflicting definitions for analysis region: ", name)
    }
    defs[[name]] <- components
  }
  if (opts$whole_tissue) defs[["whole_tissue"]] <- available_rois
  if (!length(defs)) {
    defaults <- c("lateral", "dorsomedial", "ventromedial")
    if (any(!defaults %in% available_rois)) {
      stop("Default primary ROI labels are not all present in the input")
    }
    defs <- stats::setNames(as.list(defaults), defaults)
  }
  signatures <- vapply(defs, function(x) paste(sort(unique(x)), collapse = ","), "")
  defs[!duplicated(signatures)]
}

resolve_cpus <- function(opts) {
  candidates <- c(
    opts$total_cpus,
    suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", ""))),
    suppressWarnings(as.integer(Sys.getenv("SLURM_CPUS_ON_NODE", ""))),
    parallel::detectCores(logical = FALSE)
  )
  total <- candidates[which(!is.na(candidates) & candidates >= 1L)[1]]
  total <- as.integer(total)
  jobs <- if (is.na(opts$n_jobs)) total else opts$n_jobs
  if (jobs < 1L || jobs > total) {
    stop("--n-jobs must satisfy 1 <= n-jobs <= total-cpus (", total, ")")
  }
  list(total_cpus = total, n_jobs = as.integer(jobs))
}

construct_tasks <- function(dat, opts, region_defs, input_identity,
                            crawdad_version, universe) {
  sample_depth <- unique(dat[c("sample", "z_height")])
  sample_depth <- sample_depth[order(sample_depth$z_height, sample_depth$sample), ]
  if (!identical(opts$samples, "all")) {
    absent <- setdiff(opts$samples, sample_depth$sample)
    if (length(absent)) stop("Unknown samples: ", paste(absent, collapse = ", "))
    sample_depth <- sample_depth[match(unique(opts$samples), sample_depth$sample), ]
    sample_depth <- sample_depth[order(sample_depth$z_height), ]
  }
  rows <- lapply(seq_len(nrow(sample_depth)), function(i) {
    lapply(names(region_defs), function(region) {
      components <- region_defs[[region]]
      identity <- list(
        input_sha256 = input_identity$sha256,
        sample = sample_depth$sample[[i]],
        analysis_region = region,
        component_rois = sort(components),
        cell_id_column = opts$cell_id_column,
        sample_column = opts$sample_column,
        celltype_column = opts$celltype_column,
        x_column = opts$x_column,
        y_column = opts$y_column,
        depth_column = opts$depth_column,
        roi_column = opts$roi_column,
        celltype_universe = universe,
        neighborhood_distances_um = sort(unique(opts$neighborhood_distances_um)),
        scale_input_mode = opts$scale_input_mode,
        scale_range_um = opts$scale_range_um,
        scale_interval_um = opts$scale_interval_um,
        scales_um_resolved = opts$scales_um_resolved,
        permutations = opts$permutations,
        seed = opts$seed,
        remove_dups = TRUE,
        alpha = opts$alpha,
        threshold_mode = opts$threshold_mode,
        min_reference_cells = opts$min_reference_cells,
        min_total_cells = opts$min_total_cells,
        grid_shape = opts$grid_shape,
        include_celltypes = sort(opts$include_celltypes),
        exclude_celltypes = sort(opts$exclude_celltypes),
        drop_celltypes_below_count = opts$drop_celltypes_below_count,
        crawdad_version = crawdad_version,
        plot_formats = sort(opts$plot_formats),
        z_score_limit = opts$z_score_limit,
        pairwise_plots = opts$pairwise_plots
      )
      identity_json <- as.character(jsonlite::toJSON(
        identity, auto_unbox = TRUE, null = "null", na = "null",
        digits = NA, pretty = FALSE
      ))
      identity_hash <- digest::digest(identity_json, algo = "sha256")
      data.frame(
        task_id = paste0("task_", substr(identity_hash, 1, 16)),
        task_identity_hash = identity_hash,
        task_identity_json = identity_json,
        sample = sample_depth$sample[[i]],
        depth = sample_depth$z_height[[i]],
        analysis_region = region,
        component_rois = paste(components, collapse = ","),
        stringsAsFactors = FALSE
      )
    })
  })
  dplyr::bind_rows(unlist(rows, recursive = FALSE))
}

bonferroni_z <- function(alpha, n_tests) {
  round(stats::qnorm((alpha / n_tests) / 2, lower.tail = FALSE), 2)
}

summarize_trends <- function(dat_perm) {
  dat_perm |>
    dplyr::group_by(
      donor, sample, depth, analysis_region, component_rois,
      neighborhood_distance_um,
      reference, neighbor, scale_um
    ) |>
    dplyr::summarise(
      mean_z = mean(z_score, na.rm = TRUE),
      sd_z = stats::sd(z_score, na.rm = TRUE),
      n_permutations_observed = dplyr::n(),
      n_permutation_values = dplyr::n(),
      min_z = min(z_score, na.rm = TRUE),
      max_z = max(z_score, na.rm = TRUE),
      threshold_used = dplyr::first(threshold_used),
      n_reference_cells = dplyr::first(n_reference_cells),
      n_neighbor_cells = dplyr::first(n_neighbor_cells),
      reference_eligible = dplyr::first(reference_eligible),
      .groups = "drop"
    ) |>
    dplyr::mutate(
      primary_significant = dplyr::if_else(
        reference_eligible,
        is.finite(mean_z) & abs(mean_z) >= threshold_used,
        NA
      ),
      significant_enrichment = dplyr::if_else(
        reference_eligible, mean_z >= threshold_used, NA
      ),
      significant_depletion = dplyr::if_else(
        reference_eligible, mean_z <= -threshold_used, NA
      ),
      significant_either = primary_significant,
      ineligibility_reason = dplyr::if_else(
        reference_eligible, NA_character_,
        "reference_cell_count_below_minimum"
      )
    )
}

first_significant_table <- function(summary_dat, threshold) {
  keys <- summary_dat |>
    dplyr::filter(reference_eligible) |>
    dplyr::distinct(
      donor, sample, depth, analysis_region, neighborhood_distance_um,
      reference, neighbor, reference_eligible
    )
  hits <- summary_dat |>
    dplyr::filter(reference_eligible, primary_significant %in% TRUE) |>
    dplyr::arrange(scale_um) |>
    dplyr::group_by(
      donor, sample, depth, analysis_region, neighborhood_distance_um,
      reference, neighbor
    ) |>
    dplyr::slice_head(n = 1L) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      donor, sample, depth, analysis_region, neighborhood_distance_um,
      reference, neighbor,
      first_significant_scale_um = scale_um,
      z_at_first_significant_scale = mean_z,
      relationship_at_first_significant_scale =
        ifelse(mean_z > 0, "enrichment", "depletion")
    )
  keys |>
    dplyr::left_join(
      hits,
      by = c(
        "donor", "sample", "depth", "analysis_region", "neighborhood_distance_um",
        "reference", "neighbor"
      )
    ) |>
    dplyr::mutate(
      ever_significant = !is.na(first_significant_scale_um),
      primary_significant = ever_significant,
      threshold_used = threshold
    )
}

save_plot_formats <- function(plot, base_path, formats, width, height) {
  dir.create(dirname(base_path), recursive = TRUE, showWarnings = FALSE)
  if ("png" %in% formats) {
    ggplot2::ggsave(
      paste0(base_path, ".png"), plot, width = width, height = height,
      dpi = 300, bg = "white"
    )
  }
  if ("pdf" %in% formats) {
    ggplot2::ggsave(
      paste0(base_path, ".pdf"), plot, width = width, height = height,
      device = grDevices::cairo_pdf, bg = "white"
    )
  }
}

safe_plot_output <- function(plot_factory, base_path, formats, width, height,
                             plot_name, plot_type, implementation,
                             official_function, threshold, error_log,
                             reference = NA_character_,
                             neighbor = NA_character_) {
  output_files <- paste0(
    base_path, ".", formats[formats %in% c("png", "pdf")]
  )
  error_message <- NA_character_
  status <- "completed"
  tryCatch(
    {
      plot <- plot_factory()
      save_plot_formats(plot, base_path, formats, width, height)
    },
    error = function(e) {
      status <<- "failed"
      error_message <<- conditionMessage(e)
      cat(
        format(Sys.time()), plot_name, error_message, "\n",
        file = error_log, append = TRUE
      )
    }
  )
  data.frame(
    plot_name = plot_name,
    plot_type = plot_type,
    implementation = implementation,
    official_function = official_function,
    reference = reference,
    neighbor = neighbor,
    threshold_used = threshold,
    status = status,
    output_file = paste(output_files, collapse = ";"),
    error_message = error_message,
    stringsAsFactors = FALSE
  )
}

make_no_significant_placeholder <- function(threshold, title) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text", x = 0, y = 0,
      label = paste0("No relationships passed |Z| >= ", threshold, "."),
      size = 5
    ) +
    ggplot2::xlim(-1, 1) +
    ggplot2::ylim(-1, 1) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void()
}

make_heatmap_custom_wrapper <- function(dat, threshold, z_limit, title) {
  plot_dat <- dat
  plot_dat$scale <- factor(
    as.character(plot_dat$scale),
    levels = as.character(sort(unique(plot_dat$scale))),
    ordered = TRUE
  )
  plot_dat$Z <- pmax(-z_limit, pmin(z_limit, plot_dat$Z))
  plot_dat$Z[
    is.finite(plot_dat$Z) &
      plot_dat$Z < threshold & plot_dat$Z > -threshold
  ] <- 0
  ggplot2::ggplot(
    plot_dat, ggplot2::aes(x = scale, y = neighbor, fill = Z)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~reference, ncol = 4) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-z_limit, z_limit)
    ) +
    ggplot2::labs(
      x = "Shuffle scale (µm)", y = "Neighbor", fill = "Z score",
      title = title
    ) +
    ggplot2::theme_classic()
}

safe_heatmap_output <- function(dat, base_path, formats, width, height,
                                plot_name, threshold, z_limit, title,
                                error_log, reference = NA_character_) {
  official_error <- NA_character_
  official_plot <- tryCatch(
    crawdad::vizTrends.heatmap(
      dat = dat,
      withPerms = FALSE,
      zSigThresh = threshold,
      zScoreLimit = z_limit,
      title = title
    ),
    error = function(e) {
      official_error <<- conditionMessage(e)
      NULL
    }
  )
  if (!is.null(official_plot)) {
    return(safe_plot_output(
      function() official_plot, paste0(base_path, "__official_crawdad"),
      formats, width, height, plot_name, "multi_scale_heatmap",
      "official_crawdad", "crawdad::vizTrends.heatmap", threshold,
      error_log, reference = reference
    ))
  }
  fallback_base <- paste0(base_path, "__custom_wrapper")
  fallback <- safe_plot_output(
    function() make_heatmap_custom_wrapper(
      dat, threshold, z_limit, title
    ),
    fallback_base, formats, width, height, plot_name,
    "multi_scale_heatmap", "custom_wrapper",
    "crawdad::vizTrends.heatmap attempted; dependency-incompatible",
    threshold, error_log, reference = reference
  )
  if (fallback$status == "completed") {
    fallback$status <- "completed_with_custom_fallback"
    fallback$error_message <- official_error
    cat(
      format(Sys.time()), plot_name,
      "official CRAWDAD heatmap failed; custom fallback saved:",
      official_error, "\n", file = error_log, append = TRUE
    )
  }
  fallback
}

task_required_files <- function(task_dir, distances) {
  if (length(distances) != 1L) stop("Exactly one neighborhood distance is required")
  c(
    file.path(task_dir, "status.json"),
    file.path(task_dir, "task_identity.json"),
    file.path(task_dir, "cells_sf.rds"),
    file.path(task_dir, "shuffle_list.rds"),
    file.path(task_dir, "celltype_counts.csv"),
    file.path(task_dir, "reference_eligibility.csv"),
    file.path(task_dir, "plot_manifest.csv"),
    file.path(task_dir, "find_trends.rds"),
    file.path(task_dir, "permutation_results.rds"),
    file.path(task_dir, "thresholds.csv"),
    file.path(task_dir, "trends_permutation.csv.gz"),
    file.path(task_dir, "trends_summary.csv.gz"),
    file.path(task_dir, "first_significant_relationships.csv.gz"),
    file.path(task_dir, "neighborhood_status.json")
  )
}

is_task_complete <- function(task_dir, distances, expected_identity_hash) {
  files <- task_required_files(task_dir, distances)
  if (!all(file.exists(files)) || any(file.info(files)$size <= 0)) return(FALSE)
  status <- tryCatch(jsonlite::read_json(file.path(task_dir, "status.json")), error = identity)
  !inherits(status, "error") &&
    identical(status$status, "completed") &&
    identical(status$task_identity_hash, expected_identity_hash)
}

task_subset <- function(dat, task_row) {
  components <- strsplit(task_row$component_rois, ",", fixed = TRUE)[[1]]
  ans <- dat[
    dat$sample == task_row$sample & dat$primary_roi %in% components,
    , drop = FALSE
  ]
  if (anyDuplicated(ans$cell_id)) stop("Constructed task contains duplicate cells")
  ans
}

run_one_task <- function(task_row, dat, opts, universe, thresholds, paths) {
  task_start <- Sys.time()
  sample_slug <- slugify(task_row$sample)
  region_slug <- slugify(task_row$analysis_region)
  region_output_root <- file.path(paths$output_root, region_slug)
  region_plot_root <- file.path(paths$plot_root, region_slug)
  final_dir <- file.path(region_output_root, "tasks", sample_slug)
  staging_dir <- paste0(final_dir, ".inprogress")
  plot_task_dir <- file.path(region_plot_root, sample_slug)
  log_root <- file.path(region_output_root, "logs")
  stdout_log <- file.path(log_root, paste0(task_row$task_id, ".stdout.log"))
  stderr_log <- file.path(log_root, paste0(task_row$task_id, ".stderr.log"))
  dir.create(log_root, recursive = TRUE, showWarnings = FALSE)

  if (!opts$overwrite && opts$resume &&
      is_task_complete(
        final_dir, opts$neighborhood_distances_um,
        task_row$task_identity_hash
      )) {
    previous <- jsonlite::read_json(file.path(final_dir, "status.json"))
    return(data.frame(
      task_row, status = "resumed",
      n_cells = previous$n_cells %||% NA_integer_,
      n_observed_celltypes = previous$n_observed_celltypes %||% NA_integer_,
      message = "valid completed task",
      elapsed_seconds = 0, stringsAsFactors = FALSE
    ))
  }
  if (dir.exists(staging_dir)) unlink(staging_dir, recursive = TRUE, force = TRUE)
  if (opts$overwrite && dir.exists(final_dir)) {
    unlink(final_dir, recursive = TRUE, force = TRUE)
  }
  if (dir.exists(plot_task_dir)) {
    unlink(plot_task_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plot_task_dir, recursive = TRUE, showWarnings = FALSE)

  stdout_con <- file(stdout_log, open = "wt")
  stderr_con <- file(stderr_log, open = "wt")
  sink(stdout_con, type = "output")
  sink(stderr_con, type = "message")
  on.exit({
    while (sink.number(type = "message") > 2L) sink(type = "message")
    while (sink.number(type = "output") > 0L) sink(type = "output")
    close(stdout_con)
    close(stderr_con)
  }, add = TRUE)

  result <- tryCatch({
    cat("Task:", task_row$task_id, task_row$sample, task_row$analysis_region, "\n")
    task_dat <- task_subset(dat, task_row)
    if (length(opts$include_celltypes)) {
      task_dat <- task_dat[
        task_dat$cell_type %in% opts$include_celltypes, , drop = FALSE
      ]
    }
    if (length(opts$exclude_celltypes)) {
      task_dat <- task_dat[
        !task_dat$cell_type %in% opts$exclude_celltypes, , drop = FALSE
      ]
    }
    counts <- as.data.frame(table(factor(task_dat$cell_type, levels = universe)))
    names(counts) <- c("cell_type", "n_cells")
    counts$cell_type <- as.character(counts$cell_type)
    counts$present <- counts$n_cells > 0L
    counts$reference_eligible <- counts$n_cells >= opts$min_reference_cells
    counts$sample <- task_row$sample
    counts$depth <- task_row$depth
    counts$analysis_region <- task_row$analysis_region

    observed_unknown <- setdiff(unique(task_dat$cell_type), universe)
    if (length(observed_unknown)) {
      stop("Observed cell types absent from universe: ", paste(observed_unknown, collapse = ", "))
    }
    if (!is.na(opts$drop_celltypes_below_count)) {
      keep_types <- counts$cell_type[counts$n_cells >= opts$drop_celltypes_below_count]
      task_dat <- task_dat[task_dat$cell_type %in% keep_types, , drop = FALSE]
      counts$dropped_from_task <- !counts$cell_type %in% keep_types
    } else {
      counts$dropped_from_task <- FALSE
    }
    n_cells <- nrow(task_dat)
    n_observed <- length(unique(task_dat$cell_type))
    atomic_write_csv(counts, file.path(staging_dir, "celltype_counts.csv"))

    metadata <- list(
      task_id = task_row$task_id,
      task_identity_hash = task_row$task_identity_hash,
      sample = task_row$sample,
      depth = task_row$depth,
      analysis_region = task_row$analysis_region,
      component_rois = strsplit(task_row$component_rois, ",", fixed = TRUE)[[1]],
      n_cells = n_cells,
      n_observed_celltypes = n_observed,
      coordinate_units = "micrometers",
      scale_input_mode = opts$scale_input_mode,
      scale_range_um = opts$scale_range_um,
      scale_interval_um = opts$scale_interval_um,
      scales_um_resolved = opts$scales_um_resolved,
      neighborhood_distances_um = opts$neighborhood_distances_um,
      permutations = opts$permutations,
      random_seed = opts$seed,
      crawdad_version = as.character(utils::packageVersion("crawdad")),
      stdout_log = stdout_log,
      warnings_and_stderr_log = stderr_log,
      outer_parallelism = paths$n_jobs,
      crawdad_ncores = 1L
    )
    atomic_write_json(metadata, file.path(staging_dir, "task_metadata.json"))
    atomic_write_csv(
      counts |>
        dplyr::transmute(
          sample, depth, analysis_region, cell_type,
          n_reference_cells = n_cells,
          reference_eligible,
          ineligibility_reason = dplyr::if_else(
            reference_eligible, NA_character_,
            "reference_cell_count_below_minimum"
          )
        ),
      file.path(staging_dir, "reference_eligibility.csv")
    )
    identity_resolved <- jsonlite::fromJSON(task_row$task_identity_json)
    atomic_write_json(
      list(
        task_identity_hash = task_row$task_identity_hash,
        resolved_identity = identity_resolved
      ),
      file.path(staging_dir, "task_identity.json")
    )

    if (n_cells < opts$min_total_cells || n_observed < 2L) {
      reason <- if (n_cells < opts$min_total_cells) {
        paste0("fewer than ", opts$min_total_cells, " cells")
      } else {
        "fewer than two observed cell types"
      }
      atomic_write_json(
        list(
          status = "skipped_ineligible", reason = reason,
          task_identity_hash = task_row$task_identity_hash,
          completed_at = format(Sys.time())
        ),
        file.path(staging_dir, "status.json")
      )
      if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)
      if (!file.rename(staging_dir, final_dir)) stop("Could not finalize skipped task")
      return(list(
        status = "skipped_ineligible", n_cells = n_cells,
        n_observed = n_observed, message = reason
      ))
    }

    task_input <- task_dat |>
      dplyr::transmute(
        cell_id, sample, celltype = cell_type, x = x_aligned, y = y_aligned,
        depth = z_height, primary_roi
      )
    atomic_write_csv_gz(task_input, file.path(staging_dir, "task_input.csv.gz"))
    capture.output(sessionInfo(), file = file.path(staging_dir, "session_info.txt"))

    cells <- crawdad::toSF(
      pos = task_input[, c("x", "y")],
      cellTypes = factor(
        task_input$celltype,
        levels = universe[universe %in% unique(task_input$celltype)]
      ),
      verbose = TRUE
    )
    atomic_write_rds(cells, file.path(staging_dir, "cells_sf.rds"))
    shuffle_list <- crawdad::makeShuffledCells(
      cells = cells,
      scales = opts$scales_um,
      perms = opts$permutations,
      ncores = 1,
      seed = opts$seed,
      square = identical(opts$grid_shape, "square"),
      verbose = TRUE
    )
    atomic_write_rds(shuffle_list, file.path(staging_dir, "shuffle_list.rds"))
    plot_manifest_rows <- list()
    plot_error_log <- file.path(staging_dir, "plot_errors.log")
    writeLines(character(), plot_error_log)

    for (distance in opts$neighborhood_distances_um) {
      nd_dir <- staging_dir
      nd_plot_dir <- plot_task_dir
      dir.create(nd_dir, recursive = TRUE, showWarnings = FALSE)
      dir.create(nd_plot_dir, recursive = TRUE, showWarnings = FALSE)

      results <- crawdad::findTrends(
        cells = cells,
        neighDist = distance,
        ncores = 1,
        shuffleList = shuffle_list,
        verbose = TRUE,
        removeDups = TRUE,
        returnMeans = FALSE
      )
      atomic_write_rds(results, file.path(nd_dir, "find_trends.rds"))

      dat_perm <- as.data.frame(crawdad::meltResultsList(results, withPerms = TRUE))
      names(dat_perm)[names(dat_perm) == "Z"] <- "z_score"
      names(dat_perm)[names(dat_perm) == "scale"] <- "scale_um"
      dat_perm$neighbor <- as.character(dat_perm$neighbor)
      dat_perm$reference <- as.character(dat_perm$reference)
      dat_perm$permutation <- as.integer(dat_perm$perm)
      dat_perm$perm <- NULL
      dat_perm$id <- NULL
      local_threshold <- crawdad::correctZBonferroni(dat_perm, pSigThresh = opts$alpha)
      threshold <- if (opts$threshold_mode == "global_fixed_universe") {
        thresholds$global_z_threshold[[1]]
      } else {
        local_threshold
      }
      dat_perm <- dat_perm |>
        dplyr::mutate(
          task_id = task_row$task_id,
          donor = opts$donor,
          sample = task_row$sample,
          id = task_row$sample,
          depth = task_row$depth,
          analysis_region = task_row$analysis_region,
          region = task_row$analysis_region,
          component_rois = task_row$component_rois,
          neighborhood_distance_um = distance,
          n_task_cells = n_cells,
          n_reference_cells =
            counts$n_cells[match(reference, counts$cell_type)],
          n_neighbor_cells =
            counts$n_cells[match(neighbor, counts$cell_type)],
          reference_eligible =
            counts$reference_eligible[match(reference, counts$cell_type)],
          global_z_threshold = thresholds$global_z_threshold[[1]],
          local_z_threshold = local_threshold,
          threshold_mode = opts$threshold_mode,
          threshold_used = threshold,
          scale = scale_um,
          Z = z_score,
          z = z_score,
          significant_positive = dplyr::if_else(
            reference_eligible, z_score >= threshold, NA
          ),
          significant_negative = dplyr::if_else(
            reference_eligible, z_score <= -threshold, NA
          ),
          ineligibility_reason = dplyr::if_else(
            reference_eligible, NA_character_,
            "reference_cell_count_below_minimum"
          )
        ) |>
        dplyr::select(
          task_id, donor, id, sample, depth, region, analysis_region,
          component_rois,
          neighborhood_distance_um,
          reference, neighbor, scale, scale_um, permutation, Z, z_score,
          z,
          n_task_cells, n_reference_cells, n_neighbor_cells,
          reference_eligible, ineligibility_reason,
          global_z_threshold, local_z_threshold, threshold_mode,
          threshold_used, significant_positive, significant_negative
        ) |>
        dplyr::arrange(reference, neighbor, scale_um, permutation)

      expected_rows <- length(unique(dat_perm$scale_um)) *
        length(unique(dat_perm$reference)) *
        length(unique(dat_perm$neighbor)) * opts$permutations
      if (nrow(dat_perm) != expected_rows) {
        stop("Permutation output is incomplete: expected ", expected_rows,
             " rows, observed ", nrow(dat_perm))
      }
      summary_dat <- summarize_trends(dat_perm)
      first_dat <- first_significant_table(summary_dat, threshold)
      atomic_write_csv_gz(dat_perm, file.path(nd_dir, "trends_permutation.csv.gz"))
      atomic_write_csv_gz(summary_dat, file.path(nd_dir, "trends_summary.csv.gz"))
      atomic_write_csv_gz(
        first_dat, file.path(nd_dir, "first_significant_relationships.csv.gz")
      )
      atomic_write_rds(
        dat_perm, file.path(nd_dir, "permutation_results.rds")
      )
      threshold_dat <- data.frame(
        task_id = task_row$task_id,
        sample = task_row$sample,
        depth = task_row$depth,
        analysis_region = task_row$analysis_region,
        neighborhood_distance_um = distance,
        alpha = opts$alpha,
        threshold_mode = opts$threshold_mode,
        global_celltype_count = length(universe),
        global_n_directional_tests = length(universe)^2L,
        global_z_threshold = thresholds$global_z_threshold[[1]],
        local_observed_celltype_count =
          length(unique(dat_perm$reference)),
        local_n_directional_tests =
          length(unique(dat_perm$reference))^2L,
        local_z_threshold = local_threshold,
        threshold_used = threshold,
        threshold_policy_note = paste(
          "The global fixed-universe threshold is the study-specific",
          "primary policy. The local threshold is the official CRAWDAD",
          "task-specific threshold and is saved for QC only."
        ),
        stringsAsFactors = FALSE
      )
      atomic_write_csv(threshold_dat, file.path(nd_dir, "thresholds.csv"))

      official_perm <- dat_perm |>
        dplyr::transmute(
          perm = permutation, neighbor, Z, scale, reference,
          id = task_row$task_id, reference_eligible
        )
      official_mean <- summary_dat |>
        dplyr::transmute(
          neighbor, scale = scale_um, reference, Z = mean_z,
          reference_eligible
        )
      primary_perm <- official_perm |>
        dplyr::filter(reference_eligible)
      primary_mean <- official_mean |>
        dplyr::filter(reference_eligible)
      base_title <- paste(
        task_row$sample, paste0("depth ", task_row$depth, " µm"),
        task_row$analysis_region, paste0("neighDist ", distance, " µm"),
        sep = " | "
      )

      plot_manifest_rows[[length(plot_manifest_rows) + 1L]] <-
        safe_heatmap_output(
          primary_mean |>
            dplyr::select(neighbor, scale, reference, Z),
          file.path(nd_plot_dir, "heatmap_all_references"),
          opts$plot_formats, 14, 12, "heatmap_all_references",
          threshold, opts$z_score_limit, base_title, plot_error_log
        )

      has_significant <- nrow(primary_mean) > 0L && any(
        is.finite(primary_mean$Z) & abs(primary_mean$Z) >= threshold
      )
      relationship_implementation <- if (has_significant) {
        "official_crawdad"
      } else {
        "placeholder"
      }
      relationship_function <- if (has_significant) {
        "crawdad::vizRelationships"
      } else {
        NA_character_
      }
      relationship_factory <- if (has_significant) {
        function() {
          crawdad::vizRelationships(
            primary_perm |>
              dplyr::select(perm, neighbor, Z, scale, reference, id),
            zSigThresh = threshold,
            zScoreLimit = opts$z_score_limit,
            reorder = TRUE,
            symmetrical = FALSE,
            onlySignificant = FALSE,
            dotSizes = c(2, 10)
          ) + ggplot2::labs(title = base_title)
        }
      } else {
        function() make_no_significant_placeholder(threshold, base_title)
      }
      relationship_row <- safe_plot_output(
        relationship_factory,
        file.path(nd_plot_dir, "relationship_summary"),
        opts$plot_formats, 10, 9,
        "relationship_summary", "first_significant_relationship_summary",
        relationship_implementation, relationship_function,
        threshold, plot_error_log
      )
      if (!has_significant && relationship_row$status == "completed") {
        relationship_row$status <- "no_significant_relationships"
      }
      plot_manifest_rows[[length(plot_manifest_rows) + 1L]] <- relationship_row

      eligible_references <- counts$cell_type[
        counts$reference_eligible & counts$present
      ]
      for (reference_type in eligible_references) {
        title <- paste(
          base_title,
          paste0("reference: ", reference_type), sep = " | "
        )
        ref_slug <- slugify(reference_type)
        ref_mean <- primary_mean |>
          dplyr::filter(reference == reference_type) |>
          dplyr::select(neighbor, scale, reference, Z)
        ref_perm <- primary_perm |>
          dplyr::filter(reference == reference_type) |>
          dplyr::select(perm, neighbor, Z, scale, reference, id)
        plot_manifest_rows[[length(plot_manifest_rows) + 1L]] <-
          safe_heatmap_output(
            ref_mean,
            file.path(nd_plot_dir, paste0("heatmap_ref_", ref_slug)),
            opts$plot_formats, 7, 6, paste0("heatmap_ref_", ref_slug),
            threshold, opts$z_score_limit, title, plot_error_log,
            reference = reference_type
          )
        plot_manifest_rows[[length(plot_manifest_rows) + 1L]] <-
          safe_plot_output(
            function() crawdad::vizTrends(
              ref_perm,
              id = "neighbor",
              lines = TRUE,
              points = TRUE,
              withPerms = TRUE,
              facet = FALSE,
              zSigThresh = threshold,
              title = title
            ),
            file.path(
              nd_plot_dir,
              paste0("trends_ref_", ref_slug, "__official_crawdad")
            ),
            opts$plot_formats, 10, 7,
            paste0("trends_ref_", ref_slug), "reference_neighbor_trends",
            "official_crawdad", "crawdad::vizTrends",
            threshold, plot_error_log, reference = reference_type
          )
        if (opts$pairwise_plots) {
          for (neighbor_type in unique(ref_perm$neighbor)) {
            pair_dat <- ref_perm |>
              dplyr::filter(neighbor == neighbor_type)
            pair_slug <- paste0(
              "ref_", ref_slug, "__neighbor_", slugify(neighbor_type)
            )
            plot_manifest_rows[[length(plot_manifest_rows) + 1L]] <-
              safe_plot_output(
                function() crawdad::vizTrends(
                  pair_dat,
                  lines = TRUE,
                  points = TRUE,
                  withPerms = TRUE,
                  facet = FALSE,
                  zSigThresh = threshold,
                  title = paste(
                    title, paste0("neighbor: ", neighbor_type), sep = " | "
                  )
                ),
                file.path(nd_plot_dir, "pairwise", pair_slug),
                opts$plot_formats, 6, 4,
                pair_slug, "reference_neighbor_pair_trend",
                "official_crawdad", "crawdad::vizTrends",
                threshold, plot_error_log,
                reference = reference_type, neighbor = neighbor_type
              )
          }
        }
      }
      relationship_plot_status <- relationship_row$status
      atomic_write_json(
        list(
          status = "completed",
          neighborhood_distance_um = distance,
          n_permutation_rows = nrow(dat_perm),
          threshold_used = threshold,
          relationship_plot_status = relationship_plot_status,
          completed_at = format(Sys.time())
        ),
        file.path(nd_dir, "neighborhood_status.json")
      )
    }
    plot_manifest <- dplyr::bind_rows(plot_manifest_rows)
    plot_manifest$neighborhood_distance_um <- opts$neighborhood_distances_um[[1]]
    atomic_write_csv(
      plot_manifest, file.path(staging_dir, "plot_manifest.csv")
    )

    atomic_write_json(
      list(
        status = "completed",
        task_id = task_row$task_id,
        task_identity_hash = task_row$task_identity_hash,
        completed_at = format(Sys.time()),
        runtime_seconds = as.numeric(
          difftime(Sys.time(), task_start, units = "secs")
        ),
        n_cells = n_cells,
        n_observed_celltypes = n_observed
      ),
      file.path(staging_dir, "status.json")
    )
    if (dir.exists(final_dir)) unlink(final_dir, recursive = TRUE, force = TRUE)
    dir.create(dirname(final_dir), recursive = TRUE, showWarnings = FALSE)
    if (!file.rename(staging_dir, final_dir)) stop("Could not finalize task directory")
    list(status = "completed", n_cells = n_cells, n_observed = n_observed,
         message = "")
  }, error = function(e) {
    atomic_write_json(
      list(
        status = "failed", task_id = task_row$task_id,
        task_identity_hash = task_row$task_identity_hash,
        error = conditionMessage(e), failed_at = format(Sys.time())
      ),
      file.path(staging_dir, "status.json")
    )
    list(status = "failed", n_cells = NA_integer_, n_observed = NA_integer_,
         message = conditionMessage(e))
  })

  data.frame(
    task_row,
    status = result$status,
    n_cells = result$n_cells,
    n_observed_celltypes = result$n_observed,
    message = result$message,
    elapsed_seconds = as.numeric(difftime(Sys.time(), task_start, units = "secs")),
    stringsAsFactors = FALSE
  )
}

collect_task_outputs <- function(tasks, opts, paths) {
  trend_perm <- list()
  trend_summary <- list()
  first_sig <- list()
  counts <- list()
  eligibility <- list()
  geometry <- list()
  plot_manifests <- list()
  for (i in seq_len(nrow(tasks))) {
    task <- tasks[i, ]
    task_dir <- file.path(
      paths$output_root, slugify(task$analysis_region), "tasks",
      slugify(task$sample)
    )
    status_path <- file.path(task_dir, "status.json")
    if (!file.exists(status_path)) next
    status <- tryCatch(jsonlite::read_json(status_path), error = identity)
    if (inherits(status, "error")) next
    count_path <- file.path(task_dir, "celltype_counts.csv")
    if (file.exists(count_path)) {
      count_dat <- readr::read_csv(count_path, show_col_types = FALSE)
      counts[[length(counts) + 1L]] <- count_dat
      eligibility[[length(eligibility) + 1L]] <- count_dat |>
        dplyr::select(
          sample, depth, analysis_region, cell_type, n_cells,
          present, reference_eligible, dropped_from_task
        )
    }
    if (!identical(status$status, "completed")) next
    plot_manifest_path <- file.path(task_dir, "plot_manifest.csv")
    if (file.exists(plot_manifest_path)) {
      plot_manifests[[length(plot_manifests) + 1L]] <- readr::read_csv(
        plot_manifest_path, show_col_types = FALSE
      ) |>
        dplyr::mutate(
          task_id = task$task_id,
          sample = task$sample,
          depth = task$depth,
          analysis_region = task$analysis_region,
          .before = 1
        )
    }
    input_dat <- readr::read_csv(
      file.path(task_dir, "task_input.csv.gz"), show_col_types = FALSE,
      progress = FALSE
    )
    geometry[[length(geometry) + 1L]] <- data.frame(
      sample = task$sample,
      depth = task$depth,
      analysis_region = task$analysis_region,
      n_cells = nrow(input_dat),
      n_observed_celltypes = length(unique(input_dat$celltype)),
      x_min = min(input_dat$x),
      x_max = max(input_dat$x),
      y_min = min(input_dat$y),
      y_max = max(input_dat$y),
      x_span_um = diff(range(input_dat$x)),
      y_span_um = diff(range(input_dat$y)),
      component_rois = task$component_rois,
      stringsAsFactors = FALSE
    )
    for (distance in opts$neighborhood_distances_um) {
      nd_dir <- task_dir
      trend_perm[[length(trend_perm) + 1L]] <- readr::read_csv(
        file.path(nd_dir, "trends_permutation.csv.gz"),
        show_col_types = FALSE, progress = FALSE
      )
      trend_summary[[length(trend_summary) + 1L]] <- readr::read_csv(
        file.path(nd_dir, "trends_summary.csv.gz"),
        show_col_types = FALSE, progress = FALSE
      )
      first_sig[[length(first_sig) + 1L]] <- readr::read_csv(
        file.path(nd_dir, "first_significant_relationships.csv.gz"),
        show_col_types = FALSE, progress = FALSE
      )
    }
  }
  list(
    permutation = dplyr::bind_rows(trend_perm),
    summary = dplyr::bind_rows(trend_summary),
    first = dplyr::bind_rows(first_sig),
    counts = dplyr::bind_rows(counts),
    eligibility = dplyr::bind_rows(eligibility),
    geometry = dplyr::bind_rows(geometry),
    plot_manifest = dplyr::bind_rows(plot_manifests)
  )
}

validate_auc_ready_input <- function(dat, resolved_scales) {
  required <- c("id", "reference", "neighbor", "scale", "Z")
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop("AUC-ready table is missing CRAWDAD columns: ",
         paste(missing, collapse = ", "))
  }
  if (!nrow(dat)) return(invisible(TRUE))
  if (anyNA(dat$id) || any(!nzchar(dat$id))) {
    stop("AUC-ready table contains a missing sample id")
  }
  if (!is.numeric(dat$scale) || !is.numeric(dat$Z)) {
    stop("AUC-ready scale and Z columns must be numeric")
  }
  key <- paste(
    dat$id, dat$reference, dat$neighbor, dat$scale, dat$permutation,
    sep = "\r"
  )
  if (anyDuplicated(key)) stop("AUC-ready table contains duplicate rows")
  curve_qc <- dat |>
    dplyr::group_by(id, reference, neighbor) |>
    dplyr::summarise(
      scales_match = identical(sort(unique(scale)), sort(resolved_scales)),
      finite_mean_at_every_scale = all(vapply(
        split(Z, scale), function(x) any(is.finite(x)), logical(1)
      )),
      .groups = "drop"
    )
  if (any(!curve_qc$scales_match) ||
      any(!curve_qc$finite_mean_at_every_scale)) {
    stop("AUC-ready table contains an incomplete or non-finite scale curve")
  }
  invisible(TRUE)
}

write_combined <- function(combined, task_status, opts, paths, thresholds) {
  dir.create(paths$combined_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$qc_root, recursive = TRUE, showWarnings = FALSE)
  sort_cols <- c(
    "depth", "sample", "analysis_region", "neighborhood_distance_um",
    "reference", "neighbor", "scale_um", "permutation"
  )
  tables <- list(
    all_trends_permutation = combined$permutation,
    all_trends_summary = combined$summary,
    all_first_significant_relationships = combined$first,
    all_task_qc = task_status
  )
  for (nm in names(tables)) {
    tab <- tables[[nm]]
    by <- intersect(sort_cols, names(tab))
    if (nrow(tab) && length(by)) tab <- tab[do.call(order, tab[by]), , drop = FALSE]
    arrow::write_parquet(tab, file.path(paths$combined_root, paste0(nm, ".parquet")))
    atomic_write_csv_gz(tab, file.path(paths$combined_root, paste0(nm, ".csv.gz")))
  }
  atomic_write_csv(combined$counts, file.path(paths$qc_root, "task_celltype_counts.csv"))
  task_cell_counts <- task_status |>
    dplyr::select(
      sample, depth, analysis_region, n_cells, n_observed_celltypes, status
    )
  atomic_write_csv(task_cell_counts, file.path(paths$qc_root, "task_cell_counts.csv"))
  atomic_write_csv(
    combined$eligibility, file.path(paths$qc_root, "reference_eligibility.csv")
  )
  atomic_write_csv(
    combined$geometry, file.path(paths$qc_root, "task_geometry_summary.csv")
  )
  atomic_write_csv(task_status, file.path(paths$qc_root, "task_status.csv"))
  atomic_write_csv(
    combined$plot_manifest, file.path(paths$qc_root, "plot_manifest.csv")
  )

  auc_manifest <- list()
  if (nrow(combined$permutation)) {
    auc_groups <- combined$permutation |>
      dplyr::distinct(analysis_region, neighborhood_distance_um) |>
      dplyr::arrange(analysis_region, neighborhood_distance_um)
    for (i in seq_len(nrow(auc_groups))) {
      region <- auc_groups$analysis_region[[i]]
      distance <- auc_groups$neighborhood_distance_um[[i]]
      auc_dat <- combined$permutation |>
        dplyr::filter(
          analysis_region == region,
          neighborhood_distance_um == distance,
          reference_eligible
        ) |>
        dplyr::transmute(
          id = sample,
          reference,
          neighbor,
          scale,
          Z,
          permutation,
          task_id,
          donor,
          sample,
          depth,
          region = analysis_region,
          neighborhood_distance_um,
          n_reference_cells,
          n_neighbor_cells,
          reference_eligible
        ) |>
        dplyr::group_by(id, reference, neighbor, scale) |>
        dplyr::mutate(scale_has_finite_z = any(is.finite(Z))) |>
        dplyr::ungroup() |>
        dplyr::group_by(id, reference, neighbor) |>
        dplyr::filter(
          all(scale_has_finite_z),
          dplyr::n_distinct(scale) == length(opts$scales_um_resolved)
        ) |>
        dplyr::ungroup() |>
        dplyr::select(-scale_has_finite_z) |>
        dplyr::arrange(id, reference, neighbor, scale, permutation)
      validate_auc_ready_input(auc_dat, opts$scales_um_resolved)
      auc_dir <- paths$combined_root
      dir.create(auc_dir, recursive = TRUE, showWarnings = FALSE)
      auc_path <- file.path(auc_dir, "calculate_auc_input.parquet")
      arrow::write_parquet(auc_dat, auc_path)
      atomic_write_csv_gz(
        auc_dat,
        file.path(auc_dir, "calculate_auc_input.csv.gz")
      )
      auc_manifest[[length(auc_manifest) + 1L]] <- data.frame(
        analysis_region = region,
        neighborhood_distance_um = distance,
        input_file = auc_path,
        n_rows = nrow(auc_dat),
        n_samples = dplyr::n_distinct(auc_dat$id),
        n_directional_pairs = dplyr::n_distinct(
          paste(auc_dat$reference, auc_dat$neighbor, sep = "\r")
        ),
        required_crawdad_columns = "id,reference,neighbor,scale,Z",
        split_instruction =
          "split(dat, dat$id), then crawdad::calculateAUC(..., sharedPairs=TRUE)",
        stringsAsFactors = FALSE
      )
    }
  }
  atomic_write_csv(
    dplyr::bind_rows(auc_manifest),
    file.path(paths$combined_root, "calculate_auc_input_manifest.csv")
  )
  task_thresholds <- if (nrow(combined$permutation)) {
    combined$permutation |>
      dplyr::group_by(
        donor, sample, depth, analysis_region, neighborhood_distance_um
      ) |>
      dplyr::summarise(
        alpha = opts$alpha,
        threshold_mode = dplyr::first(threshold_mode),
        global_celltype_count = 20L,
        local_observed_celltype_count = dplyr::n_distinct(reference),
        global_n_directional_tests = 400L,
        local_n_directional_tests = local_observed_celltype_count^2L,
        global_z_threshold = dplyr::first(global_z_threshold),
        local_z_threshold = dplyr::first(local_z_threshold),
        threshold_used = dplyr::first(threshold_used),
        threshold_policy_note = paste(
          "The global fixed-universe threshold is the study-specific",
          "primary policy. The local threshold is the official CRAWDAD",
          "task-specific threshold and is saved for QC only."
        ),
        .groups = "drop"
      ) |>
      dplyr::arrange(depth, sample, analysis_region, neighborhood_distance_um)
  } else {
    data.frame(
      donor = character(), sample = character(), depth = numeric(),
      analysis_region = character(), neighborhood_distance_um = numeric(),
      alpha = numeric(), threshold_mode = character(),
      global_celltype_count = integer(),
      local_observed_celltype_count = integer(),
      global_n_directional_tests = integer(),
      local_n_directional_tests = integer(),
      global_z_threshold = numeric(), local_z_threshold = numeric(),
      threshold_used = numeric(), threshold_policy_note = character()
    )
  }
  atomic_write_csv(
    task_thresholds, file.path(paths$output_root, "significance_thresholds.csv")
  )
}

input_identity <- function(path) {
  normalized <- normalizePath(path, mustWork = TRUE)
  info <- file.info(normalized)
  list(
    path = normalized,
    size_bytes = unname(info$size),
    modification_time = format(info$mtime, tz = "UTC", usetz = TRUE),
    sha256 = digest::digest(file = normalized, algo = "sha256")
  )
}

git_commit <- function(script_dir) {
  out <- tryCatch(
    system2("git", c("-C", shQuote(script_dir), "rev-parse", "HEAD"),
            stdout = TRUE, stderr = FALSE),
    error = function(e) NA_character_
  )
  if (!length(out)) NA_character_ else out[[1]]
}

write_run_metadata <- function(opts, paths, identity, universe, region_defs, cpus,
                               crawdad_version, thresholds) {
  dir.create(paths$output_root, recursive = TRUE, showWarnings = FALSE)
  dir.create(paths$plot_root, recursive = TRUE, showWarnings = FALSE)
  atomic_write_json(
    list(
      input = identity,
      git_commit = git_commit(dirname(normalizePath(sys.frame(1)$ofile %||% "."))),
      r_version = R.version.string,
      crawdad_version = crawdad_version,
      package_versions = lapply(
        c("sf", "BiocParallel", "dplyr", "ggplot2", "arrow"),
        function(x) as.character(utils::packageVersion(x))
      ) |> setNames(c("sf", "BiocParallel", "dplyr", "ggplot2", "arrow")),
      coordinate_reference = "aligned Xenium Cartesian coordinates",
      coordinate_units = "micrometers",
      total_cpus = cpus$total_cpus,
      n_jobs = cpus$n_jobs,
      inner_crawdad_ncores = 1L,
      scale_input_mode = opts$scale_input_mode,
      scale_range_um = opts$scale_range_um,
      scale_interval_um = opts$scale_interval_um,
      scales_um_resolved = opts$scales_um_resolved,
      neighborhood_distances_um = opts$neighborhood_distances_um,
      permutations = opts$permutations,
      requested_seed = opts$seed,
      effective_permutation_seeds =
        if (opts$permutations > 1L) seq_len(opts$permutations) else opts$seed,
      seed_behavior_note =
        "CRAWDAD 1.0.1 uses permutation index as seed when perms > 1",
      celltype_universe = universe,
      threshold_mode = opts$threshold_mode,
      thresholds = thresholds,
      analysis_regions = lapply(names(region_defs), function(x) {
        list(analysis_region = x, component_rois = region_defs[[x]])
      }),
      blas_thread_environment = as.list(THREAD_ENV),
      created_at = format(Sys.time(), tz = "UTC", usetz = TRUE)
    ),
    file.path(paths$output_root, "run_metadata.json")
  )
  resolved <- opts
  resolved$input_assignment_parquet <- identity$path
  resolved$total_cpus <- cpus$total_cpus
  resolved$n_jobs <- cpus$n_jobs
  atomic_write_json(resolved, file.path(paths$output_root, "resolved_arguments.json"))
  yaml::write_yaml(
    resolved, file.path(paths$output_root, "resolved_config.yaml")
  )
}

`%||%` <- function(x, y) if (is.null(x)) y else x

validate_options <- function(opts) {
  if (is.null(opts$input_assignment_parquet)) {
    stop("--input-assignment-parquet is required")
  }
  if (!file.exists(opts$input_assignment_parquet)) stop("Input file does not exist")
  if (opts$permutations < 1L) stop("--permutations must be >= 1")
  if (opts$min_reference_cells < 1L || opts$min_total_cells < 1L) {
    stop("Cell-count thresholds must be >= 1")
  }
  if (!opts$threshold_mode %in% c("global_fixed_universe", "local")) {
    stop("--threshold-mode must be global_fixed_universe or local")
  }
  if (!(opts$alpha > 0 && opts$alpha < 1)) stop("--alpha must be between 0 and 1")
  if (!identical(opts$grid_shape, "square")) {
    stop("CRAWDAD Module 02 currently supports --grid-shape square only")
  }
  if (!all(opts$plot_formats %in% c("png", "pdf"))) {
    stop("--plot-formats may contain only png and pdf")
  }
  if (is.null(opts$output_dir) || is.null(opts$plot_dir)) {
    stop("--output-dir and --plot-dir are required")
  }
}

run_self_tests <- function() {
  opts <- parse_cli(c(
    "--regions", "lateral", "dorsomedial",
    "--combine-rois", "medial=dorsomedial,ventromedial",
    "--whole-tissue", "--scale-range", "100", "400",
    "--scale-interval", "100"
  ))
  defs <- parse_region_definitions(
    opts, c("lateral", "dorsomedial", "ventromedial", "outside_global_roi")
  )
  stopifnot(
    identical(defs$medial, c("dorsomedial", "ventromedial")),
    "outside_global_roi" %in% defs$whole_tissue,
    identical(opts$scales_um, c(100, 200, 300, 400))
  )
  test_dat <- data.frame(
    sample = rep(c("late", "early"), each = 4),
    z_height = rep(c(200, 100), each = 4),
    primary_roi = rep(c("lateral", "dorsomedial"), 4),
    cell_id = paste0("c", seq_len(8)),
    cell_type = rep(c("A", "B"), 4),
    x_aligned = seq_len(8),
    y_aligned = seq_len(8)
  )
  test_opts <- opts
  test_opts$samples <- "all"
  test_opts$permutations <- 2L
  test_opts$neighborhood_distances_um <- 50
  test_opts$seed <- 1L
  test_opts$drop_celltypes_below_count <- NA_integer_
  task_tab <- construct_tasks(
    test_dat, test_opts, defs, list(sha256 = "test"), "1.0.1",
    CELLTYPE_UNIVERSE_20
  )
  stopifnot(
    identical(unique(task_tab$sample), c("early", "late")),
    nrow(task_subset(test_dat, task_tab[1, ])) ==
      length(unique(task_subset(test_dat, task_tab[1, ])$cell_id))
  )
  s <- data.frame(
    donor = "d", sample = "s", depth = 1, analysis_region = "r",
    component_rois = "lateral",
    neighborhood_distance_um = 50, reference = "A", neighbor = "B",
    scale_um = rep(c(100, 200), each = 2),
    permutation = rep(1:2, 2), z_score = c(1, 3, 4, 6),
    threshold_used = 3, n_reference_cells = 30L,
    n_neighbor_cells = 40L, reference_eligible = TRUE
  )
  sm <- summarize_trends(s)
  fs <- first_significant_table(sm, 3)
  stopifnot(
    identical(sm$mean_z, c(2, 5)),
    isTRUE(all.equal(sm$sd_z, c(sqrt(2), sqrt(2)))),
    fs$first_significant_scale_um == 200
  )
  s_ineligible <- s
  s_ineligible$reference_eligible <- FALSE
  sm_ineligible <- summarize_trends(s_ineligible)
  stopifnot(
    all(is.na(sm_ineligible$primary_significant)),
    nrow(first_significant_table(sm_ineligible, 3)) == 0L
  )
  endpoint_opts <- parse_cli(c(
    "--scale-range", "100", "450", "--scale-interval", "100"
  ))
  stopifnot(identical(endpoint_opts$scales_um_resolved, c(100, 200, 300, 400)))
  stopifnot(inherits(
    try(parse_cli(c(
      "--scale-range", "100", "400", "--scale-interval", "100",
      "--scales-um", "100", "200", "400"
    )), silent = TRUE),
    "try-error"
  ))
  changed_opts <- test_opts
  changed_opts$alpha <- 0.01
  changed_task <- construct_tasks(
    test_dat, changed_opts, defs, list(sha256 = "test"), "1.0.1",
    CELLTYPE_UNIVERSE_20
  )
  stopifnot(task_tab$task_identity_hash[[1]] !=
              changed_task$task_identity_hash[[1]])
  identity_variants <- list(
    within(test_opts, threshold_mode <- "local"),
    within(test_opts, min_reference_cells <- min_reference_cells + 1L),
    within(test_opts, {
      scale_interval_um <- 50
      scales_um_resolved <- c(100, 150, 200, 250, 300, 350, 400)
      scales_um <- scales_um_resolved
    })
  )
  for (variant in identity_variants) {
    variant_task <- construct_tasks(
      test_dat, variant, defs, list(sha256 = "test"), "1.0.1",
      CELLTYPE_UNIVERSE_20
    )
    stopifnot(
      task_tab$task_identity_hash[[1]] !=
        variant_task$task_identity_hash[[1]]
    )
  }
  universe_task <- construct_tasks(
    test_dat, test_opts, defs, list(sha256 = "test"), "1.0.1",
    rev(CELLTYPE_UNIVERSE_20)
  )
  stopifnot(
    task_tab$task_identity_hash[[1]] != universe_task$task_identity_hash[[1]]
  )
  old <- Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA_character_)
  Sys.setenv(SLURM_CPUS_PER_TASK = "4")
  cpu_opts <- default_options()
  cpu_opts$total_cpus <- NA_integer_
  cpu_opts$n_jobs <- 4L
  stopifnot(resolve_cpus(cpu_opts)$n_jobs == 4L)
  cpu_opts$n_jobs <- 5L
  stopifnot(inherits(try(resolve_cpus(cpu_opts), silent = TRUE), "try-error"))
  if (is.na(old)) Sys.unsetenv("SLURM_CPUS_PER_TASK") else
    Sys.setenv(SLURM_CPUS_PER_TASK = old)
  stopifnot(all(Sys.getenv(names(THREAD_ENV)) == "1"))
  cat("All dependency-free Module 02 self-tests passed.\n")
}

main <- function() {
  opts <- parse_cli(commandArgs(trailingOnly = TRUE))
  if (opts$self_test) {
    run_self_tests()
    return(invisible(NULL))
  }
  validate_options(opts)
  cpus <- resolve_cpus(opts)
  crawdad_version <- as.character(utils::packageVersion("crawdad"))
  identity <- input_identity(opts$input_assignment_parquet)
  raw <- read_assignment(opts$input_assignment_parquet)
  dat <- validate_assignment(canonicalize_assignment(raw, opts))
  rm(raw)
  invisible(gc())
  if ("donor" %in% names(dat) && !is.null(opts$donor)) {
    dat <- dat[dat$donor == opts$donor, , drop = FALSE]
    if (!nrow(dat)) stop("No input rows remain for donor ", opts$donor)
  }
  universe <- read_celltype_universe(opts$celltype_universe_file)
  if (length(universe) != 20L || anyDuplicated(universe) ||
      anyNA(universe) || any(!nzchar(universe))) {
    stop(
      "Primary Module 02 analysis requires 20 ordered, unique cell-type labels"
    )
  }
  missing_from_universe <- setdiff(unique(dat$cell_type), universe)
  if (length(missing_from_universe)) {
    stop("Input contains cell types outside the universe: ",
         paste(missing_from_universe, collapse = ", "))
  }
  requested_types <- union(
    opts$include_celltypes,
    union(opts$exclude_celltypes, opts$highlight_neighbor_celltypes)
  )
  unknown_requested <- setdiff(requested_types, universe)
  if (length(unknown_requested)) {
    stop("Requested cell types outside the universe: ",
         paste(unknown_requested, collapse = ", "))
  }
  if (length(intersect(opts$include_celltypes, opts$exclude_celltypes))) {
    stop("A cell type cannot be both included and excluded")
  }
  region_defs <- parse_region_definitions(opts, sort(unique(dat$primary_roi)))
  thresholds <- data.frame(
    threshold_mode = c("global_fixed_universe", "local_definition"),
    alpha = opts$alpha,
    n_celltypes = c(length(universe), NA_integer_),
    n_directional_tests = c(length(universe)^2, NA_integer_),
    global_z_threshold = c(
      bonferroni_z(opts$alpha, length(universe)^2), NA_real_
    ),
    description = c(
      "One run-wide threshold based on all 20 x 20 directional pairs",
      "Task-specific crawdad::correctZBonferroni threshold; saved per result"
    ),
    stringsAsFactors = FALSE
  )
  tasks <- construct_tasks(
    dat, opts, region_defs, identity, crawdad_version, universe
  )
  paths <- list(
    output_root = file.path(
      normalizePath(opts$output_dir, mustWork = FALSE),
      paste0("neighdist_", opts$neighborhood_distances_um[[1]])
    ),
    plot_root = file.path(
      normalizePath(opts$plot_dir, mustWork = FALSE),
      paste0("neighdist_", opts$neighborhood_distances_um[[1]])
    )
  )
  paths$n_jobs <- min(cpus$n_jobs, nrow(tasks))
  software_versions <- c(
    R.version.string,
    paste("crawdad", crawdad_version),
    paste("sf", packageVersion("sf")),
    paste("BiocParallel", packageVersion("BiocParallel")),
    paste("dplyr", packageVersion("dplyr")),
    paste("ggplot2", packageVersion("ggplot2"))
  )
  for (region in names(region_defs)) {
    region_slug <- slugify(region)
    region_paths <- list(
      output_root = file.path(paths$output_root, region_slug),
      plot_root = file.path(paths$plot_root, region_slug),
      task_root = file.path(paths$output_root, region_slug, "tasks"),
      combined_root = file.path(paths$output_root, region_slug, "combined"),
      qc_root = file.path(paths$output_root, region_slug, "qc"),
      log_root = file.path(paths$output_root, region_slug, "logs"),
      n_jobs = paths$n_jobs
    )
    for (path in unlist(region_paths[names(region_paths) != "n_jobs"])) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
    }
    region_opts <- opts
    region_opts$output_dir <- region_paths$output_root
    region_opts$plot_dir <- region_paths$plot_root
    write_run_metadata(
      region_opts, region_paths, identity, universe, region_defs[region], cpus,
      crawdad_version, thresholds
    )
    atomic_write_csv(
      data.frame(order = seq_along(universe), cell_type = universe),
      file.path(region_paths$output_root, "celltype_universe.csv")
    )
    region_manifest <- tasks |>
      dplyr::filter(analysis_region == region)
    region_manifest$status <- "pending"
    atomic_write_csv(
      region_manifest, file.path(region_paths$output_root, "task_manifest.csv")
    )
    writeLines(
      software_versions,
      file.path(region_paths$output_root, "software_versions.txt")
    )
  }

  cat(
    "Running", nrow(tasks), "sample-region tasks with", paths$n_jobs,
    "parallel workers; CRAWDAD ncores=1 per worker.\n"
  )
  task_rows <- lapply(seq_len(nrow(tasks)), function(i) tasks[i, , drop = FALSE])
  runner <- function(row) run_one_task(
    row, dat, opts, universe, thresholds, paths
  )
  if (paths$n_jobs == 1L || opts$fail_fast) {
    results <- list()
    for (i in seq_along(task_rows)) {
      results[[i]] <- runner(task_rows[[i]])
      if (opts$fail_fast && results[[i]]$status == "failed") break
    }
  } else {
    results <- parallel::mclapply(
      task_rows, runner, mc.cores = paths$n_jobs, mc.preschedule = FALSE
    )
  }
  task_status <- dplyr::bind_rows(results)
  manifest <- tasks |>
    dplyr::left_join(
      task_status |>
        dplyr::select(task_id, status, n_cells, n_observed_celltypes, message),
      by = "task_id"
    )
  combined <- collect_task_outputs(tasks, opts, paths)
  for (region in names(region_defs)) {
    region_slug <- slugify(region)
    region_paths <- list(
      output_root = file.path(paths$output_root, region_slug),
      plot_root = file.path(paths$plot_root, region_slug),
      combined_root = file.path(paths$output_root, region_slug, "combined"),
      qc_root = file.path(paths$output_root, region_slug, "qc")
    )
    region_manifest <- manifest |>
      dplyr::filter(analysis_region == region)
    atomic_write_csv(
      region_manifest, file.path(region_paths$output_root, "task_manifest.csv")
    )
    region_combined <- lapply(combined, function(tab) {
      if ("analysis_region" %in% names(tab)) {
        tab[tab$analysis_region == region, , drop = FALSE]
      } else {
        tab
      }
    })
    region_status <- task_status |>
      dplyr::filter(analysis_region == region)
    write_combined(
      region_combined, region_status, opts, region_paths, thresholds
    )
  }

  failed <- task_status[task_status$status == "failed", ]
  cat(
    "Finished:", sum(task_status$status %in% c("completed", "resumed")),
    "completed/resumed,",
    sum(task_status$status == "skipped_ineligible"), "skipped,",
    nrow(failed), "failed.\n"
  )
  if (nrow(failed)) {
    stop("One or more tasks failed; see qc/task_status.csv and logs/")
  }
  invisible(task_status)
}

if (sys.nframe() == 0L) main()
