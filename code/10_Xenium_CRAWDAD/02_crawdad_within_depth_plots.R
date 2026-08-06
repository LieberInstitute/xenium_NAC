#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(crawdad)
  library(dplyr)
  library(ggplot2)
  library(readr)
})

CELLTYPE_COLORS <- c(
  Excitatory = "#DB1C5F",
  Microglia_A = "#0D87E4",
  D1_Island_B = "#00FF0D",
  Astro_A = "#FD00FD",
  Fibroblast_A = "#FFA7E2",
  Fibroblast_B = "#2AFECA",
  DRD2_MSN = "#7AAA16",
  DRD1_MSN = "#9400FF",
  Astro_B = "#823526",
  MSN_Oligo = "#F5DEC0",
  Inh_PVALB = "#93E5FF",
  OPC = "#7A1699",
  Microglia_B = "#FF0DBC",
  Inh_SST = "#C4B3FB",
  Astrocyte_Oligo = "#F80D2A",
  D1_Island_A = "#224B82",
  CHAT = "#FBA475",
  Ependymal = "#73690D",
  Microglia_Oligo = "#B62A7A",
  WM = "orange"
)

usage_text <- function() {
  cat(
"Usage:
  Rscript 02_crawdad_within_depth_plots.R [options]

Required:
  --input-root PATH                 Module 02 processed-data root
  --plot-root PATH                  Existing Module 02 plot root
  --regions NAME [NAME ...]         Regions to redraw

Selection:
  --samples all|SAMPLE [...]        Samples to redraw [all]
  --references all|TYPE [...]       Eligible references to redraw [all]
  --neighborhood-distances-um all|NUM [...]
                                    Neighborhood distances to redraw [all]
  --scale-range-um MIN MAX          Inclusive shuffle-scale range [all]
  --scale-interval-um NUM           Required interval with --scale-range-um

Execution:
  --workers INT                     Parallel task workers [1]
  --overwrite                       Required to replace existing PNG files
  --dry-run                         Validate and report targets without writing
  --help                            Show this message

Examples:
  Rscript 02_crawdad_within_depth_plots.R \\
    --input-root ../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --plot-root ../plots/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --regions lateral \\
    --workers 10 --overwrite

  Rscript 02_crawdad_within_depth_plots.R \\
    --input-root ../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --plot-root ../plots/10_Xenium_CRAWDAD/module02_crawdad_within_depth \\
    --regions dorsomedial ventromedial \\
    --samples Br6660_NAc1_580 Br6660_NAc2_1090 \\
    --references DRD1_MSN DRD2_MSN \\
    --dry-run
", sep = "")
}

default_options <- function() {
  list(
    input_root = NULL,
    plot_root = NULL,
    regions = character(),
    samples = "all",
    references = "all",
    neighborhood_distances_um = "all",
    scale_range_um = NULL,
    scale_interval_um = NULL,
    workers = 1L,
    overwrite = FALSE,
    dry_run = FALSE
  )
}

parse_cli <- function(argv) {
  if ("--help" %in% argv || "-h" %in% argv) {
    usage_text()
    quit(save = "no", status = 0)
  }
  opts <- default_options()
  boolean_keys <- c("overwrite", "dry_run")
  multi_keys <- c(
    "regions", "samples", "references", "neighborhood_distances_um",
    "scale_range_um"
  )
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (!startsWith(token, "--")) {
      stop("Unexpected positional argument: ", token)
    }
    key <- gsub("-", "_", substring(token, 3), fixed = TRUE)
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
      opts[[key]] <- values
      i <- j
    } else {
      opts[[key]] <- argv[[j]]
      i <- j + 1L
    }
  }
  opts$workers <- suppressWarnings(as.integer(opts$workers))
  if (is.na(opts$workers) || opts$workers < 1L) {
    stop("--workers must be a positive integer")
  }
  if (is.null(opts$input_root)) stop("--input-root is required")
  if (is.null(opts$plot_root)) stop("--plot-root is required")
  if (!length(opts$regions)) stop("--regions requires at least one region")
  if (xor(is.null(opts$scale_range_um), is.null(opts$scale_interval_um))) {
    stop("--scale-range-um and --scale-interval-um must be used together")
  }
  if (!is.null(opts$scale_range_um)) {
    opts$scale_range_um <- suppressWarnings(as.numeric(opts$scale_range_um))
    opts$scale_interval_um <- suppressWarnings(as.numeric(opts$scale_interval_um))
    if (length(opts$scale_range_um) != 2L ||
        any(!is.finite(opts$scale_range_um)) ||
        opts$scale_range_um[[1]] > opts$scale_range_um[[2]]) {
      stop("--scale-range-um requires two finite values: MIN MAX")
    }
    if (length(opts$scale_interval_um) != 1L ||
        !is.finite(opts$scale_interval_um) || opts$scale_interval_um <= 0) {
      stop("--scale-interval-um must be one positive finite value")
    }
    opts$selected_scales_um <- seq(
      opts$scale_range_um[[1]], opts$scale_range_um[[2]],
      by = opts$scale_interval_um
    )
    if (!isTRUE(all.equal(
      tail(opts$selected_scales_um, 1), opts$scale_range_um[[2]]
    ))) {
      stop("Scale range endpoints must align to --scale-interval-um")
    }
  } else {
    opts$selected_scales_um <- NULL
  }
  if (!opts$overwrite && !opts$dry_run) {
    stop("Use --overwrite to replace existing plots, or --dry-run to inspect targets")
  }
  opts$input_root <- normalizePath(opts$input_root, mustWork = TRUE)
  opts$plot_root <- normalizePath(opts$plot_root, mustWork = TRUE)
  opts
}

slugify <- function(x) {
  ans <- gsub("[^A-Za-z0-9._-]+", "_", x)
  ans <- gsub("^_+|_+$", "", ans)
  ifelse(nchar(ans) == 0L, "unnamed", ans)
}

select_requested <- function(available, requested, label) {
  available <- unique(as.character(available))
  if (identical(requested, "all")) return(available)
  if ("all" %in% requested) {
    stop(label, ": use either all or explicit values")
  }
  unknown <- setdiff(requested, available)
  if (length(unknown)) {
    stop(label, ": unavailable values: ", paste(unknown, collapse = ", "))
  }
  available[available %in% requested]
}

discover_task_specs <- function(opts) {
  nd_roots <- list.dirs(
    opts$input_root, recursive = FALSE, full.names = FALSE
  )
  nd_roots <- nd_roots[grepl("^neighdist_[0-9.]+$", nd_roots)]
  if (!length(nd_roots)) {
    stop("No distance-first neighdist_<D> directories found under ", opts$input_root)
  }
  available_distances <- sub("^neighdist_", "", nd_roots)
  if (identical(opts$neighborhood_distances_um, "all")) {
    selected_distances <- available_distances
  } else {
    requested_numeric <- suppressWarnings(
      as.numeric(opts$neighborhood_distances_um)
    )
    if (anyNA(requested_numeric)) {
      stop("--neighborhood-distances-um must be all or numeric values")
    }
    selected_distances <- as.character(requested_numeric)
    missing_distances <- setdiff(selected_distances, available_distances)
    if (length(missing_distances)) {
      stop("Unavailable neighborhood distances: ", paste(missing_distances, collapse = ", "))
    }
  }

  specs <- list()
  for (distance in selected_distances) {
    distance_name <- paste0("neighdist_", distance)
    distance_root <- file.path(opts$input_root, distance_name)
    available_regions <- list.dirs(
      distance_root, recursive = FALSE, full.names = FALSE
    )
    missing_regions <- setdiff(opts$regions, available_regions)
    if (length(missing_regions)) {
      stop(
        distance_name, " is missing requested regions: ",
        paste(missing_regions, collapse = ", ")
      )
    }
    for (region in opts$regions) {
      task_root <- file.path(distance_root, region, "tasks")
      if (!dir.exists(task_root)) stop("Missing task directory: ", task_root)
      available_samples <- list.dirs(
        task_root, recursive = FALSE, full.names = FALSE
      )
      samples <- select_requested(
        available_samples, opts$samples, paste0("Region ", region, " samples")
      )
      for (sample in samples) {
        task_dir <- file.path(task_root, sample)
        input_file <- file.path(task_dir, "trends_permutation.csv.gz")
        if (!file.exists(input_file) || file.info(input_file)$size <= 0L) {
          stop("Missing or empty CRAWDAD trend input: ", input_file)
        }
        specs[[length(specs) + 1L]] <- list(
          region = region,
          sample = sample,
          distance = as.numeric(distance),
          input_file = input_file,
          identity_file = file.path(task_dir, "task_identity.json"),
          plot_dir = file.path(
            opts$plot_root, distance_name, region, sample
          )
        )
      }
    }
  }
  specs
}

atomic_replace_png <- function(plot, target, width = 10, height = 7) {
  tmp <- file.path(
    dirname(target),
    paste0(".", tools::file_path_sans_ext(basename(target)),
           ".tmp.", Sys.getpid(), ".png")
  )
  on.exit(if (file.exists(tmp)) unlink(tmp), add = TRUE)
  ggplot2::ggsave(
    tmp, plot, width = width, height = height, dpi = 300, bg = "white"
  )
  if (!file.exists(tmp) || file.info(tmp)$size <= 1000L) {
    stop("Temporary plot was not created correctly: ", tmp)
  }
  if (!file.rename(tmp, target)) {
    if (!file.copy(tmp, target, overwrite = TRUE)) {
      stop("Could not replace plot: ", target)
    }
    unlink(tmp)
  }
  invisible(TRUE)
}

make_heatmap_plot <- function(dat, threshold, z_limit, title, ncol = 4L) {
  plot_dat <- dat
  plot_dat$scale <- factor(
    as.character(plot_dat$scale),
    levels = as.character(sort(unique(plot_dat$scale))),
    ordered = TRUE
  )
  plot_dat$Z <- pmax(-z_limit, pmin(z_limit, plot_dat$Z))
  plot_dat$Z[
    is.finite(plot_dat$Z) & abs(plot_dat$Z) < threshold
  ] <- 0
  ggplot2::ggplot(
    plot_dat, ggplot2::aes(x = scale, y = neighbor, fill = Z)
  ) +
    ggplot2::geom_tile() +
    ggplot2::facet_wrap(~reference, ncol = ncol) +
    ggplot2::scale_fill_gradient2(
      low = "blue", mid = "white", high = "red", midpoint = 0,
      limits = c(-z_limit, z_limit)
    ) +
    ggplot2::labs(
      x = "Shuffle scale (µm)", y = "Neighbor", fill = "Z score",
      title = title
    ) +
    ggplot2::theme_classic() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        size = if (ncol == 1L) 10 else 12
      ),
      axis.text.x = ggplot2::element_text(
        angle = 90, hjust = 1, vjust = 0.5
      )
    )
}

redraw_task <- function(spec, opts) {
  required <- c(
    "task_id", "sample", "depth", "analysis_region",
    "neighborhood_distance_um", "reference", "neighbor", "scale",
    "permutation", "Z", "reference_eligible", "threshold_used"
  )
  dat <- suppressMessages(
    readr::read_csv(spec$input_file, show_col_types = FALSE, progress = FALSE)
  )
  missing <- setdiff(required, names(dat))
  if (length(missing)) {
    stop(spec$input_file, " is missing columns: ", paste(missing, collapse = ", "))
  }
  if (!identical(unique(dat$sample), spec$sample)) {
    stop("Sample mismatch in ", spec$input_file)
  }
  if (!identical(unique(dat$analysis_region), spec$region)) {
    stop("Region mismatch in ", spec$input_file)
  }
  distances <- unique(dat$neighborhood_distance_um)
  if (length(distances) != 1L ||
      !isTRUE(all.equal(as.numeric(distances), spec$distance))) {
    stop("Neighborhood-distance mismatch in ", spec$input_file)
  }
  if (!is.null(opts$selected_scales_um)) {
    available_scales <- sort(unique(dat$scale[is.finite(dat$scale)]))
    missing_scales <- opts$selected_scales_um[
      !vapply(
        opts$selected_scales_um,
        function(x) any(abs(available_scales - x) < 1e-8),
        logical(1)
      )
    ]
    if (length(missing_scales)) {
      stop(
        "Requested scales are absent from ", spec$input_file, ": ",
        paste(missing_scales, collapse = ", ")
      )
    }
    dat <- dat |>
      dplyr::filter(
        vapply(
          scale,
          function(x) any(abs(opts$selected_scales_um - x) < 1e-8),
          logical(1)
        )
      )
  }
  thresholds <- unique(dat$threshold_used[is.finite(dat$threshold_used)])
  if (length(thresholds) != 1L) {
    stop("Expected one finite Z threshold in ", spec$input_file)
  }
  if (!file.exists(spec$identity_file)) {
    stop("Missing task identity metadata: ", spec$identity_file)
  }
  identity <- jsonlite::read_json(spec$identity_file, simplifyVector = TRUE)
  z_score_limit <- as.numeric(identity$resolved_identity$z_score_limit)
  if (length(z_score_limit) != 1L || !is.finite(z_score_limit) ||
      z_score_limit <= 0) {
    stop("Invalid z_score_limit in ", spec$identity_file)
  }

  eligible <- dat |>
    dplyr::filter(reference_eligible) |>
    dplyr::pull(reference) |>
    unique()
  references <- if (identical(opts$references, "all")) {
    eligible
  } else {
    unknown <- setdiff(opts$references, names(CELLTYPE_COLORS))
    if (length(unknown)) {
      stop("References missing from the standard palette: ",
           paste(unknown, collapse = ", "))
    }
    intersect(eligible, opts$references)
  }

  if (!length(references)) {
    return(data.frame(
      region = spec$region,
      sample = spec$sample,
      neighborhood_distance_um = spec$distance,
      reference = NA_character_,
      target = NA_character_,
      status = "no_selected_eligible_references",
      message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  depth <- unique(dat$depth)
  task_id <- unique(dat$task_id)
  if (length(depth) != 1L || length(task_id) != 1L) {
    stop("Expected one depth and task_id in ", spec$input_file)
  }
  base_title <- paste(
    spec$sample, paste0("depth ", depth, " µm"),
    spec$region, paste0("neighDist ", spec$distance, " µm"),
    sep = " | "
  )

  relationship_target <- file.path(spec$plot_dir, "relationship_summary.png")
  if (!file.exists(relationship_target)) {
    stop(
      "Refusing to create a new plot; expected replacement target is missing: ",
      relationship_target
    )
  }
  primary_perm <- dat |>
    dplyr::filter(reference_eligible) |>
    dplyr::transmute(
      perm = permutation, neighbor, Z, scale, reference, id = task_id
    )
  primary_mean <- primary_perm |>
    dplyr::group_by(neighbor, scale, reference) |>
    dplyr::summarise(Z = mean(Z), .groups = "drop")
  relationship_status <- "validated_dry_run"
  relationship_message <- NA_character_
  if (!opts$dry_run) {
    tryCatch(
      {
        if (nrow(primary_mean) > 0L && any(
          is.finite(primary_mean$Z) & abs(primary_mean$Z) >= thresholds[[1]]
        )) {
          relationship_plot <- crawdad::vizRelationships(
            primary_perm,
            zSigThresh = thresholds[[1]],
            zScoreLimit = z_score_limit,
            reorder = TRUE,
            symmetrical = FALSE,
            onlySignificant = FALSE,
            dotSizes = c(2, 10)
          ) +
            ggplot2::labs(title = base_title) +
            ggplot2::theme(
              plot.title = ggplot2::element_text(size = 13),
              axis.text.x.top = ggplot2::element_text(
                angle = 90, hjust = 0, vjust = 0.5, size = 8,
                margin = ggplot2::margin(b = 3)
              ),
              plot.margin = ggplot2::margin(t = 12, r = 15, b = 10, l = 10)
            )
        } else {
          relationship_plot <- ggplot2::ggplot() +
            ggplot2::annotate(
              "text", x = 0, y = 0,
              label = paste0(
                "No relationships passed |Z| >= ", thresholds[[1]], "."
              ),
              size = 5
            ) +
            ggplot2::xlim(-1, 1) +
            ggplot2::ylim(-1, 1) +
            ggplot2::labs(title = base_title) +
            ggplot2::theme_void() +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 13))
        }
        atomic_replace_png(
          relationship_plot, relationship_target, width = 11, height = 10
        )
        relationship_status <- "replaced"
      },
      error = function(e) {
        relationship_status <<- "failed"
        relationship_message <<- conditionMessage(e)
      }
    )
  }
  relationship_row <- data.frame(
    region = spec$region,
    sample = spec$sample,
    neighborhood_distance_um = spec$distance,
    reference = NA_character_,
    target = relationship_target,
    status = relationship_status,
    message = relationship_message,
    stringsAsFactors = FALSE
  )

  heatmap_all_target <- file.path(
    spec$plot_dir, "heatmap_all_references__custom_wrapper.png"
  )
  if (!file.exists(heatmap_all_target)) {
    stop(
      "Refusing to create a new plot; expected replacement target is missing: ",
      heatmap_all_target
    )
  }
  heatmap_all_status <- "validated_dry_run"
  heatmap_all_message <- NA_character_
  if (!opts$dry_run) {
    tryCatch(
      {
        heatmap_all_plot <- make_heatmap_plot(
          primary_mean |>
            dplyr::filter(reference %in% references),
          thresholds[[1]], z_score_limit, base_title, ncol = 4L
        )
        atomic_replace_png(
          heatmap_all_plot, heatmap_all_target, width = 14, height = 12
        )
        heatmap_all_status <- "replaced"
      },
      error = function(e) {
        heatmap_all_status <<- "failed"
        heatmap_all_message <<- conditionMessage(e)
      }
    )
  }
  heatmap_all_row <- data.frame(
    region = spec$region,
    sample = spec$sample,
    neighborhood_distance_um = spec$distance,
    reference = NA_character_,
    target = heatmap_all_target,
    status = heatmap_all_status,
    message = heatmap_all_message,
    stringsAsFactors = FALSE
  )

  rows <- list()
  for (i in seq_along(references)) {
    reference_type <- references[[i]]
    ref_slug <- slugify(reference_type)
    target <- file.path(
      spec$plot_dir,
      paste0("trends_ref_", ref_slug, "__official_crawdad.png")
    )
    if (!file.exists(target)) {
      stop(
        "Refusing to create a new plot; expected replacement target is missing: ",
        target
      )
    }
    ref_perm <- dat |>
      dplyr::filter(
        reference_eligible,
        reference == reference_type
      ) |>
      dplyr::transmute(
        perm = permutation, neighbor, Z, scale, reference,
        id = task_id
      )
    neighbors <- unique(ref_perm$neighbor)
    missing_colors <- setdiff(neighbors, names(CELLTYPE_COLORS))
    if (length(missing_colors)) {
      stop(
        "Neighbor cell types are missing from the standard palette: ",
        paste(missing_colors, collapse = ", ")
      )
    }
    title <- paste(
      base_title, paste0("reference: ", reference_type), sep = " | "
    )
    status <- "validated_dry_run"
    message <- NA_character_
    if (!opts$dry_run) {
      tryCatch(
        {
          plot <- crawdad::vizTrends(
            ref_perm,
            id = "neighbor",
            lines = TRUE,
            points = TRUE,
            withPerms = TRUE,
            facet = FALSE,
            zSigThresh = thresholds[[1]],
            colors = CELLTYPE_COLORS[neighbors],
            title = title
          )
          x_scale <- plot$scales$get_scales("x")
          if (!is.null(x_scale)) {
            x_scale$expand <- ggplot2::expansion(mult = c(0.035, 0.055))
          }
          plot <- plot +
            ggplot2::theme(
              plot.title = ggplot2::element_text(size = 12),
              axis.text.x = ggplot2::element_text(
                angle = 90, hjust = 1, vjust = 0.5,
                margin = ggplot2::margin(t = 5)
              ),
              plot.margin = ggplot2::margin(t = 10, r = 18, b = 18, l = 12)
            )
          atomic_replace_png(plot, target)
          status <- "replaced"
        },
        error = function(e) {
          status <<- "failed"
          message <<- conditionMessage(e)
        }
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      region = spec$region,
      sample = spec$sample,
      neighborhood_distance_um = spec$distance,
      reference = reference_type,
      target = target,
      status = status,
      message = message,
      stringsAsFactors = FALSE
    )

    heatmap_target <- file.path(
      spec$plot_dir,
      paste0("heatmap_ref_", ref_slug, "__custom_wrapper.png")
    )
    if (!file.exists(heatmap_target)) {
      stop(
        "Refusing to create a new plot; expected replacement target is missing: ",
        heatmap_target
      )
    }
    heatmap_status <- "validated_dry_run"
    heatmap_message <- NA_character_
    if (!opts$dry_run) {
      tryCatch(
        {
          heatmap_plot <- make_heatmap_plot(
            primary_mean |>
              dplyr::filter(reference == reference_type),
            thresholds[[1]], z_score_limit,
            paste0(base_title, "\nreference: ", reference_type),
            ncol = 1L
          )
          atomic_replace_png(
            heatmap_plot, heatmap_target, width = 7, height = 6
          )
          heatmap_status <- "replaced"
        },
        error = function(e) {
          heatmap_status <<- "failed"
          heatmap_message <<- conditionMessage(e)
        }
      )
    }
    rows[[length(rows) + 1L]] <- data.frame(
      region = spec$region,
      sample = spec$sample,
      neighborhood_distance_um = spec$distance,
      reference = reference_type,
      target = heatmap_target,
      status = heatmap_status,
      message = heatmap_message,
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(relationship_row, heatmap_all_row, rows)
}

main <- function(argv = commandArgs(trailingOnly = TRUE)) {
  opts <- parse_cli(argv)
  specs <- discover_task_specs(opts)
  if (!length(specs)) stop("No Module 02 task outputs matched the selection")
  cat(
    "Selected ", length(specs), " sample-region-distance tasks; workers=",
    opts$workers, if (opts$dry_run) "; mode=dry-run\n" else "; mode=overwrite\n",
    sep = ""
  )

  worker <- function(spec) {
    tryCatch(
      redraw_task(spec, opts),
      error = function(e) {
        data.frame(
          region = spec$region,
          sample = spec$sample,
          neighborhood_distance_um = spec$distance,
          reference = NA_character_,
          target = NA_character_,
          status = "failed",
          message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      }
    )
  }
  results <- if (opts$workers > 1L && length(specs) > 1L) {
    parallel::mclapply(
      specs, worker, mc.cores = min(opts$workers, length(specs)),
      mc.preschedule = FALSE
    )
  } else {
    lapply(specs, worker)
  }
  result <- dplyr::bind_rows(results)
  counts <- table(result$status, useNA = "ifany")
  cat(
    paste(names(counts), as.integer(counts), sep = "="),
    sep = "; "
  )
  cat("\n")
  failed <- result |>
    dplyr::filter(status == "failed")
  if (nrow(failed)) {
    print(as.data.frame(failed), row.names = FALSE)
    stop(nrow(failed), " plot replacement operation(s) failed")
  }
  invisible(result)
}

tryCatch(
  main(),
  error = function(e) {
    message("ERROR: ", conditionMessage(e))
    quit(save = "no", status = 1)
  }
)
