#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(crawdad)
  library(dplyr)
  library(ggplot2)
  library(jsonlite)
  library(readr)
})

OVERLAY_GROUPS <- list(
  c("Fibroblast_A", "Fibroblast_B", "Microglia_B"),
  c("DRD1_MSN", "DRD2_MSN", "MSN_Oligo", "D1_Island_A", "Astro_B"),
  c("Astro_A", "WM", "D1_Island_B", "Astrocyte_Oligo")
)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- normalizePath(sub("^--file=", "", script_arg[[1]]))
code_root <- dirname(dirname(dirname(script_path)))
project_root <- dirname(code_root)

default_task_dir <- file.path(
  project_root,
  "processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth",
  "neighdist_50/dorsomedial/tasks/Br6660_Nac10_4080"
)
default_output <- file.path(
  project_root,
  "plots/10_Xenium_CRAWDAD/figure_plots",
  "relationship_summary_Br6660_Nac10_4080_dorsomedial_neighdist_50",
  "relationship_summary.png"
)

usage <- function() {
  cat(
"Usage:
  Rscript 10_Xenium_CRAWDAD/figure_plots/03_relationship_summary_overlay_groups.R [options]

Options:
  --task-dir PATH              Module 02 task directory
  --output PATH                Output PNG path
  --scale-range-um MIN MAX     Included null scales [200 1000]
  --scale-interval-um NUM      Null-scale interval [100]
  --title TEXT                 Optional centered plot title
  --no-boxes                   Do not draw overlay-group boxes
  --overwrite                  Replace an existing output
  --help                       Show this message

The defaults select Br6660, depth 4080 um, dorsomedial, and neighDist 50 um.
", sep = "")
}

parse_args <- function(argv) {
  opts <- list(
    task_dir = default_task_dir,
    output = default_output,
    scale_range_um = c(200, 1000),
    scale_interval_um = 100,
    title = NULL,
    boxes = TRUE,
    overwrite = FALSE
  )
  if ("--help" %in% argv || "-h" %in% argv) {
    usage()
    quit(save = "no", status = 0)
  }
  i <- 1L
  while (i <= length(argv)) {
    token <- argv[[i]]
    if (token == "--overwrite") {
      opts$overwrite <- TRUE
      i <- i + 1L
    } else if (token == "--no-boxes") {
      opts$boxes <- FALSE
      i <- i + 1L
    } else if (token == "--task-dir" || token == "--output" ||
               token == "--scale-interval-um" || token == "--title") {
      if (i == length(argv)) stop(token, " requires a value")
      key <- gsub("-", "_", substring(token, 3), fixed = TRUE)
      opts[[key]] <- argv[[i + 1L]]
      i <- i + 2L
    } else if (token == "--scale-range-um") {
      if (i + 2L > length(argv)) stop(token, " requires MIN and MAX")
      opts$scale_range_um <- argv[c(i + 1L, i + 2L)]
      i <- i + 3L
    } else {
      stop("Unknown argument: ", token)
    }
  }
  opts$scale_range_um <- suppressWarnings(as.numeric(opts$scale_range_um))
  opts$scale_interval_um <- suppressWarnings(as.numeric(opts$scale_interval_um))
  if (length(opts$scale_range_um) != 2L ||
      any(!is.finite(opts$scale_range_um)) ||
      opts$scale_range_um[[1]] > opts$scale_range_um[[2]]) {
    stop("--scale-range-um requires two finite values: MIN MAX")
  }
  if (length(opts$scale_interval_um) != 1L ||
      !is.finite(opts$scale_interval_um) || opts$scale_interval_um <= 0) {
    stop("--scale-interval-um must be positive")
  }
  opts$task_dir <- normalizePath(opts$task_dir, mustWork = TRUE)
  opts$output <- path.expand(opts$output)
  opts
}

select_scales <- function(dat, bounds, interval) {
  requested <- seq(bounds[[1]], bounds[[2]], by = interval)
  if (!isTRUE(all.equal(tail(requested, 1), bounds[[2]]))) {
    stop("Scale range must be exactly divisible by --scale-interval-um")
  }
  available <- sort(unique(dat$scale[is.finite(dat$scale)]))
  missing <- requested[!vapply(
    requested,
    function(value) any(abs(available - value) < 1e-8),
    logical(1)
  )]
  if (length(missing)) {
    stop("Requested scales are absent: ", paste(missing, collapse = ", "))
  }
  dat |>
    filter(vapply(
      scale,
      function(value) any(abs(requested - value) < 1e-8),
      logical(1)
    ))
}

group_box <- function(group, x_order, y_order) {
  x_positions <- match(group, x_order)
  y_positions <- match(group, y_order)
  if (anyNA(x_positions) || anyNA(y_positions)) {
    stop("Overlay group is absent from the CRAWDAD relationship axes")
  }
  if (diff(range(x_positions)) + 1L != length(group) ||
      diff(range(y_positions)) + 1L != length(group)) {
    stop(
      "Overlay group is not contiguous after CRAWDAD reordering: ",
      paste(group, collapse = ", ")
    )
  }
  data.frame(
    xmin = min(x_positions) - 0.5,
    xmax = max(x_positions) + 0.5,
    ymin = min(y_positions) - 0.5,
    ymax = max(y_positions) + 0.5
  )
}

save_png <- function(plot, output) {
  dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
  temporary <- file.path(
    dirname(output),
    paste0(".", tools::file_path_sans_ext(basename(output)),
           ".tmp.", Sys.getpid(), ".png")
  )
  on.exit(if (file.exists(temporary)) unlink(temporary), add = TRUE)
  ggsave(temporary, plot, width = 11, height = 10, dpi = 300, bg = "white")
  if (!file.rename(temporary, output)) {
    if (!file.copy(temporary, output, overwrite = TRUE)) {
      stop("Could not save output: ", output)
    }
    unlink(temporary)
  }
}

main <- function() {
  opts <- parse_args(commandArgs(trailingOnly = TRUE))
  if (file.exists(opts$output) && !opts$overwrite) {
    stop("Output exists; use --overwrite: ", opts$output)
  }

  input_path <- file.path(opts$task_dir, "trends_permutation.csv.gz")
  identity_path <- file.path(opts$task_dir, "task_identity.json")
  if (!file.exists(input_path)) stop("Missing input: ", input_path)
  if (!file.exists(identity_path)) stop("Missing input: ", identity_path)

  dat <- read_csv(input_path, show_col_types = FALSE)
  unique_fields <- c(
    "sample", "depth", "analysis_region", "neighborhood_distance_um"
  )
  for (field in unique_fields) {
    observed <- unique(dat[[field]])
    if (length(observed) != 1L) {
      stop("Expected one ", field, " value in ", input_path)
    }
  }
  if (!isTRUE(all.equal(unique(dat$neighborhood_distance_um), 50))) {
    stop("Expected neighborhood_distance_um = 50 in ", input_path)
  }
  dat <- select_scales(dat, opts$scale_range_um, opts$scale_interval_um)

  threshold <- unique(dat$threshold_used[is.finite(dat$threshold_used)])
  if (length(threshold) != 1L) stop("Expected one finite Z threshold")
  identity <- read_json(identity_path, simplifyVector = TRUE)
  z_limit <- as.numeric(identity$resolved_identity$z_score_limit)
  if (length(z_limit) != 1L || !is.finite(z_limit)) {
    stop("Invalid z_score_limit in task_identity.json")
  }

  primary_perm <- dat |>
    filter(reference_eligible) |>
    transmute(
      perm = permutation, neighbor, Z, scale, reference, id = task_id
    )
  plot <- crawdad::vizRelationships(
    primary_perm,
    zSigThresh = threshold,
    zScoreLimit = z_limit,
    reorder = TRUE,
    symmetrical = FALSE,
    onlySignificant = FALSE,
    dotSizes = c(2, 10)
  )
  n_boxes <- 0L
  if (opts$boxes) {
    built <- ggplot_build(plot)
    x_order <- built$layout$panel_params[[1]]$x$get_labels()
    y_order <- built$layout$panel_params[[1]]$y$get_labels()
    boxes <- bind_rows(lapply(OVERLAY_GROUPS, group_box,
                              x_order = x_order, y_order = y_order))
    plot <- plot + geom_rect(
      data = boxes,
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
      inherit.aes = FALSE,
      fill = NA,
      color = "black",
      linewidth = 1.15
    )
    n_boxes <- nrow(boxes)
  }

  plot <- plot +
    labs(title = opts$title, x = "Reference", y = "Neighbor") +
    theme(
      plot.title = if (is.null(opts$title)) {
        element_blank()
      } else {
        element_text(size = 14, face = "bold", hjust = 0.5)
      },
      axis.text.x.top = element_text(
        angle = 90, hjust = 0, vjust = 0.5, size = 10,
        face = "bold",
        margin = margin(b = 3)
      ),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      plot.margin = margin(t = 8, r = 15, b = 10, l = 10)
    )

  save_png(plot, opts$output)
  cat("Saved:", opts$output, "\n")
  cat("Scales:", paste(sort(unique(primary_perm$scale)), collapse = ", "), "um\n")
  cat("Black boxes:", n_boxes, "\n")
}

main()

# Example usage from the code directory:
# /jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript \
#   10_Xenium_CRAWDAD/figure_plots/03_relationship_summary_overlay_groups.R
# Add --overwrite when replacing the existing manuscript plot.
