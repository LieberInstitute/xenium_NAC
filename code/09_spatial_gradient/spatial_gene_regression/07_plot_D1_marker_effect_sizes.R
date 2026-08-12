#!/usr/bin/env Rscript

argv <- commandArgs(FALSE)
script_dir <- dirname(normalizePath(
  sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1])
))

source(file.path(script_dir, "common.R"))

args <- commandArgs(TRUE)

value_after <- function(option) {
  index <- match(option, args)
  if (is.na(index) || index == length(args)) {
    stop("Missing ", option)
  }
  args[index + 1]
}

output_dir <- value_after("--output-dir")
plot_dir <- value_after("--plot-dir")

production_dir <- file.path(output_dir, "production")
production_plot_dir <- file.path(plot_dir, "production")
dir.create(production_plot_dir, recursive = TRUE, showWarnings = FALSE)

results_path <- file.path(production_dir, "01_all_M0_M1_results.csv")
marker_path <- file.path(output_dir, "11_frozen_D1_A_B_marker_panel.csv")

if (!file.exists(results_path)) {
  stop("Production result table is missing: ", results_path)
}
if (!file.exists(marker_path)) {
  stop("Frozen D1 marker panel is missing: ", marker_path)
}

results <- read.csv(results_path, stringsAsFactors = FALSE, check.names = FALSE)
markers <- read.csv(marker_path, stringsAsFactors = FALSE, check.names = FALSE)

required_result_columns <- c(
  "Gene",
  "CellType",
  "model_id",
  "model_status",
  "beta_AP_per_mm",
  "CI_AP_lower",
  "CI_AP_upper",
  "AP_rate_ratio_per_mm",
  "AP_percent_rate_change_per_mm",
  "p_AP_nominal",
  "p_BH_tested_model_family",
  "p_Bonferroni_planned_model_family_7320",
  "p_Bonferroni_planned_combined_family_14640"
)
required_marker_columns <- c("marker_order", "source_label", "Gene")

missing_result_columns <- setdiff(required_result_columns, names(results))
missing_marker_columns <- setdiff(required_marker_columns, names(markers))

if (length(missing_result_columns)) {
  stop(
    "Production result table lacks: ",
    paste(missing_result_columns, collapse = ", ")
  )
}
if (length(missing_marker_columns)) {
  stop(
    "Frozen marker panel lacks: ",
    paste(missing_marker_columns, collapse = ", ")
  )
}

markers <- markers[order(markers$marker_order), required_marker_columns]
marker_key <- paste(markers$source_label, markers$Gene, sep = "\r")

own_label_results <- results[
  results$model_id %in% c("M0", "M1") &
    paste(results$CellType, results$Gene, sep = "\r") %in% marker_key,
]

marker_results <- merge(
  markers,
  own_label_results,
  by.x = c("source_label", "Gene"),
  by.y = c("CellType", "Gene"),
  all.x = TRUE,
  sort = FALSE
)

marker_results <- marker_results[
  order(marker_results$marker_order, match(marker_results$model_id, c("M0", "M1"))),
]

expected_rows <- 2L * nrow(markers)
if (nrow(marker_results) != expected_rows) {
  stop(
    "Expected two own-label model rows per frozen marker; found ",
    nrow(marker_results),
    " of ",
    expected_rows
  )
}
if (anyDuplicated(marker_results[c("source_label", "Gene", "model_id")])) {
  stop("Duplicated own-label marker/model result")
}
if (any(marker_results$model_status != "tested")) {
  failed <- marker_results[marker_results$model_status != "tested", ]
  stop(
    "Cannot plot untested marker/model rows: ",
    paste(paste(failed$source_label, failed$Gene, failed$model_id), collapse = "; ")
  )
}

finite_columns <- c("beta_AP_per_mm", "CI_AP_lower", "CI_AP_upper")
if (any(!is.finite(as.matrix(marker_results[finite_columns])))) {
  stop("Marker effect-size estimates or confidence limits are non-finite")
}

marker_results$significant_primary_M0_Bonferroni_7320_p05 <- (
  marker_results$model_id == "M0" &
    marker_results$p_Bonferroni_planned_model_family_7320 < 0.05
)
marker_results$display_model <- ifelse(
  marker_results$model_id == "M0",
  "M0: primary",
  "M1: density sensitivity"
)
marker_results$effect_interpretation <- ifelse(
  marker_results$beta_AP_per_mm > 0,
  "increases toward posterior",
  ifelse(
    marker_results$beta_AP_per_mm < 0,
    "increases toward anterior",
    "no estimated AP direction"
  )
)
marker_results$multiplicity_display_rule <- paste(
  "asterisk denotes M0 primary 7,320-test Bonferroni adjusted p < 0.05; M1 is sensitivity only"
)

output_columns <- c(
  "marker_order",
  "source_label",
  "Gene",
  "model_id",
  "display_model",
  "model_status",
  "beta_AP_per_mm",
  "CI_AP_lower",
  "CI_AP_upper",
  "AP_rate_ratio_per_mm",
  "AP_percent_rate_change_per_mm",
  "effect_interpretation",
  "p_AP_nominal",
  "p_BH_tested_model_family",
  "p_Bonferroni_planned_model_family_7320",
  "p_Bonferroni_planned_combined_family_14640",
  "significant_primary_M0_Bonferroni_7320_p05",
  "multiplicity_display_rule"
)

axyd_csv(
  marker_results[output_columns],
  file.path(production_dir, "09_D1_marker_AP_effect_sizes.csv")
)

plot_file <- file.path(
  production_plot_dir,
  "04_D1_Island_A_B_marker_AP_effect_sizes.png"
)

png(plot_file, width = 2400, height = 1500, res = 220)

layout(matrix(c(1, 2), nrow = 1), widths = c(1.15, 1))
par(
  mar = c(5.8, 6.4, 4.8, 1.2),
  oma = c(2.2, 0, 0.5, 0),
  mgp = c(3.1, 0.9, 0),
  tcl = -0.25,
  las = 1,
  family = "sans"
)

model_colors <- c(M0 = "#0072B2", M1 = "#D55E00")
model_symbols <- c(M0 = 16, M1 = 17)
model_offsets <- c(M0 = 0.13, M1 = -0.13)

for (celltype_index in seq_along(c("D1_Island_A", "D1_Island_B"))) {
  celltype <- c("D1_Island_A", "D1_Island_B")[celltype_index]
  panel <- marker_results[marker_results$source_label == celltype, ]
  marker_order <- markers$Gene[markers$source_label == celltype]
  marker_order <- rev(marker_order)
  marker_labels <- as.expression(
    lapply(marker_order, function(gene) bquote(italic(.(gene))))
  )
  y_position <- match(panel$Gene, marker_order)
  y_upper <- length(marker_order) + if (celltype_index == 2L) 1.5 else 0.5

  x_limits <- range(c(panel$CI_AP_lower, panel$CI_AP_upper, 0), finite = TRUE)
  x_padding <- max(diff(x_limits) * 0.12, 0.05)
  x_limits <- x_limits + c(-x_padding, x_padding)

  plot(
    NA,
    xlim = x_limits,
    ylim = c(0.5, y_upper),
    xlab = expression(beta[AP] ~ "(log rate ratio per mm; A " %->% " P)"),
    ylab = "",
    yaxt = "n",
    main = "",
    bty = "l"
  )
  axis(2, at = seq_along(marker_order), labels = marker_labels, tick = FALSE)
  top_ticks <- pretty(x_limits, n = 5)
  axis(
    3,
    at = top_ticks,
    labels = format(round(exp(top_ticks), 2), trim = TRUE),
    cex.axis = 0.75
  )
  mtext("Rate ratio per mm", side = 3, line = 2.1, cex = 0.72)
  mtext(
    gsub("_", " ", celltype),
    side = 3,
    line = 3.8,
    cex = 1.15,
    font = 2
  )
  abline(v = 0, lty = 2, col = "#666666", lwd = 1.2)
  abline(h = seq_along(marker_order), col = "#EEEEEE", lwd = 0.8)

  for (model_id in c("M0", "M1")) {
    model_rows <- panel$model_id == model_id
    model_y <- y_position[model_rows] + model_offsets[[model_id]]

    segments(
      panel$CI_AP_lower[model_rows],
      model_y,
      panel$CI_AP_upper[model_rows],
      model_y,
      col = model_colors[[model_id]],
      lwd = 2.2
    )
    points(
      panel$beta_AP_per_mm[model_rows],
      model_y,
      pch = model_symbols[[model_id]],
      col = model_colors[[model_id]],
      cex = 1.05
    )

    significant <- model_rows & panel$significant_primary_M0_Bonferroni_7320_p05
    if (any(significant)) {
      text(
        panel$CI_AP_upper[significant] + 0.02 * diff(x_limits),
        y_position[significant] + model_offsets[[model_id]],
        labels = "*",
        col = model_colors[[model_id]],
        cex = 1.25,
        font = 2
      )
    }
  }

  if (celltype_index == 2L) {
    legend(
      "topright",
      legend = c("M0: primary", "M1: density sensitivity"),
      col = model_colors,
      pch = model_symbols,
      lwd = 2,
      bty = "o",
      bg = "white",
      box.col = "white",
      cex = 0.68
    )
  }
}

mtext(
  paste0(
    "Points: AP effect sizes; lines: 95% CR2/Satterthwaite CI. ",
    "* M0 primary: 7,320-test Bonferroni adjusted p < 0.05; M1 is sensitivity only."
  ),
  side = 1,
  outer = TRUE,
  line = 0.4,
  cex = 0.78
)

dev.off()

cat("D1 Island A/B own-label marker AP effect-size plot complete\n")
