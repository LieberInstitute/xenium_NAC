#!/usr/bin/env Rscript
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

dir.create(plots, FALSE, TRUE)

dir.create(file.path(plots, "marker_fits"), FALSE, TRUE)

unlink(file.path(
  plots,
  paste0(
    "01_",
    c("D1_Island_A", "D1_Island_B"),
    "_absolute_common_support_by_slice_2000um.png"
  )
))

pngplot <- function(path, w = 2200, h = 1600, expr) {
  png(path, width = w, height = h, res = 200)
  on.exit(dev.off())
  force(expr)
  dev.off()
  on.exit(NULL)
}

support_rows <- read.csv(
  file.path(out, "06a_retained_repeated_tile_rows_2000um.csv"),
  stringsAsFactors = FALSE
)

ret <- read.csv(file.path(out, "07_support_retention_by_slice_2000um.csv"), stringsAsFactors = FALSE)

area <- read.csv(file.path(out, "08_valid_tissue_area_by_slice_tile_2000um.csv"), stringsAsFactors = FALSE)

density <- read.csv(file.path(out, "09_celltype_density_by_slice_tile_2000um.csv"), stringsAsFactors = FALSE)

decomp <- read.csv(file.path(out, "10_density_within_between_decomposition.csv"), stringsAsFactors = FALSE)

m0 <- read.csv(file.path(out, "13_pilot_M0_spatially_adjusted_results.csv"), stringsAsFactors = FALSE)

m1 <- read.csv(file.path(out, "14_pilot_M1_density_adjusted_results.csv"), stringsAsFactors = FALSE)

cmp <- read.csv(file.path(out, "15_pilot_M0_vs_M1_comparison.csv"), stringsAsFactors = FALSE)

loo <- read.csv(file.path(out, "17_pilot_leave_one_slice_out.csv"), stringsAsFactors = FALSE)

shift <- read.csv(file.path(out, "18_pilot_shifted_grid_results.csv"), stringsAsFactors = FALSE)

obs <- read.csv(file.path(out, "pilot_observed_and_fitted_rows.csv"), stringsAsFactors = FALSE)

masksens <- read.csv(file.path(out, "22_pilot_mask_parameter_sensitivity.csv"), stringsAsFactors = FALSE)

cts <- c("D1_Island_A", "D1_Island_B")

s <- support_rows[support_rows$grid_variant == "primary" & support_rows$CellType %in% cts, ]

for (ct in cts) {
  z <- s[s$CellType == ct, ]
  samps <- unique(z$Sample[order(z$AP_um)])
  pngplot(
    file.path(plots, paste0("01_", ct, "_unbalanced_repeated_tile_support_by_slice_2000um.png")),
    2800,
    2200,
    {
      par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
      for (sm in samps) {
        q <- z[z$Sample == sm, ]
        cols <- hcl.colors(50, "Viridis")[1 + pmin(49, floor(49 * (q$n_cells - min(q$n_cells)) / max(
          diff(range(q$n_cells)),
          .Machine$double.eps
        )))]
        plot(q$ML_center_um, q$DV_center_um,
          pch = 15, cex = 2.2, col = cols, asp = 1, xlab = "Aligned ML (µm)",
          ylab = "Aligned DV (µm)", main = paste(sub("^Br6660_", "", sm), unique(q$AP_um), "µm")
        )
        grid()
      }
      plot.new()
      legend(
        "center",
        "Points are retained repeated tiles present in this slice; color = CellType cells",
        bty = "n"
      )
    }
  )
}

if (all(c("valid_tissue_area_mm2", "Sample", "bx", "by") %in% names(area)) && all(c(
  "density_cells_per_mm2",
  "CellType"
) %in% names(density))) {
  for (ct in cts) {
    z <- density[density$grid_variant == "primary" & density$CellType == ct, ]
    samps <- unique(z$Sample[order(z$AP_um)])
    for (kind in c("valid_tissue_area_mm2", "density_cells_per_mm2")) {
      pngplot(file.path(plots, paste0("02_", ct, "_", kind, "_maps.png")), 2800, 2200, {
        par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))
        lim <- range(z[[kind]])
        for (sm in samps) {
          q <- z[z$Sample == sm, ]
          v <- q[[kind]]
          cols <- hcl.colors(50, "Viridis")[1 + pmin(49, floor(49 * (v - lim[1]) / max(
            diff(lim),
            .Machine$double.eps
          )))]
          plot(q$ML_center_um, q$DV_center_um,
            pch = 15, cex = 2.2, col = cols, asp = 1, xlab = "Aligned ML (µm)",
            ylab = "Aligned DV (µm)", main = paste(sub("^Br6660_", "", sm), kind)
          )
          grid()
        }
        plot.new()
      })
    }
  }
  pngplot(file.path(plots, "03_AP_vs_between_slice_density.png"), 2200, 1100, {
    par(mfrow = c(1, 2), mar = c(5, 5, 4, 2))
    for (ct in cts) {
      z <- unique(decomp[decomp$grid_variant == "primary" & decomp$CellType == ct, c(
        "Sample",
        "AP_um", "density_between"
      )])
      plot(z$AP_um / 1000, z$density_between,
        pch = 16, xlab = "AP (mm; A to P)", ylab = "Between-slice centered mean log density",
        main = ct
      )
      text(z$AP_um / 1000, z$density_between, sub("^Br6660_", "", z$Sample), pos = 3, cex = 0.65)
      grid()
    }
  })
}

ok <- obs$grid_variant == "primary"

pairs <- unique(obs[ok, c("Gene", "CellType")])

for (i in seq_len(nrow(pairs))) {
  gn <- pairs$Gene[i]
  ct <- pairs$CellType[i]
  z <- obs[ok & obs$Gene == gn & obs$CellType == ct, ]
  p <- file.path(plots, "marker_fits", paste0(
    gsub("[^A-Za-z0-9._-]", "_", paste(ct, gn, sep = "__")),
    ".png"
  ))
  pngplot(p, 1800, 1300, {
    par(mfrow = c(2, 1), mar = c(5, 5, 4, 2))
    q <- z[z$model_id == "M0", ]
    cols <- as.integer(factor(q$absolute_tile))
    plot(q$AP_um / 1000, q$observed_rate,
      pch = 16, col = cols, xlab = "AP (mm; A to P)", ylab = "Raw count / exposure",
      main = paste(ct, gn, "observed rates by absolute tile")
    )
    grid()
    shown <- character()
    for (m in c("M0", "M1")) {
      u <- z[z$model_id == m, ]
      if (nrow(u)) {
        a <- aggregate(cbind(raw_count, fitted_count, exposure) ~ AP_um, u, sum)
        lines(a$AP_um / 1000, a$fitted_count / a$exposure, type = "b", pch = if (m == "M0") {
          1
        } else {
          16
        }, col = if (m == "M0") {
          "#0072B2"
        } else {
          "#D55E00"
        })
        shown <- c(shown, m)
      }
    }
    legend("topright", paste(shown, "fitted"), pch = ifelse(shown == "M0", 1, 16), col = ifelse(shown ==
      "M0", "#0072B2", "#D55E00"), bty = "n")
    q0 <- z[z$model_id == "M0", ]
    plot(q0$AP_um / 1000, q0$raw_count,
      pch = 16, col = cols, xlab = "AP (mm; A to P)", ylab = "Raw count",
      main = "Zero counts and tile-level observations retained"
    )
    abline(h = 0, lty = 2, col = "red")
    grid()
  })
}

tested <- cmp$M0_model_status == "tested" & cmp$M1_model_status == "tested"

z <- cmp[tested, ]

lab <- paste(z$CellType, z$Gene, sep = " / ")

yy <- rev(seq_len(nrow(z)))

pngplot(file.path(plots, "04_M0_vs_M1_AP_coefficient_forest.png"), 2600, 2200, {
  if (!nrow(z)) {
    plot.new()
    text(0.5, 0.5, "M1 unavailable: no defensible tissue-area denominator")
  } else {
    par(mar = c(5, 20, 4, 2))
    xr <- range(c(z$M0_CI_AP_lower, z$M0_CI_AP_upper, z$M1_CI_AP_lower, z$M1_CI_AP_upper))
    plot(z$M0_beta_AP_per_mm, yy + 0.12,
      pch = 1, col = "#0072B2", xlim = xr, yaxt = "n", xlab = "AP beta per mm (CR2/Satterthwaite 95% CI)",
      ylab = "", main = "Spatially adjusted M0 versus density-adjusted M1"
    )
    segments(z$M0_CI_AP_lower, yy + 0.12, z$M0_CI_AP_upper, yy + 0.12, col = "#0072B2")
    points(z$M1_beta_AP_per_mm, yy - 0.12, pch = 16, col = "#D55E00")
    segments(z$M1_CI_AP_lower, yy - 0.12, z$M1_CI_AP_upper, yy - 0.12, col = "#D55E00")
    axis(2, yy, lab, las = 2, cex.axis = 0.65)
    abline(v = 0, lty = 2, col = "red")
    legend("topright", c("M0", "M1"), pch = c(1, 16), col = c("#0072B2", "#D55E00"), bty = "n")
  }
})

lp <- loo[loo$grid_origin_id == "primary" & loo$status == "tested", ]

pngplot(file.path(plots, "05_LOSO_AP_coefficients.png"), 2600, 2200, {
  if (!nrow(lp)) {
    plot.new()
    text(0.5, 0.5, "No successful LOSO estimates")
  } else {
    labs <- unique(paste(lp$CellType, lp$Gene, lp$model_id, sep = " / "))
    y <- match(paste(lp$CellType, lp$Gene, lp$model_id, sep = " / "), labs)
    par(mar = c(5, 20, 4, 2))
    plot(lp$LOSO_beta_AP_per_mm, rev(y),
      pch = 16, cex = 0.5, yaxt = "n", xlab = "LOSO AP beta per mm",
      ylab = "", main = "Leave-one-physical-slice-out estimates"
    )
    axis(2, rev(seq_along(labs)), labs, las = 2, cex.axis = 0.55)
    abline(v = 0, lty = 2, col = "red")
    grid()
  }
})

sc <- shift[shift$model_status_primary == "tested" & shift$model_status_shifted == "tested", ]

pngplot(file.path(plots, "06_primary_vs_shifted_grid_coefficients.png"), 1800, 1600, {
  if (!nrow(sc)) {
    plot.new()
    text(0.5, 0.5, "No comparable shifted-grid estimates")
  } else {
    plot(sc$beta_AP_per_mm_primary, sc$beta_AP_per_mm_shifted,
      pch = ifelse(sc$model_id == "M1",
        16, 1
      ), col = ifelse(sc$model_id == "M1", "#D55E00", "#0072B2"), xlab = "Primary grid AP beta/mm",
      ylab = "+1000-µm shifted grid AP beta/mm", main = "Prespecified grid-boundary sensitivity"
    )
    abline(0, 1, lty = 2, col = "red")
    grid()
    legend("topleft", c("M0", "M1"), pch = c(1, 16), col = c("#0072B2", "#D55E00"), bty = "n")
  }
})

ms <- masksens[masksens$grid_origin_id == "primary" & masksens$model_status == "tested", ]

pngplot(file.path(plots, "07_M1_mask_parameter_sensitivity.png"), 1800, 1600, {
  if (!nrow(ms)) {
    plot.new()
    text(0.5, 0.5, "No successful mask-parameter sensitivity estimates")
  } else {
    plot(ms$primary_mask_beta_AP_per_mm, ms$beta_AP_per_mm,
      pch = ifelse(grepl("75um", ms$mask_parameter_set),
        1, 16
      ), col = ifelse(ms$CellType == "D1_Island_A", "#0072B2", "#D55E00"), xlab = "M1 AP beta/mm: primary 100-µm mask",
      ylab = "M1 AP beta/mm: sensitivity mask", main = "Cell-coordinate support-mask parameter sensitivity"
    )
    abline(0, 1, lty = 2, col = "red")
    grid()
    legend("topleft", c("75-µm buffer", "125-µm buffer", "D1 Island A", "D1 Island B"), pch = c(
      1,
      16, 16, 16
    ), col = c("black", "black", "#0072B2", "#D55E00"), bty = "n")
  }
})

manifest <- data.frame(
  file = sort(list.files(plots, recursive = TRUE)), label = "required pilot diagnostic",
  format = "PNG", stringsAsFactors = FALSE
)

write.csv(manifest, file.path(plots, "figure_manifest.csv"), row.names = FALSE)

cat("Pilot PNG figures complete\n")
