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

input <- normalizePath(val("--input-dir"))

out <- val("--output-dir")

cfg_path <- normalizePath(val("--config"))

marker_path <- normalizePath(val("--marker-panel"))

axis_path <- normalizePath(val("--axis-metadata"))

dir.create(out, FALSE, TRUE)

dir.create(file.path(out, "checkpoints"), FALSE, TRUE)

dir.create(file.path(out, "logs"), FALSE, TRUE)

unlink(file.path(
  out,
  "checkpoints",
  c(
    "PREPARE_M0_READY.txt",
    "PREPARE_COMPLETE.txt",
    "PILOT_ACCEPTED.txt",
    "PRODUCTION_APPROVED.txt",
    "PRODUCTION_COMPLETE.txt"
  )
))

unlink(file.path(
  out,
  c(
    "06_common_absolute_tiles_2000um.csv",
    "06a_retained_common_tile_rows_2000um.csv"
  )
))

cfg <- yaml::read_yaml(cfg_path)

axyd_validate_config(cfg)

if (!file.copy(cfg_path, file.path(out, "00_config_frozen.yaml"), overwrite = TRUE)) {
  stop("Failed to freeze the active configuration")
}

active_files <- list.files(script_dir, recursive = TRUE, full.names = TRUE)

active_files <- active_files[file.info(active_files)$isdir %in% FALSE & !grepl(
  "/(Archive|__pycache__)/",
  active_files
) & grepl("[.](R|py|sh|yaml|csv|md)$", active_files)]

active_manifest <- data.frame(
  relative_path = substring(active_files, nchar(script_dir) + 2), size_bytes = file.info(active_files)$size,
  md5 = vapply(active_files, axyd_hash_file, character(1)), role = "active self-contained pipeline source",
  stringsAsFactors = FALSE
)

axyd_csv(active_manifest, file.path(out, "01_active_source_manifest.csv"))

if (!file.exists(file.path(input, "checkpoints", "GRID_INPUT_COMPLETE.txt"))) stop("Run 00_build_inputs.py first")

primary <- read.csv(gzfile(file.path(input, "primary_metadata.csv.gz")), stringsAsFactors = FALSE)

shifted <- read.csv(gzfile(file.path(input, "shifted_xy_plus1000um_metadata.csv.gz")), stringsAsFactors = FALSE)

genes <- read.delim(file.path(input, "genes.tsv"), stringsAsFactors = FALSE)$Gene

axis <- yaml::read_yaml(axis_path)

expected_ap <- as.numeric(unlist(cfg$samples))

expected_samples <- names(cfg$samples)

observed <- unique(primary[c("Sample", "AP_um")])

observed <- observed[match(expected_samples, observed$Sample), ]

axis_pass <- nrow(observed) == 11 && !anyNA(observed$AP_um) && identical(
  as.numeric(observed$AP_um),
  expected_ap
) && isTRUE(axis$mapping_verified) && grepl("anterior.*posterior", axis$axis_sign$AP,
  ignore.case = TRUE
)

audit <- data.frame(
  Sample = expected_samples, biological_order = seq_along(expected_samples), AP_um = expected_ap,
  AP_mm = expected_ap / 1000, observed_AP_um = observed$AP_um, AP_matches = observed$AP_um == expected_ap,
  ML_column = cfg$coordinate$ML_column, DV_column = cfg$coordinate$DV_column, coordinate_version = cfg$coordinate$coordinate_version,
  ML_direction = cfg$coordinate$orientation$ML, DV_direction = cfg$coordinate$orientation$DV, AP_direction = cfg$coordinate$orientation$AP,
  mapping_verified = axis$mapping_verified, audit_status = if (axis_pass) "PASS" else "FAIL"
)

axyd_csv(audit, file.path(out, "02_input_ap_coordinate_audit.csv"))

if (!axis_pass) stop("Biological AP/coordinate audit failed")

invisible(file.copy(file.path(input, "03_global_grid_definition_2000um.csv"), file.path(out, "03_global_grid_definition_2000um.csv"),
  overwrite = TRUE
))

invisible(file.copy(file.path(input, "04_cell_to_tile_reconciliation_2000um.csv"), file.path(out, "04_cell_to_tile_reconciliation_2000um.csv"),
  overwrite = TRUE
))

whole <- aggregate(cbind(n_cells, exposure) ~ Sample + CellType + AP_um, primary, sum)

whole$Sample <- as.character(whole$Sample)

whole$CellType <- as.character(whole$CellType)

sp <- axyd_build_support(primary, whole, cfg, "primary")

ss <- axyd_build_support(shifted, whole, cfg, "shifted_xy_plus1000um")

dec <- rbind(sp$decisions, ss$decisions)

candidates <- rbind(transform(sp$candidates, grid_variant = "primary"), transform(ss$candidates, grid_variant = "shifted_xy_plus1000um"))

axyd_csv(dec, file.path(out, "05_celltype_interval_and_support_2000um.csv"))

axyd_csv(candidates, file.path(out, "05a_candidate_celltype_intervals_2000um.csv"))

axyd_csv(
  rbind(sp$tile_recurrence, ss$tile_recurrence),
  file.path(out, "06_absolute_tile_recurrence_2000um.csv")
)

axyd_csv(rbind(sp$retention, ss$retention), file.path(out, "07_support_retention_by_slice_2000um.csv"))

axyd_csv(
  rbind(sp$rows, ss$rows),
  file.path(out, "06a_retained_repeated_tile_rows_2000um.csv")
)

area_path <- file.path(input, "cell_support_mask_area_2000um.csv.gz")

if (!file.exists(file.path(input, "checkpoints", "CELL_SUPPORT_MASK_COMPLETE.txt"))) stop("Cell-coordinate support mask checkpoint missing")

area_primary <- axyd_validate_area_and_density(sp$rows, area_path, cfg)

area_shifted <- axyd_validate_area_and_density(ss$rows, area_path, cfg)

area <- axyd_bind(list(area_primary$area, area_shifted$area))

density <- axyd_bind(list(area_primary$density, area_shifted$density))

decomp <- axyd_bind(list(area_primary$decomposition, area_shifted$decomposition))

axyd_csv(area, file.path(out, "08_valid_tissue_area_by_slice_tile_2000um.csv"))

axyd_csv(density, file.path(out, "09_celltype_density_by_slice_tile_2000um.csv"))

axyd_csv(decomp, file.path(out, "10_density_within_between_decomposition.csv"))

all_area <- read.csv(gzfile(area_path), stringsAsFactors = FALSE)

axyd_csv(all_area, file.path(out, "08a_cell_coordinate_support_area_all_parameter_sets.csv"))

invisible(file.copy(file.path(input, "cell_support_mask_parameters.csv"), file.path(out, "08b_cell_support_mask_parameters.csv"),
  overwrite = TRUE
))

invisible(file.copy(file.path(input, "cell_support_mask_parameter_sensitivity.csv"), file.path(out, "08c_cell_support_mask_area_sensitivity.csv"),
  overwrite = TRUE
))

sensitivity_sets <- setdiff(unique(all_area$mask_parameter_set), as.character(cfg$tissue_area$primary_parameter_set))

sensitivity_rows <- axyd_bind(lapply(sensitivity_sets, function(ps) {
  a <- axyd_validate_area_and_density(sp$rows, area_path, cfg, ps)
  b <- axyd_validate_area_and_density(ss$rows, area_path, cfg, ps)
  z <- rbind(a$rows, b$rows)
  z$mask_parameter_set <- ps
  z
}))

if (nrow(sensitivity_rows)) {
  sensitivity_rows$Sample <- as.character(sensitivity_rows$Sample)
  sensitivity_rows$absolute_tile <- factor(sensitivity_rows$absolute_tile)
  sensitivity_rows$AP_centered_mm <- ave(sensitivity_rows$AP_um, interaction(
    sensitivity_rows$CellType,
    sensitivity_rows$grid_variant
  ), FUN = function(x) (x - mean(unique(x))) / 1000)
}

panel <- read.csv(marker_path, stringsAsFactors = FALSE)

req <- c("marker_order", "source_label", "Gene", "source_file", "source_definition", "in_366_gene_panel_required")

if (!all(req %in% names(panel)) || anyDuplicated(panel[c("source_label", "Gene")]) || !setequal(
  unique(panel$source_label),
  c("D1_Island_A", "D1_Island_B")
)) {
  stop("Frozen marker panel is invalid")
}

if (any(table(panel$source_label) > as.integer(cfg$pilot$maximum_markers_per_label))) stop("Frozen marker cap exceeded")

source_paths <- file.path(script_dir, panel$source_file)

if (any(!file.exists(source_paths))) {
  stop("Self-contained marker source snapshot missing: ", paste(panel$source_file[!file.exists(source_paths)],
    collapse = ","
  ))
}

for (p in unique(source_paths)) {
  snapshot <- read.csv(p, stringsAsFactors = FALSE)
  if (!all(c("source_label", "Gene") %in% names(snapshot)) || !setequal(paste(panel$source_label[source_paths ==
    p], panel$Gene[source_paths == p]), paste(snapshot$source_label, snapshot$Gene))) {
    stop("Marker panel does not match self-contained source snapshot: ", p)
  }
}

panel$in_366_gene_panel <- panel$Gene %in% genes

if (any(axyd_bool(panel$in_366_gene_panel_required) & !panel$in_366_gene_panel)) {
  stop(
    "Frozen marker absent from 366-gene panel: ",
    paste(panel$Gene[!panel$in_366_gene_panel], collapse = ",")
  )
}

panel$marker_panel_md5 <- axyd_hash_file(marker_path)

panel$source_snapshot_md5 <- vapply(source_paths, axyd_hash_file, character(1))

panel$ordering_rule <- cfg$pilot$marker_order_rule

axyd_csv(panel, file.path(out, "11_frozen_D1_A_B_marker_panel.csv"))

union_genes <- unique(panel$Gene[order(panel$marker_order)])

manifest <- data.frame(task_index = seq_along(union_genes), Gene = union_genes, checkpoint_stem = sprintf(
  "%02d_%s",
  seq_along(union_genes), gsub("[^A-Za-z0-9._-]", "_", union_genes)
), stringsAsFactors = FALSE)

axyd_csv(manifest, file.path(out, "checkpoints", "pilot_gene_array_manifest.csv"))

all_rows <- rbind(area_primary$rows, area_shifted$rows)

if (any(
  all_rows$n_cells < as.integer(cfg$support$minimum_cells_per_slice_celltype_tile) |
    all_rows$exposure <= as.numeric(cfg$support$minimum_exposure_per_slice_celltype_tile)
)) {
  stop("Retained repeated-tile rows violate the per-slice tile qualification threshold")
}

retained_recurrence <- aggregate(
  Sample ~ CellType + grid_variant + absolute_tile,
  all_rows,
  function(x) length(unique(x))
)
names(retained_recurrence)[names(retained_recurrence) == "Sample"] <- "n_slices_present"

if (any(
  retained_recurrence$n_slices_present <
    as.integer(cfg$support$minimum_slices_per_absolute_tile)
)) {
  stop("Retained absolute tile violates the minimum four-slice recurrence contract")
}

all_rows$Sample <- as.character(all_rows$Sample)

all_rows$CellType <- as.character(all_rows$CellType)

all_rows$absolute_tile <- factor(all_rows$absolute_tile)

all_rows$AP_centered_mm <- ave(all_rows$AP_um, interaction(all_rows$CellType, all_rows$grid_variant),
  FUN = function(x) (x - mean(unique(x))) / 1000
)

rowaudit_base <- do.call(rbind, lapply(c("primary", "shifted_xy_plus1000um"), function(v) {
  do.call(rbind, lapply(c("D1_Island_A", "D1_Island_B"), function(ct) {
    z <- all_rows[all_rows$grid_variant == v & all_rows$CellType == ct, ]
    data.frame(
      CellType = ct, grid_variant = v, support_design = cfg$support$design,
      minimum_slices_per_absolute_tile = as.integer(cfg$support$minimum_slices_per_absolute_tile),
      n_slices = length(unique(z$Sample)), n_tiles = length(unique(z$absolute_tile)),
      n_rows = nrow(z), M0_rows_frozen = nrow(z), M1_rows_frozen = if (all(c(
        "density_within",
        "density_between"
      ) %in% names(z))) {
        nrow(z)
      } else {
        0
      }, M0_M1_identical_rows = area_primary$available && area_shifted$available, support_status = dec$status[match(paste(
        ct,
        v
      ), paste(dec$CellType, dec$grid_variant))], density_status = if ((v == "primary" &&
        area_primary$available) || (v != "primary" && area_shifted$available)) {
        "available"
      } else {
        "blocked"
      }
    )
  }))
}))

rowaudit <- merge(
  expand.grid(Gene = union_genes, CellType = c("D1_Island_A", "D1_Island_B"), grid_variant = c(
    "primary",
    "shifted_xy_plus1000um"
  ), stringsAsFactors = FALSE), rowaudit_base,
  by = c("CellType", "grid_variant"),
  all.x = TRUE, sort = FALSE
)

rowaudit$M0_M1_same_row_key_contract <- rowaudit$M0_M1_identical_rows

rowaudit$row_key <- "Sample|absolute_tile|CellType|Gene"

axyd_csv(rowaudit, file.path(out, "12_pilot_model_rows_audit.csv"))

sources <- c(file.path(script_dir, "common.R"), file.path(script_dir, "00_build_inputs.py"), file.path(
  script_dir,
  "01_build_cell_support_mask.py"
), normalizePath(sub("^--file=", "", grep("^--file=", argv, value = TRUE)[1])))

identity <- axyd_identity(cfg_path, input, marker_path, sources)

cr2_ok <- TRUE

cr2_reason <- ""

invisible(tryCatch(axyd_require_cr2(), error = function(e) {
  cr2_ok <<- FALSE
  cr2_reason <<- conditionMessage(e)
}))

blockers <- data.frame(gate = c(
  "cell_coordinate_support_area_primary", "cell_coordinate_support_area_shifted",
  "CR2_Satterthwaite_quasipoisson_validation"
), status = c(
  if (area_primary$available) "PASS" else "BLOCKED",
  if (area_shifted$available) "PASS" else "BLOCKED", if (cr2_ok) "PASS" else "BLOCKED"
), reason = c(
  area_primary$reason,
  area_shifted$reason, cr2_reason
), stringsAsFactors = FALSE)

axyd_csv(blockers, file.path(out, "checkpoints", "prepare_blockers.csv"))

prepared <- list(
  identity = identity, config = cfg, genes = genes, panel = panel, manifest = manifest,
  rows = all_rows, mask_sensitivity_rows = sensitivity_rows, support_decisions = dec, area_available = area_primary$available &&
    area_shifted$available, cr2_validated = cr2_ok, blockers = blockers, input_dir = input
)

axyd_rds(prepared, file.path(out, "checkpoints", "prepared_data.rds"))

area_ok <- area_primary$available && area_shifted$available

if (cr2_ok) writeLines("PASS", file.path(out, "checkpoints", "PREPARE_M0_READY.txt"))

if (cr2_ok && area_ok) {
  unlink(file.path(out, "PREPARE_BLOCKER_REPORT.md"))
  writeLines("PASS", file.path(out, "checkpoints", "PREPARE_COMPLETE.txt"))
  cat("2000-um absolute-XY density preparation passed all M0 and M1 hard gates using the cell-coordinate-derived analyzed-tissue support mask\n")
} else {
  writeLines(
    c(
      "# Preparation blocker report", "", paste0(
        "- ", blockers$gate, ": ", blockers$status,
        " — ", blockers$reason
      ), "", if (cr2_ok) "M0 is ready; M1 is not testable until the cell-coordinate support-area construction passes." else "No AP inference may run until CR2/Satterthwaite validation passes.",
      "No DAPI, morphology image, CellType-specific footprint, gene expression, or regression result was used to define analyzed-tissue support."
    ),
    file.path(out, "PREPARE_BLOCKER_REPORT.md")
  )
  if (!cr2_ok) {
    stop("Preparation audits completed but CR2/Satterthwaite is BLOCKED; inspect PREPARE_BLOCKER_REPORT.md")
  }
  cat("Preparation complete for M0; M1 density inference is BLOCKED. Inspect PREPARE_BLOCKER_REPORT.md\n")
}
