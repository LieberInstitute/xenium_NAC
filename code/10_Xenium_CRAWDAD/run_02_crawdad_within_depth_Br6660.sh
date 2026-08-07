#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
R_SCRIPT="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-10}"
N_JOBS="${CRAWDAD_N_JOBS:-${TOTAL_CPUS}}"
NEIGHBORHOOD_DISTANCE_UM="${CRAWDAD_NEIGHBORHOOD_DISTANCE_UM:-50}"
SCALE_MIN_UM="${CRAWDAD_SCALE_MIN_UM:-200}"
SCALE_MAX_UM="${CRAWDAD_SCALE_MAX_UM:-1000}"
SCALE_INTERVAL_UM="${CRAWDAD_SCALE_INTERVAL_UM:-100}"
PERMUTATIONS="${CRAWDAD_PERMUTATIONS:-3}"

export R_LIBS_USER="${R_LIBS_USER:-/users/jyao/R/4.5}"

"${R_SCRIPT}" "${SCRIPT_DIR}/02_crawdad_within_depth.R" \
  --input-assignment-parquet \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois/cell_assignments/all_cells_roi_assignment.parquet" \
  --donor Br6660 \
  --samples all \
  --scale-range "${SCALE_MIN_UM}" "${SCALE_MAX_UM}" \
  --scale-interval "${SCALE_INTERVAL_UM}" \
  --neighborhood-distances-um "${NEIGHBORHOOD_DISTANCE_UM}" \
  --permutations "${PERMUTATIONS}" \
  --seed 1 \
  --threshold-mode global_fixed_universe \
  --alpha 0.05 \
  --min-reference-cells 20 \
  --min-total-cells 100 \
  --highlight-neighbor-celltypes D1_Island_A D1_Island_B \
  --total-cpus "${TOTAL_CPUS}" \
  --n-jobs "${N_JOBS}" \
  --output-dir \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --plot-dir \
    "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --plot-formats png \
  --overwrite \
  "$@"

# Examples from the code directory:
# bash 10_Xenium_CRAWDAD/run_02_crawdad_within_depth_Br6660.sh --regions lateral
# bash 10_Xenium_CRAWDAD/run_02_crawdad_within_depth_Br6660.sh \
#   --combine-rois medial=dorsomedial,ventromedial
# bash 10_Xenium_CRAWDAD/run_02_crawdad_within_depth_Br6660.sh --whole-tissue
# bash 10_Xenium_CRAWDAD/run_02_crawdad_within_depth_Br6660.sh \
#   --regions lateral dorsomedial ventromedial outside_global_roi \
#   --combine-rois medial=dorsomedial,ventromedial --whole-tissue
