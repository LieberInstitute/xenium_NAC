#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
PYTHON="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/bin/python"

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-4}"
PLOT_WORKERS="${CRAWDAD_PLOT_WORKERS:-${TOTAL_CPUS}}"
NEIGHBORHOOD_DISTANCE_UM="${CRAWDAD_NEIGHBORHOOD_DISTANCE_UM:-50}"
SCALE_MIN_UM="${CRAWDAD_SCALE_MIN_UM:-200}"
SCALE_MAX_UM="${CRAWDAD_SCALE_MAX_UM:-1000}"
SCALE_INTERVAL_UM="${CRAWDAD_SCALE_INTERVAL_UM:-100}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-crawdad-region-depth-${USER}"

"${PYTHON}" "${SCRIPT_DIR}/05_heatmap_region_depth.py" \
  --input-root \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --regions lateral dorsomedial ventromedial outside_global_roi \
  --neighborhood-distance-um "${NEIGHBORHOOD_DISTANCE_UM}" \
  --z-score-limit 10 \
  --scale-range-um "${SCALE_MIN_UM}" "${SCALE_MAX_UM}" \
  --scale-interval-um "${SCALE_INTERVAL_UM}" \
  --workers "${PLOT_WORKERS}" \
  --output-dir \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module05_heatmap_region_depth" \
  --plot-dir \
    "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module05_heatmap_region_depth" \
  --overwrite \
  "$@"

# Full 20 × 20 run:
# bash 10_Xenium_CRAWDAD/run_05_heatmap_region_depth_Br6660.sh
#
# Preview one directional pair:
# bash 10_Xenium_CRAWDAD/run_05_heatmap_region_depth_Br6660.sh \
#   --references DRD1_MSN --neighbors DRD2_MSN
