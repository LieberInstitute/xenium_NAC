#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
PYTHON="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/bin/python"
N_JOBS="${SLURM_CPUS_PER_TASK:-11}"
if [[ "$#" -lt 1 ]]; then
  echo "Usage: bash ${BASH_SOURCE[0]} REGION [options]" >&2
  echo "REGION: lateral | dorsomedial | ventromedial | outside_global_roi" >&2
  exit 2
fi
REGION="$1"
shift
case "${REGION}" in
  lateral|dorsomedial|ventromedial|outside_global_roi) ;;
  *) echo "Unsupported region: ${REGION}" >&2; exit 2 ;;
esac

OUTPUT_DIR="${CRAWDAD_MODULE01_OUTPUT_DIR:-${REPO_ROOT}/processed-data/10_Xenium_CRAWDAD/module01_neighborhood_distance_diagnostics_Br6660}"
PLOT_DIR="${CRAWDAD_MODULE01_PLOT_DIR:-${REPO_ROOT}/plots/10_Xenium_CRAWDAD/module01_neighborhood_distance_diagnostics_Br6660}"
RUN_MODE=(--resume)
EXTRA_ARGS=("$@")
for ((i = 0; i < ${#EXTRA_ARGS[@]}; i++)); do
  if [[ "${EXTRA_ARGS[$i]}" == "--overwrite" || "${EXTRA_ARGS[$i]}" == "--resume" ]]; then
    RUN_MODE=()
  elif [[ "${EXTRA_ARGS[$i]}" == "--output-dir" ]]; then
    if ((i + 1 >= ${#EXTRA_ARGS[@]})); then
      echo "--output-dir requires a value" >&2
      exit 2
    fi
    OUTPUT_DIR="${EXTRA_ARGS[$((i + 1))]}"
  elif [[ "${EXTRA_ARGS[$i]}" == "--plot-dir" ]]; then
    if ((i + 1 >= ${#EXTRA_ARGS[@]})); then
      echo "--plot-dir requires a value" >&2
      exit 2
    fi
    PLOT_DIR="${EXTRA_ARGS[$((i + 1))]}"
  fi
done
TIMESTAMP="$(date '+%Y%m%d_%H%M%S')"
LOG_DIR="${OUTPUT_DIR}_logs"
LOG_FILE="${LOG_DIR}/run_01_${REGION}_${TIMESTAMP}.log"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-crawdad-distance-qc-${USER}"
mkdir -p "${MPLCONFIGDIR}" "${LOG_DIR}"

"${PYTHON}" "${SCRIPT_DIR}/01_diagnose_crawdad_neighborhood_distance.py" \
  --input-cell-table \
    "${REPO_ROOT}/processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois/cell_assignments/all_cells_roi_assignment.parquet" \
  --donor Br6660 \
  --samples all \
  --regions "${REGION}" \
  --regional-job \
  --reference-celltype D1_Island_A \
  --target-celltype D1_Island_B \
  --candidate-distances-um 25 50 75 100 150 200 300 \
  --island-bin-size-um 25 \
  --island-connectivity 8 \
  --island-closing-bins 1 \
  --min-island-cells 20 \
  --min-cells-per-type 20 \
  --n-jobs "${N_JOBS}" \
  --seed 1 \
  --output-dir "${OUTPUT_DIR}" \
  --plot-dir "${PLOT_DIR}" \
  "${RUN_MODE[@]}" \
  "${EXTRA_ARGS[@]}" 2>&1 | tee "${LOG_FILE}"

echo "Shared output: ${OUTPUT_DIR}"
echo "Shared plots: ${PLOT_DIR}"
echo "Region completed: ${REGION}"
echo "Global recommendation after all four jobs: ${OUTPUT_DIR}/neighborhood_distance_recommendation.md"
echo "Log: ${LOG_FILE}"

# Pilot example for one region and three representative AP slices.
# bash 10_Xenium_CRAWDAD/run_01_diagnose_crawdad_neighborhood_distance_Br6660.sh lateral \
#   --samples Br6660_NAc2_1090 Br6660_NAc6_3080 Br6660_NAc9_5080 \
#   --n-jobs 3 \
#   --output-dir /tmp/module01_neighborhood_distance_pilot \
#   --plot-dir /tmp/module01_neighborhood_distance_pilot_plots \
#   --overwrite
