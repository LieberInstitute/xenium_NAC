#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
R_SCRIPT="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript"

if [[ "$#" -eq 0 ]]; then
  echo "Usage: bash ${BASH_SOURCE[0]} --regions REGION [REGION ...] [options]" >&2
  echo "Example: bash ${BASH_SOURCE[0]} --regions lateral" >&2
  exit 2
fi

TOTAL_CPUS="${SLURM_CPUS_PER_TASK:-1}"
PLOT_WORKERS="${CRAWDAD_PLOT_WORKERS:-${TOTAL_CPUS}}"
NEIGHBORHOOD_DISTANCE_UM="${CRAWDAD_NEIGHBORHOOD_DISTANCE_UM:-50}"
SCALE_MIN_UM="${CRAWDAD_SCALE_MIN_UM:-200}"
SCALE_MAX_UM="${CRAWDAD_SCALE_MAX_UM:-1000}"
SCALE_INTERVAL_UM="${CRAWDAD_SCALE_INTERVAL_UM:-100}"

export R_LIBS_USER="${R_LIBS_USER:-/users/jyao/R/4.5}"
export XDG_CACHE_HOME="${TMPDIR:-/tmp}/crawdad-plot-cache-${USER}"
mkdir -p "${XDG_CACHE_HOME}/fontconfig"

"${R_SCRIPT}" "${SCRIPT_DIR}/02_crawdad_within_depth_plots.R" \
  --input-root \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --plot-root \
    "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --neighborhood-distances-um "${NEIGHBORHOOD_DISTANCE_UM}" \
  --scale-range-um "${SCALE_MIN_UM}" "${SCALE_MAX_UM}" \
  --scale-interval-um "${SCALE_INTERVAL_UM}" \
  --workers "${PLOT_WORKERS}" \
  --overwrite \
  "$@"

# Replace all eligible reference-trend plots for one region:
# bash run_02_crawdad_within_depth_plots.sh \
#   --regions lateral
#
# Replace plots for several regions:
# bash run_02_crawdad_within_depth_plots.sh \
#   --regions lateral dorsomedial ventromedial outside_global_roi
#
# Validate selected targets without replacing them:
# bash run_02_crawdad_within_depth_plots.sh \
#   --regions lateral --samples Br6660_NAc1_580 \
#   --references DRD1_MSN --dry-run
