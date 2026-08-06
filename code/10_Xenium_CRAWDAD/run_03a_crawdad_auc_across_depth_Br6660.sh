#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
R_SCRIPT="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript"
export R_LIBS_USER="${R_LIBS_USER:-/users/jyao/R/4.5}"
NEIGHBORHOOD_DISTANCE_UM="${CRAWDAD_NEIGHBORHOOD_DISTANCE_UM:-50}"

"${R_SCRIPT}" "${SCRIPT_DIR}/03a_crawdad_auc_across_depth.R" \
  --input-root "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --donor Br6660 \
  --regions lateral dorsomedial ventromedial outside_global_roi \
  --samples all \
  --neighborhood-distance-um "${NEIGHBORHOOD_DISTANCE_UM}" \
  --expected-scales-um 200 300 400 500 600 700 800 900 1000 \
  --expected-permutations 3 \
  --top-n-pairs 20 \
  --highlight-pairs D1_Island_A:D1_Island_B D1_Island_B:D1_Island_A \
  --output-dir "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module03_crawdad_auc_across_depth" \
  --plot-dir "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module03_crawdad_auc_across_depth" \
  --plot-formats png \
  "$@"

# Pilot example:
# bash 10_Xenium_CRAWDAD/run_03a_crawdad_auc_across_depth_Br6660.sh \
#   --regions lateral --samples Br6660_NAc1_580 Br6660_NAc2_1090 Br6660_NAc3_1580 --plot-formats png --overwrite
