#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
R_SCRIPT="/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript"
export R_LIBS_USER="${R_LIBS_USER:-/users/jyao/R/4.5}"
NEIGHBORHOOD_DISTANCE_UM="${CRAWDAD_NEIGHBORHOOD_DISTANCE_UM:-50}"

"${R_SCRIPT}" "${SCRIPT_DIR}/03b_crawdad_auc_across_region.R" \
  --input-root "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth" \
  --donor Br6660 \
  --comparison-sets full_exclusive_partition \
  --samples all \
  --neighborhood-distance-um "${NEIGHBORHOOD_DISTANCE_UM}" \
  --expected-scales-um 200 300 400 500 600 700 800 900 1000 \
  --expected-permutations 3 \
  --top-n-pairs 20 \
  --highlight-pairs D1_Island_A:D1_Island_B D1_Island_B:D1_Island_A \
  --output-dir "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module03_crawdad_auc_across_region" \
  --plot-dir "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module03_crawdad_auc_across_region" \
  --plot-formats png \
  "$@"

# Pilot example:
# bash 10_Xenium_CRAWDAD/run_03b_crawdad_auc_across_region_Br6660.sh \
#   --samples Br6660_NAc1_580 --plot-formats png --overwrite
