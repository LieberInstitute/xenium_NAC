#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

bash "${SCRIPT_DIR}/run_03a_crawdad_auc_across_depth_Br6660.sh" "$@"
bash "${SCRIPT_DIR}/run_03b_crawdad_auc_across_region_Br6660.sh" "$@"
