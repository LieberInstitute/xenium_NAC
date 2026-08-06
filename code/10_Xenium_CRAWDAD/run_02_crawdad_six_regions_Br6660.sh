#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/run_02_crawdad_within_depth_Br6660.sh"

# 1. Lateral
bash "${RUNNER}" \
  --regions lateral

# 2. Dorsomedial
bash "${RUNNER}" \
  --regions dorsomedial

# 3. Ventromedial
bash "${RUNNER}" \
  --regions ventromedial

# 4. Medial: dorsomedial and ventromedial combined
bash "${RUNNER}" \
  --combine-rois medial=dorsomedial,ventromedial

# 5. Outside the global ROI
bash "${RUNNER}" \
  --regions outside_global_roi

# 6. Entire Xenium tissue, including outside_global_roi
bash "${RUNNER}" \
  --whole-tissue
