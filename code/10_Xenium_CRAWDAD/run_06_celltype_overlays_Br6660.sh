#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd -- "${SCRIPT_DIR}/.." && pwd)"
PYTHON="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/bin/python"

export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib-crawdad-celltype-overlays-${USER}"

"${PYTHON}" "${SCRIPT_DIR}/06_celltype_overlays.py" \
  --input-root \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois" \
  --samples \
    Br6660_NAc1_580 Br6660_NAc2_1090 Br6660_NAc3_1580 \
    Br6660_NAc4_2080 Br6660_NAc5_2580 Br6660_NAc6_3080 \
    Br6660_NAc7_3580 Br6660_Nac10_4080 Br6660_NAc8_4580 \
    Br6660_NAc9_5080 Br6660_Nac11_5580 \
  --regions \
    lateral dorsomedial ventromedial outside_global_roi \
  --output-dir \
    "${PROJECT_ROOT}/../processed-data/10_Xenium_CRAWDAD/module06_celltype_overlays" \
  --plot-dir \
    "${PROJECT_ROOT}/../plots/10_Xenium_CRAWDAD/module06_celltype_overlays" \
  --overwrite \
  "$@"

# Default run:
# bash run_06_celltype_overlays_Br6660.sh
