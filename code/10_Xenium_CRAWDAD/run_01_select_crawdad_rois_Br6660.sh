#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"

export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_crawdad_rois_${USER}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_crawdad_rois_${USER}"

"${CONDA_ENV}/bin/python" \
    "${SCRIPT_DIR}/01_select_crawdad_rois.py" \
    --donor Br6660 \
    --xenium-coordinate-csv "${REPO_ROOT}/processed-data/05_Xenium_alignment/module02_Spateo_second_pass_Br6660/Br6660_second_pass_coordinates.csv" \
    --vhd-coordinate-dir "${REPO_ROOT}/processed-data/11_Xenium_VisiumHD_alignment/module_02_Spateo_second_pass" \
    --lateral-roi-name H1_M3TCP9V_D1 \
    --hull-method concave_rectangle \
    --hull-ratio 0.05 \
    --coordinate-unit micrometer \
    --geometry-tolerance-um2 1e-6 \
    --d1a-label D1_Island_A \
    --d1b-label D1_Island_B \
    --neighborhood-distances-um 50 100 \
    --tile-sizes-um 100 200 400 \
    --output-dir "${REPO_ROOT}/processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois" \
    --plot-dir "${REPO_ROOT}/plots/10_Xenium_CRAWDAD/module01_select_crawdad_rois" \
    --overwrite \
    "$@"
