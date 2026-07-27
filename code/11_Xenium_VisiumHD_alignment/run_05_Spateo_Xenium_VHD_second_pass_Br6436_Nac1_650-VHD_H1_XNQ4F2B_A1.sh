#!/usr/bin/env bash
set -euo pipefail

# Fixed reference: Br6436_Nac1_650 (Xenium, 650 um)
# Moving slice:    VHD_H1_XNQ4F2B_A1 (VisiumHD, 630 um)
# Edit this array after inspecting the first-pass overlay.
REALIGN_CELLTYPES=(D1_Island_A D1_Island_B)

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"

export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_spateo_${USER}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_spateo_${USER}"

"${CONDA_ENV}/bin/python" "${SCRIPT_DIR}/05_Spateo_Xenium_VHD_second_pass.py" \
    --donor Br6436 \
    --xenium-sample Br6436_Nac1_650 \
    --vhd-sample VHD_H1_XNQ4F2B_A1 \
    --vhd-z-height 630 \
    --realign-celltypes "${REALIGN_CELLTYPES[@]}" \
    --first-pass-coordinates "${REPO_ROOT}/processed-data/11_Xenium_VisiumHD_alignment/module_01_Spateo_first_pass/Br6436_Nac1_650-VHD_H1_XNQ4F2B_A1/Br6436_Nac1_650-VHD_H1_XNQ4F2B_A1_first_pass_coordinates.csv" \
    --device cuda \
    --max-iter 200 \
    --partial-robust-level 100 \
    --overwrite
