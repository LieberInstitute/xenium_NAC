#!/usr/bin/env bash
set -euo pipefail

# Fixed reference: Br6660_Nac10_4080 (Xenium, 4080 um)
# Moving slice:    VHD_H1_XKYDCP3_D1 (VisiumHD, 4120 um)
# Manual transform parameters derived from the first-pass overlay.
REALIGN_CELLTYPES=(D1_Island_A D1_Island_B)

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"

export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_spateo_${USER}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_spateo_${USER}"

"${CONDA_ENV}/bin/python" "${SCRIPT_DIR}/05_Spateo_Xenium_VHD_second_pass.py" \
    --donor Br6660 \
    --xenium-sample Br6660_Nac10_4080 \
    --vhd-sample VHD_H1_XKYDCP3_D1 \
    --vhd-z-height 4120 \
    --realign-celltypes "${REALIGN_CELLTYPES[@]}" \
    --first-pass-coordinates "${REPO_ROOT}/processed-data/11_Xenium_VisiumHD_alignment/module_01_Spateo_first_pass/Br6660_Nac10_4080-VHD_H1_XKYDCP3_D1/Br6660_Nac10_4080-VHD_H1_XKYDCP3_D1_first_pass_coordinates.csv" \
    --device cuda \
    --max-iter 200 \
    --partial-robust-level 100 \
    --overwrite
