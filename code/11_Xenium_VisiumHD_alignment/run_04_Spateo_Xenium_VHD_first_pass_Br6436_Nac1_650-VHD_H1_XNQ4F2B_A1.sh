#!/usr/bin/env bash
set -euo pipefail

# Run 3: fixed Br6436_Nac1_650 (Xenium, 650 um);
# moving VHD_H1_XNQ4F2B_A1 (VisiumHD, 630 um).
SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"
export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_spateo_${USER}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_spateo_${USER}"

"${CONDA_ENV}/bin/python" "${SCRIPT_DIR}/04_Spateo_Xenium_VHD_first_pass.py" \
    --donor Br6436 \
    --xenium-sample Br6436_Nac1_650 \
    --vhd-sample VHD_H1_XNQ4F2B_A1 \
    --vhd-z-height 630 \
    --xenium-h5ad "${REPO_ROOT}/processed-data/02_build_spe/h5ad/spe_NormCounts_nucleus_normcounts.h5ad" \
    --celltype-csv "${REPO_ROOT}/processed-data/05_Clustering/Banksy_cell_types.csv" \
    --aligned-coordinates "${REPO_ROOT}/processed-data/05_Xenium_alignment/module02_Spateo_second_pass_Br6436/Br6436_second_pass_coordinates.csv" \
    --vhd-h5ad "${REPO_ROOT}/processed-data/HD_Full_Analysis/h5ad/VHD_H1_XNQ4F2B_A1_spaceranger_square008um_spatial_um.h5ad" \
    --output-dir "${REPO_ROOT}/processed-data/11_Xenium_VisiumHD_alignment/module_01_Spateo_first_pass" \
    --plot-dir "${REPO_ROOT}/plots/11_Xenium_VisiumHD_alignment/module_01_Spateo_first_pass" \
    --device cuda \
    --max-iter 500 \
    --partial-robust-level 100 \
    --chunk-capacity 2 \
    --overwrite
