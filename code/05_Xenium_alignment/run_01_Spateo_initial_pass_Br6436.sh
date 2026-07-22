#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"
export PYTHONNOUSERSITE=1
export NUMBA_CACHE_DIR="${TMPDIR:-/tmp}/numba_spateo_${USER}"
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_spateo_${USER}"

"${CONDA_ENV}/bin/python" "${SCRIPT_DIR}/01_Spateo_initial_pass.py" \
    --input-h5ad "${REPO_ROOT}/processed-data/02_build_spe/h5ad/spe_NormCounts_nucleus_normcounts.h5ad" \
    --celltype-csv "${REPO_ROOT}/processed-data/05_Clustering/Banksy_cell_types.csv" \
    --output-dir "${REPO_ROOT}/processed-data/05_Xenium_alignment/module01_Spateo_initial_pass_Br6436" \
    --plot-dir "${REPO_ROOT}/plots/05_Xenium_alignment/module01_Spateo_initial_pass_Br6436" \
    --donor Br6436 \
    --donor-column Donor \
    --sample-column Sample \
    --sample-order \
        Br6436_Nac1_650 \
        Br6436_Nac2_1150 \
        Br6436_Nac3_1650 \
        Br6436_Nac_4_2150 \
        Br6436_Nac_5_2650 \
        Br6436_Nac_6_3150 \
        Br6436_Nac_7_3660 \
        Br6436_Nac_8_4150 \
        Br6436_Nac_9_4650 \
        Br6436_Nac_10_5150 \
        Br6436_Nac_11_5650 \
    --adata-cell-id-column cell_id \
    --celltype-csv-id-column cell_id \
    --celltype-csv-column CellType \
    --celltype-column CellType \
    --spatial-key spatial \
    --key-added align_spatial \
    --device cuda \
    --max-iter 200 \
    --partial-robust-level 30 \
    --chunk-capacity 2 \
    --sparse-calculation-mode \
    --use-chunk \
    --flip-y \
    --overwrite
