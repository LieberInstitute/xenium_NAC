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
    --output-dir "${REPO_ROOT}/processed-data/05_Xenium_alignment/module01_Spateo_initial_pass_Br6660" \
    --plot-dir "${REPO_ROOT}/plots/05_Xenium_alignment/module01_Spateo_initial_pass_Br6660" \
    --donor Br6660 \
    --donor-column Donor \
    --sample-column Sample \
    --sample-order \
        Br6660_NAc1_580 \
        Br6660_NAc2_1090 \
        Br6660_NAc3_1580 \
        Br6660_NAc4_2080 \
        Br6660_NAc5_2580 \
        Br6660_NAc6_3080 \
        Br6660_NAc7_3580 \
        Br6660_Nac10_4080 \
        Br6660_NAc8_4580 \
        Br6660_NAc9_5080 \
        Br6660_Nac11_5580 \
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
