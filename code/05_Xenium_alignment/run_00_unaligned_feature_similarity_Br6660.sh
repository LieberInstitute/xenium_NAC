#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"

export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_unaligned_feature_${USER}"

"${CONDA_ENV}/bin/python" "${SCRIPT_DIR}/00_unaligned_feature_similarity.py" \
    --input-h5ad "${REPO_ROOT}/processed-data/02_build_spe/h5ad/spe_NormCounts_nucleus_normcounts.h5ad" \
    --celltype-csv "${REPO_ROOT}/processed-data/05_Clustering/Banksy_cell_types.csv" \
    --output-dir "${REPO_ROOT}/processed-data/05_Xenium_alignment/module00_unaligned_feature_similarity_Br6660" \
    --plot-dir "${REPO_ROOT}/plots/05_Xenium_alignment/module00_unaligned_feature_similarity_Br6660" \
    --donor Br6660 \
    --donor-column Donor \
    --sample-column Sample \
    --adata-cell-id-column cell_id \
    --celltype-csv-id-column cell_id \
    --celltype-csv-column CellType \
    --celltype-column CellType \
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
    --spatial-key spatial \
    --coordinate-unit um \
    --flip-y \
    --grid-sizes-um 200 300 400 \
    --min-cells-reference 5 \
    --min-cells-moving 5 \
    --min-overlap-grids 10 \
    --min-gene-nz-grids 5 \
    --min-gene-valid-offsets 1 \
    --min-gene-valid-scales 2 \
    --n-jobs 10 \
    --random-state 0 \
    --histogram-bins 30 \
    --plot-dpi 300 \
    --overwrite
