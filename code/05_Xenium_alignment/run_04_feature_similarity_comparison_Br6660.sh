#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(git -C "${SCRIPT_DIR}" rev-parse --show-toplevel)"
CONDA_ENV="/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo"

export PYTHONNOUSERSITE=1
export MPLCONFIGDIR="${TMPDIR:-/tmp}/matplotlib_compare_feature_${USER}"

"${CONDA_ENV}/bin/python" \
    "${SCRIPT_DIR}/04_feature_similarity_comparison.py" \
    --unaligned-pair-summary-csv "${REPO_ROOT}/processed-data/05_Xenium_alignment/module00_unaligned_feature_similarity_Br6660/Br6660_unaligned_feature_similarity_pair_summary.csv" \
    --first-pass-pair-summary-csv "${REPO_ROOT}/processed-data/05_Xenium_alignment/module03_Spateo_feature_similarity_first_pass_Br6660/Br6660_first_pass_feature_similarity_pair_summary.csv" \
    --second-pass-pair-summary-csv "${REPO_ROOT}/processed-data/05_Xenium_alignment/module03_Spateo_feature_similarity_second_pass_Br6660/Br6660_second_pass_feature_similarity_pair_summary.csv" \
    --output-dir "${REPO_ROOT}/processed-data/05_Xenium_alignment/module04_feature_similarity_comparison_Br6660" \
    --plot-dir "${REPO_ROOT}/plots/05_Xenium_alignment/module04_feature_similarity_comparison_Br6660" \
    --donor Br6660 \
    --summary-statistic multiscale_median \
    --plot-dpi 300 \
    --overwrite
