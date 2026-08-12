#!/usr/bin/env bash
# Active Br6660 2000-µm unbalanced repeated-absolute-tile + density entry point.
set -euo pipefail
CODE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${CODE}/../../.." && pwd)"
BASE="${ROOT}/processed-data/09_spatial_gradient/spatial_gene_regression"
OUT="${BASE}/absolute_xy_2000um_density_v1"
CACHE="${OUT}/input_cache"
PLOTS="${ROOT}/plots/09_spatial_gradient/spatial_gene_regression/absolute_xy_2000um_density_v1"
CFG="${CODE}/config_Br6660.yaml"
MARKERS="${CODE}/pilot_D1_marker_panel.csv"
PY="${PYTHON_BIN:-/dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/bin/python}"
RS="${RSCRIPT_BIN:-/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript}"

stage=""
gene=""
resume=0

while (( $# )); do
  case "$1" in
    --stage)
      stage="$2"
      shift 2
      ;;
    --gene)
      gene="$2"
      shift 2
      ;;
    --resume)
      resume=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 2
      ;;
  esac
done

if [[ -z "${stage}" ]]; then
  echo "Missing --stage" >&2
  exit 2
fi

mkdir -p "${OUT}/logs" "${CACHE}" "${PLOTS}"
trap 's=$?; echo "ERROR: 2000um absolute-XY density stage='"${stage}"' line=${LINENO} status=${s}" >&2; exit $s' ERR

case "${stage}" in
  prepare)
    "${RS}" -e '
      if (!requireNamespace("clubSandwich", quietly = TRUE)) {
        stop("CR2_INFERENCE_BLOCKED: install clubSandwich in this R environment before phase1")
      }
    '

    if [[ "${REUSE_GRID_CACHE:-0}" == 1 ]]; then
      required_cache_files=(
        checkpoints/GRID_INPUT_COMPLETE.txt
        INPUT_CACHE_MANIFEST.json
        primary_metadata.csv.gz
        primary_counts.csv.gz
        shifted_xy_plus1000um_metadata.csv.gz
        shifted_xy_plus1000um_counts.csv.gz
      )

      for file in "${required_cache_files[@]}"; do
        if [[ ! -f "${CACHE}/${file}" ]]; then
          echo "Cannot reuse incomplete grid cache: ${CACHE}/${file}" >&2
          exit 1
        fi
      done

      echo "Reusing explicitly requested completed 2000-um count/grid cache"
    else
      "${PY}" "${CODE}/00_build_inputs.py" \
        --metadata "${BASE}/input_validation/merged_cell_metadata.csv.gz" \
        --h5ad "${ROOT}/processed-data/02_build_spe/h5ad/spe_counts_nucleus_cell_area_sf.h5ad" \
        --gene-filtering "${BASE}/pseudobulk_500um/gene_filtering.csv" \
        --validated-grid-config "${BASE}/pseudobulk_500um/voxel_config.json" \
        --output-dir "${CACHE}" \
        2>&1 | tee "${OUT}/logs/00_build_inputs.log"
    fi

    mkdir -p "${PLOTS}/cell_support_mask_overlays" "${OUT}/logs/matplotlib"

    MPLCONFIGDIR="${OUT}/logs/matplotlib" "${PY}" "${CODE}/01_build_cell_support_mask.py" \
      --metadata "${BASE}/input_validation/merged_cell_metadata.csv.gz" \
      --validated-grid-config "${BASE}/pseudobulk_500um/voxel_config.json" \
      --output-dir "${CACHE}" \
      --plot-dir "${PLOTS}/cell_support_mask_overlays" \
      2>&1 | tee "${OUT}/logs/01_build_cell_support_mask.log"

    "${RS}" "${CODE}/01_prepare.R" \
      --input-dir "${CACHE}" \
      --output-dir "${OUT}" \
      --config "${CFG}" \
      --marker-panel "${MARKERS}" \
      --axis-metadata "${BASE}/input_validation/axis_mapping_metadata.yaml" \
      2>&1 | tee "${OUT}/logs/01_prepare.log"
    ;;

  fit_gene)
    if [[ -z "${gene}" ]]; then
      echo "Missing --gene" >&2
      exit 2
    fi

    extra=()
    if (( resume )); then
      extra+=(--resume)
    fi

    "${RS}" "${CODE}/02_fit_pilot_gene.R" \
      --output-dir "${OUT}" \
      --gene "${gene}" \
      "${extra[@]}"
    ;;

  summarize_pilot)
    mkdir -p "${OUT}/logs/cache"
    XDG_CACHE_HOME="${OUT}/logs/cache" \
      "${RS}" "${CODE}/03_summarize_pilot.R" \
      --output-dir "${OUT}" \
      --plot-dir "${PLOTS}" \
      2>&1 | tee "${OUT}/logs/03_summarize_pilot.log"
    ;;

  prepare_production)
    "${RS}" "${CODE}/05_prepare_production.R" --output-dir "${OUT}"
    ;;

  fit_production_gene)
    if [[ -z "${gene}" ]]; then
      echo "Missing --gene" >&2
      exit 2
    fi

    extra=()
    if (( resume )); then
      extra+=(--resume)
    fi

    "${RS}" "${CODE}/05_fit_production_gene.R" \
      --output-dir "${OUT}" \
      --gene "${gene}" \
      "${extra[@]}"
    ;;

  summarize_production)
    mkdir -p "${OUT}/logs/cache"
    XDG_CACHE_HOME="${OUT}/logs/cache" \
      "${RS}" "${CODE}/06_summarize_production.R" \
      --output-dir "${OUT}" \
      --plot-dir "${PLOTS}" \
      2>&1 | tee "${OUT}/logs/06_summarize_production.log"
    ;;

  *)
    echo "Unknown stage: ${stage}" >&2
    exit 2
    ;;
esac
