#!/usr/bin/env bash
# After explicit pilot approval: fit one production gene per single-core process.
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
OUT="${ROOT}/processed-data/09_spatial_gradient/spatial_gene_regression/absolute_xy_2000um_density_v1"
MAN="${OUT}/checkpoints/production_gene_manifest.csv"

if [[ "${1:-}" == worker ]]; then
  GENE="$(awk -F, -v n="${SLURM_ARRAY_TASK_ID:?}" 'NR==n+1{gsub(/\r|"/,"",$2);print $2}' "${AP_MANIFEST:-$MAN}")"
  bash "${HERE}/run_pipeline.sh" \
    --stage fit_production_gene \
    --gene "${GENE}" \
    --resume
  exit
fi

if [[ "${1:-}" == local_worker ]]; then
  GENE="${2:?}"
  SAFE="${GENE//[^A-Za-z0-9._-]/_}"
  bash "${HERE}/run_pipeline.sh" \
    --stage fit_production_gene \
    --gene "${GENE}" \
    --resume \
    >"${OUT}/logs/local_production/${SAFE}.out" \
    2>"${OUT}/logs/local_production/${SAFE}.err"
  exit
fi

bash "${HERE}/run_pipeline.sh" --stage prepare_production
N="$(( $(wc -l < "${MAN}") - 1 ))"

if [[ "${PRODUCTION_LOCAL_JOBS:-0}" =~ ^[1-9][0-9]*$ ]]; then
  mkdir -p "${OUT}/logs/local_production"
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  awk -F, 'NR>1{gsub(/\r|"/,"",$2);print $2}' "${MAN}" |
    xargs -r -n1 -P "${PRODUCTION_LOCAL_JOBS}" \
      bash "${HERE}/run_04_fit_production.sh" local_worker

  echo "Completed ${N}-gene production with ${PRODUCTION_LOCAL_JOBS} local workers"
  exit
fi

mkdir -p "${OUT}/logs/slurm_production"

JOB="$(
  sbatch \
    --parsable \
    --export=ALL,AP_MANIFEST="${MAN}" \
    --chdir="${HERE}" \
    --partition="${PRODUCTION_PARTITION:-shared}" \
    --array="1-${N}" \
    --cpus-per-task=1 \
    --mem="${PRODUCTION_MEM:-4G}" \
    --time="${PRODUCTION_TIME:-02:00:00}" \
    --job-name=Br6660_2kXY_prod \
    --output="${OUT}/logs/slurm_production/%A_%a.out" \
    --error="${OUT}/logs/slurm_production/%A_%a.err" \
    "${HERE}/run_04_fit_production.sh" worker
)"

printf '%s\n' "${JOB}" | tee "${OUT}/logs/last_production_array_job_id.txt"
echo "Submitted production array ${JOB}; run phase5 after it leaves squeue."
