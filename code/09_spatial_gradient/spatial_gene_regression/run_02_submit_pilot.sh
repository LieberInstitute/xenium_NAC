#!/usr/bin/env bash
# Submit the frozen 14-gene pilot; each task fits both D1 annotations and both nested models.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${HERE}/../../.." && pwd)"
OUT="${ROOT}/processed-data/09_spatial_gradient/spatial_gene_regression/absolute_xy_2000um_density_v1"
MAN="${OUT}/checkpoints/pilot_gene_array_manifest.csv"

if [[ "${1:-}" == worker ]]; then
  if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "SLURM_ARRAY_TASK_ID missing" >&2
    exit 2
  fi

  GENE="$(awk -F, -v n="${SLURM_ARRAY_TASK_ID}" 'NR==n+1{gsub(/\r/,"",$2);print $2}' "${AP_MANIFEST:-$MAN}")"
  if [[ -z "${GENE}" ]]; then
    echo "No gene for task ${SLURM_ARRAY_TASK_ID}" >&2
    exit 2
  fi

  bash "${HERE}/run_pipeline.sh" --stage fit_gene --gene "${GENE}" --resume
  exit
fi

if [[ "${1:-}" == local_worker ]]; then
  GENE="${2:?Gene missing for local worker}"
  SAFE="${GENE//[^A-Za-z0-9._-]/_}"
  bash "${HERE}/run_pipeline.sh" \
    --stage fit_gene \
    --gene "${GENE}" \
    --resume \
    >"${OUT}/logs/local_pilot/${SAFE}.out" \
    2>"${OUT}/logs/local_pilot/${SAFE}.err"
  exit
fi

if [[ ! -f "${OUT}/checkpoints/PREPARE_M0_READY.txt" ]]; then
  echo "M0 preparation/CR2 hard gates have not passed; inspect ${OUT}/PREPARE_BLOCKER_REPORT.md" >&2
  exit 1
fi

N="$(( $(wc -l < "${MAN}") - 1 ))"
mkdir -p "${OUT}/logs/slurm_pilot"

if [[ "${PILOT_LOCAL_JOBS:-0}" =~ ^[1-9][0-9]*$ ]]; then
  if [[ ! -f "${OUT}/checkpoints/PREPARE_COMPLETE.txt" ]]; then
    echo "Full M0/M1 preparation has not passed" >&2
    exit 1
  fi

  mkdir -p "${OUT}/logs/local_pilot"
  export OMP_NUM_THREADS=1
  export OPENBLAS_NUM_THREADS=1
  export MKL_NUM_THREADS=1
  export VECLIB_MAXIMUM_THREADS=1
  export NUMEXPR_NUM_THREADS=1

  awk -F, 'NR>1{gsub(/\r/,"",$2);print $2}' "${MAN}" |
    xargs -r -n1 -P "${PILOT_LOCAL_JOBS}" \
      bash "${HERE}/run_02_submit_pilot.sh" local_worker

  echo "Completed ${N}-gene pilot with ${PILOT_LOCAL_JOBS} local workers"
  exit
fi

JOB="$(
  sbatch \
    --parsable \
    --export=ALL,AP_MANIFEST="${MAN}" \
    --chdir="${HERE}" \
    --partition="${PILOT_PARTITION:-shared}" \
    --array="1-${N}" \
    --cpus-per-task=1 \
    --mem="${PILOT_MEM:-4G}" \
    --time="${PILOT_TIME:-02:00:00}" \
    --job-name=Br6660_2kXY_pilot \
    --output="${OUT}/logs/slurm_pilot/%A_%a.out" \
    --error="${OUT}/logs/slurm_pilot/%A_%a.err" \
    "${HERE}/run_02_submit_pilot.sh" worker
)"

printf '%s\n' "${JOB}" | tee "${OUT}/logs/last_pilot_array_job_id.txt"
echo "Submitted frozen 14-gene array ${JOB}. After every task leaves squeue, run phase3."
