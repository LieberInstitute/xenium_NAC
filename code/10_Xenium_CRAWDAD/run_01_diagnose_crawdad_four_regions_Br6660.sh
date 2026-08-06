#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
GENERIC_RUNNER="${SCRIPT_DIR}/run_01_diagnose_crawdad_neighborhood_distance_Br6660.sh"
REGIONS=(lateral dorsomedial ventromedial outside_global_roi)

if [[ -n "${SLURM_ARRAY_TASK_ID:-}" ]]; then
  INDEX="${SLURM_ARRAY_TASK_ID}"
  if ((INDEX < 0 || INDEX >= ${#REGIONS[@]})); then
    echo "SLURM_ARRAY_TASK_ID must be 0-3; got ${INDEX}" >&2
    exit 2
  fi
  REGION="${REGIONS[$INDEX]}"
elif [[ "$#" -ge 1 ]]; then
  REGION="$1"
  shift
else
  echo "Run one of these commands in each of the four 11-CPU jobs:" >&2
  for region in "${REGIONS[@]}"; do
    echo "bash ${BASH_SOURCE[0]} ${region}" >&2
  done
  echo "Alternatively submit this script as a SLURM array with indices 0-3." >&2
  exit 2
fi

case "${REGION}" in
  lateral|dorsomedial|ventromedial|outside_global_roi) ;;
  *) echo "Unsupported region: ${REGION}" >&2; exit 2 ;;
esac

exec bash "${GENERIC_RUNNER}" "${REGION}" "$@"
