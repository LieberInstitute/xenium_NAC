#!/usr/bin/env bash
# Ordered orchestration record for the active pilot.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

case "${1:-}" in
  phase1)
    bash "${HERE}/run_01_prepare.sh"
    ;;
  phase2)
    bash "${HERE}/run_02_submit_pilot.sh"
    ;;
  phase3)
    bash "${HERE}/run_03_finalize_pilot.sh"
    ;;
  phase4 | production)
    bash "${HERE}/run_04_fit_production.sh"
    ;;
  phase5)
    bash "${HERE}/run_05_finalize_production.sh"
    ;;
  *)
    echo "Usage: bash run_full.sh {phase1|phase2|phase3|phase4|phase5}" >&2
    exit 2
    ;;
esac
