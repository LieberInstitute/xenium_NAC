#!/usr/bin/env bash
# Interactive, one CPU. This performs all gene-independent audits and hard-gate checks.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bash "${HERE}/run_pipeline.sh" --stage prepare

