#!/usr/bin/env bash
# Combine the completed pilot only; this never launches production.
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bash "${HERE}/run_pipeline.sh" --stage summarize_pilot

