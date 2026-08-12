#!/usr/bin/env bash
# Combine the complete frozen production universe and apply global corrections.
set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"

bash "${HERE}/run_pipeline.sh" --stage summarize_production
