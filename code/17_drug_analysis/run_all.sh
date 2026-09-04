#!/bin/bash
# Run the drug2cell pipeline as a SLURM dependency chain.
#
# Usage (from this directory, on a login node):
#   bash run_all.sh
#
# Submits:
#   01_drug2cell -> 02_drug2cell_enrichment -> 03_curated_drug_sets -> 04_drug_viz

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

RUN_LOG="logs/run_all_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "${RUN_LOG}") 2>&1

echo "********* Drug analysis pipeline submission *********"
date
echo "User: ${USER}; Host: ${HOSTNAME}"
echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"

# Recreate the conda env to ensure it's up to date.
# Skip if any scDRS job that uses this env is queued or running. The list must
# cover EVERY 16_scDRS step that runs `source activate scdrs_nac` -- steps 01,
# 02 and 08 use conda_R instead and are unaffected. Destroying the env from
# under a running job would fail it mid-flight.
SCDRS_ENV_JOBS=05_make_gs,06_scdrs_score,07_scdrs_down,07_scdrs_merge,09_risk_gene_drivers
module load conda/3-24.3.0
if squeue -u "${USER}" -n "${SCDRS_ENV_JOBS}" -h 2>/dev/null | grep -q .; then
    echo "WARNING: scDRS jobs are still running — skipping env recreation to avoid breaking them."
    echo "Using existing 'scdrs_nac' environment."
else
    if conda env list | grep -q "scdrs_nac"; then
        echo "Removing existing conda env 'scdrs_nac' ..."
        conda env remove -y -n scdrs_nac
    fi
    echo "Creating conda env 'scdrs_nac' via ../16_scDRS/00_setup_scdrs_env.sh ..."
    bash ../16_scDRS/00_setup_scdrs_env.sh
fi

j_d2c=$(sbatch --parsable 01_drug2cell.sh)
echo "01_drug2cell            -> job ${j_d2c}"

j_enrich=$(sbatch --parsable --dependency=afterok:${j_d2c} 02_drug2cell_enrichment.sh)
echo "02_drug2cell_enrichment -> job ${j_enrich} (after ${j_d2c})"

j_curated=$(sbatch --parsable --dependency=afterok:${j_enrich} 03_curated_drug_sets.sh)
echo "03_curated_drug_sets    -> job ${j_curated} (after ${j_enrich})"

j_viz=$(sbatch --parsable --dependency=afterok:${j_curated} 04_drug_viz.sh)
echo "04_drug_viz             -> job ${j_viz} (after ${j_curated})"

echo
echo "All SLURM jobs submitted. Monitor with: squeue -u ${USER}"
date
