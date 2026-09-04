#!/bin/bash
# Run the entire scDRS pipeline as a SLURM dependency chain.
#
# Usage (from this directory, on a login node):
#   bash run_all.sh
#
# Each step is submitted with --dependency=afterok on the previous step, so
# the whole workflow runs unattended; a step failing cancels everything
# downstream. Per-job output goes to logs/<step>.log (see #SBATCH headers);
# this script's own submission record goes to logs/run_all_<timestamp>.log.
#
# Prerequisite (one-time): bash 00_setup_scdrs_env.sh  (builds conda env)

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p logs

RUN_LOG="logs/run_all_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "${RUN_LOG}") 2>&1

echo "********* scDRS pipeline submission *********"
date
echo "User: ${USER}; Host: ${HOSTNAME}"
echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"

# Recreate the conda env to ensure it's up to date
module load conda/3-24.3.0
if conda env list | grep -q "scdrs_nac"; then
    echo "Removing existing conda env 'scdrs_nac' ..."
    conda env remove -y -n scdrs_nac
fi
echo "Creating conda env 'scdrs_nac' via 00_setup_scdrs_env.sh ..."
bash 00_setup_scdrs_env.sh

j_manifest=$(sbatch --parsable 01_trait_manifest.sh)
echo "01_trait_manifest       -> job ${j_manifest}"

j_prep=$(sbatch --parsable --dependency=afterok:${j_manifest} 02_prep_magma_input.sh)
echo "02_prep_magma_input     -> job ${j_prep} (after ${j_manifest})"

j_annot=$(sbatch --parsable --dependency=afterok:${j_prep} 03_magma_annotate.sh)
echo "03_magma_annotate       -> job ${j_annot} (after ${j_prep})"

j_genes=$(sbatch --parsable --dependency=afterok:${j_annot} 04_magma_genes.sh)
echo "04_magma_genes (array)  -> job ${j_genes} (after ${j_annot})"

j_gs=$(sbatch --parsable --dependency=afterok:${j_genes} 05_make_gs_file.sh)
echo "05_make_gs_file         -> job ${j_gs} (after ${j_genes})"

j_score=$(sbatch --parsable --dependency=afterok:${j_gs} 06_scdrs_score.sh)
echo "06_scdrs_score (array)  -> job ${j_score} (after ${j_gs})"

j_down=$(sbatch --parsable --dependency=afterok:${j_score} 07_scdrs_downstream.sh)
echo "07_scdrs_downstream (array) -> job ${j_down} (after ${j_score})"

j_merge=$(sbatch --parsable --dependency=afterok:${j_down} 07_scdrs_downstream_merge.sh)
echo "07_scdrs_merge          -> job ${j_merge} (after ${j_down})"

j_viz=$(sbatch --parsable --dependency=afterok:${j_merge} 08_scdrs_viz.sh)
echo "08_scdrs_viz            -> job ${j_viz} (after ${j_merge})"

# 09 needs only the merged norm scores and the gene sets, not the plots, so it
# depends on the merge rather than on 08 and runs in parallel with it.
j_drivers=$(sbatch --parsable --dependency=afterok:${j_merge} 09_risk_gene_drivers.sh)
echo "09_risk_gene_drivers    -> job ${j_drivers} (after ${j_merge})"

# 10 plots the driver and expression tables 09 writes, so it must wait for 09
# specifically -- afterok, since on a 09 failure the expression file would be
# absent and 10 would silently emit a partial figure set.
j_driverplots=$(sbatch --parsable --dependency=afterok:${j_drivers} 10_driver_gene_plots.sh)
echo "10_driver_gene_plots    -> job ${j_driverplots} (after ${j_drivers})"

echo
echo "All jobs submitted. Monitor with: squeue -u ${USER}"
date
