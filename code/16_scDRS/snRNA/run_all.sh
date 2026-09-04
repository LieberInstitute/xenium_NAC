#!/bin/bash
# Run the snRNA-seq scDRS pipeline as a SLURM dependency chain.
#
# Usage (from this directory):
#   bash run_all.sh
#
# Prerequisite: the scdrs_nac conda env must exist (shared with parent pipeline).

set -eo pipefail
cd "$(dirname "$0")"
mkdir -p logs

RUN_LOG="logs/run_all_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "${RUN_LOG}") 2>&1

echo "********* snRNA-seq scDRS pipeline submission *********"
date
echo "User: ${USER}; Host: ${HOSTNAME}"
echo "Git commit: $(git rev-parse HEAD 2>/dev/null || echo 'not a git repo')"

# Check that scdrs_nac env exists (don't recreate — shared with parent)
module load conda/3-24.3.0
if ! conda env list | grep -q "scdrs_nac"; then
    echo "ERROR: scdrs_nac conda env not found. Run the parent 16_scDRS/run_all.sh first."
    exit 1
fi
echo "Using existing scdrs_nac environment."

j_h5ad=$(sbatch --parsable 01_prep_h5ad.sh)
echo "01_prep_h5ad            -> job ${j_h5ad}"

j_prep=$(sbatch --parsable --dependency=afterok:${j_h5ad} 02_prep_magma_input.sh)
echo "02_prep_magma_input     -> job ${j_prep} (after ${j_h5ad})"

j_annot=$(sbatch --parsable --dependency=afterok:${j_prep} 03_magma_annotate.sh)
echo "03_magma_annotate       -> job ${j_annot} (after ${j_prep})"

j_genes=$(sbatch --parsable --dependency=afterok:${j_annot} 04_magma_genes.sh)
echo "04_magma_genes (array)  -> job ${j_genes} (after ${j_annot})"

j_gs=$(sbatch --parsable --dependency=afterok:${j_genes} 05_make_gs_file.sh)
echo "05_make_gs_file         -> job ${j_gs} (after ${j_genes})"

j_score=$(sbatch --parsable --dependency=afterok:${j_gs} 06_scdrs_score.sh)
echo "06_scdrs_score (array)  -> job ${j_score} (after ${j_gs})"

j_down=$(sbatch --parsable --dependency=afterok:${j_score} 07_scdrs_downstream.sh)
echo "07_scdrs_downstream     -> job ${j_down} (after ${j_score})"

j_merge=$(sbatch --parsable --dependency=afterok:${j_down} 07_scdrs_downstream_merge.sh)
echo "07_scdrs_merge          -> job ${j_merge} (after ${j_down})"

j_viz=$(sbatch --parsable --dependency=afterok:${j_merge} 08_scdrs_viz.sh)
echo "08_scdrs_viz            -> job ${j_viz} (after ${j_merge})"

echo
echo "All jobs submitted. Monitor with: squeue -u ${USER}"
date
