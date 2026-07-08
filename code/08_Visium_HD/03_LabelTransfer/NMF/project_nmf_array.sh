#!/bin/bash
#SBATCH --output=logs/nmf_proj_array.log
#SBATCH --error=logs/nmf_proj_array.log
# ============================================================================
# Submit the NMF-projection pipeline as three dependent SLURM stages:
#   1. presplit   (one job)            split the SPE by sample_id
#   2. project    (array, 1 per sample) project each sample onto the loadings
#   3. recombine  (one job)            merge results back into the SPE
#
# Run with:  bash run_projectr.sh
# ============================================================================

module load conda_R/4.5

set -euo pipefail

# ---- edit these ------------------------------------------------------------
N_SAMPLES=8                      # number of distinct sample_id values (= lines in the manifest)
RSCRIPT=project_nmf_array.R
# ----------------------------------------------------------------------------

# Stage 1 — split the SPE (holds the full object, so give it generous memory)
jid1=$(sbatch --parsable \
       --job-name=nmf_presplit --mem=64G --cpus-per-task=2 --time=24:00:00 \
       --output=logs/presplit_%j.out \
       --wrap="Rscript $RSCRIPT presplit")
echo "presplit submitted:  $jid1"

# Stage 2 — project each sample (array); each task holds ONE sample, so less memory
jid2=$(sbatch --parsable --dependency=afterok:$jid1 \
       --job-name=nmf_project --array=1-${N_SAMPLES} \
       --mem=32G --cpus-per-task=2 --time=48:00:00 \
       --output=logs/project_%A_%a.out \
       --wrap="Rscript $RSCRIPT project \$SLURM_ARRAY_TASK_ID")
echo "project array submitted: $jid2  (--array=1-${N_SAMPLES})"

# Stage 3 — recombine (loads the full SPE again to write the reducedDim)
jid3=$(sbatch --parsable --dependency=afterok:$jid2 \
       --job-name=nmf_recombine --mem=64G --cpus-per-task=2 --time=24:00:00 \
       --output=logs/recombine_%j.out \
       --wrap="Rscript $RSCRIPT recombine")
echo "recombine submitted: $jid3"

echo
echo "If a single array task fails (e.g. OOM), rerun just that index, then resubmit recombine:"
echo "  sbatch --array=<i> --mem=48G --wrap=\"Rscript $RSCRIPT project \\\$SLURM_ARRAY_TASK_ID\""
echo "  sbatch --wrap=\"Rscript $RSCRIPT recombine\""
