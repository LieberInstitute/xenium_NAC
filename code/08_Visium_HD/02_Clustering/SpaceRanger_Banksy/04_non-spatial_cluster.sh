#!/bin/bash
#SBATCH --job-name=non-spatial_cluster
#SBATCH --array=0-4
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=96:00:00
#SBATCH --output=logs/%x_%A_%a.log

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5

# Five resolutions (indices 0..4)
RES_LIST=(0.4 0.6 0.8 1.0 1.2)

RES="${RES_LIST[$SLURM_ARRAY_TASK_ID]}"

echo ">>> $(date) : Running clusterBanksy with resolution=${RES}"
Rscript 04_non-spatial_cluster.R --res "${RES}"
