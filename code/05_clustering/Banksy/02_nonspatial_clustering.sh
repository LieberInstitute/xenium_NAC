#!/bin/bash
#SBATCH --job-name=non_spatial_banksy_resgrid
#SBATCH --array=0-5
#SBATCH --cpus-per-task=1
#SBATCH --mem=150G
#SBATCH --time=48:00:00
#SBATCH --output=logs/%x_%A_%a_nonspatial.log

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5

# Define exactly six resolutions (indices 0..5)
RES_LIST=(0.25 0.5 0.75 1.0 1.25 1.5)

# Select the resolution corresponding to this array task
RES="${RES_LIST[$SLURM_ARRAY_TASK_ID]}"

echo ">>> $(date) : Running clusterBanksy with resolution=${RES}"
Rscript 02_nonspatial_clustering.R --res "${RES}"
