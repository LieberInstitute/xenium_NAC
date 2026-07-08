#!/bin/bash
#SBATCH --job-name=plot_nonspatial_cluster
#SBATCH --array=0-8
#SBATCH --cpus-per-task=4
#SBATCH --mem=50G
#SBATCH --time=24:00:00
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

RES_LIST=(0.4 0.6 0.8 1.0 1.2 1.4 1.6 1.8 2.0)

RES="${RES_LIST[$SLURM_ARRAY_TASK_ID]}"

echo ">>> $(date) : Plotting non-spatial clusters with resolution=${RES}"
Rscript 04.5_plot_nonspatial.R --res "${RES}"
