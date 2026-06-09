#!/bin/bash
#SBATCH --job-name=11_non-spatial_plotting
#SBATCH --cpus-per-task=4
#SBATCH --mem=12G
#SBATCH --time=96:00:00
#SBATCH --output=logs/11_non-spatial_plotting.log
#SBATCH --error=logs/11_non-spatial_plotting.log

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5

Rscript 11_non-spatial_plotting_saving.R
