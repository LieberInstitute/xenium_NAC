#!/bin/bash
#SBATCH --job-name=presplit
#SBATCH --mem=64G
#SBATCH --cpus-per-task=2
#SBATCH --output=logs/presplit.log
#SBATCH --error=logs/presplit.log
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
Rscript 01_split_objects.R

echo "********* Job Ends *********"
date


