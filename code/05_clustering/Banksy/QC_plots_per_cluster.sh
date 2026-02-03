#!/bin/bash

#SBATCH --job-name=QC_plots
#SBATCH --output=logs/QC_plots.log
#SBATCH --error=logs/QC_plots.log 
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
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
Rscript QC_plots_per_cluster.R

echo "********* Job Ends *********"
date
