#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=01_qc
#SBATCH --output=logs/02_after_qc_plot.log
#SBATCH --error=logs/02_after_qc_plot.err 
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=jyao37@jhmi.edu

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
Rscript 02_after_qc_plot.R

echo "********* Job Ends *********"
date