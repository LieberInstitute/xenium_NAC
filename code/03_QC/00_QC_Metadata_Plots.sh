#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=00_plot
#SBATCH --output=logs/00_plot.log
#SBATCH --error=logs/00_plot.log 
#SBATCH --mem=25G
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4.x
Rscript 00_QC_Metadata_Plots.R

echo "********* Job Ends *********"
date
