#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=02_snRNA_plot
#SBATCH --output=logs/02_plot.log
#SBATCH --error=logs/02_plot.log 
#SBATCH --mem=50G
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --time=7-00:00:00


echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
Rscript 02_plot_snRNA_results.R

echo "********* Job Ends *********"
date
