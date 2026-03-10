#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=06_plot
#SBATCH --output=logs/06_plot.log
#SBATCH --error=logs/06_plot.log 
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
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
Rscript 06_plot_cell_level_spatial_clusters.R

echo "********* Job Ends *********"
date
