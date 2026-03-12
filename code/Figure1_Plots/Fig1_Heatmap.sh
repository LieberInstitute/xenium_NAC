#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=Fig1_Heatmap
#SBATCH --output=logs/Fig1_Heatmap.log
#SBATCH --error=logs/Fig1_Heatmap.log 
#SBATCH --mem=30G
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
Rscript Fig1_Heatmap.R

echo "********* Job Ends *********"
date
