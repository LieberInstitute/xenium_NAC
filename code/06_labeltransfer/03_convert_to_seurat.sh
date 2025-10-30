#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=03_seurat
#SBATCH --output=logs/03_seurat.log
#SBATCH --error=logs/03_seurat.log 
#SBATCH --mem=100G
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

module load conda_R/4.5
Rscript 03_convert_to_seurat.R

echo "********* Job Ends *********"
date
