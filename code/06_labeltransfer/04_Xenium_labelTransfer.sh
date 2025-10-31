#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=04_xenium
#SBATCH --output=logs/04_xenium.log
#SBATCH --error=logs/04_xenium.log 
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
Rscript 04_Xenium_labelTransfer.R

echo "********* Job Ends *********"
date
