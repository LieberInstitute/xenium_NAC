#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=01_qc
#SBATCH --output=logs/01_qc.log
#SBATCH --error=logs/01_qc.log 
#SBATCH --mem=50G
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
Rscript 01_qc.R

echo "********* Job Ends *********"
date
