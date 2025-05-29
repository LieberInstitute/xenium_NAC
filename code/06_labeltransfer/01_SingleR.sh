#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=01_SingleR
#SBATCH --output=logs/SingleR.%a.log
#SBATCH --error=logs/SingleR.%a.log 
#SBATCH --mem=25G
#SBATCH --array=1-22
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
Rscript 01_SingleR.R

echo "********* Job Ends *********"
date
