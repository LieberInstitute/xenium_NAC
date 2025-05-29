#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=Banksy
#SBATCH --output=logs/Banksy.log
#SBATCH --error=logs/Banksy.log 
#SBATCH --mem=250G
#SBATCH --cpus-per-task=1
#SBATCH --time=5-00:00:00
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
Rscript Banksy_clustering.R

echo "********* Job Ends *********"
date
