#!/bin/bash

#SBATCH --job-name=Annotate
#SBATCH --output=logs/Annotate.log
#SBATCH --error=logs/Annotate.log 
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
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
Rscript Banksy_Annotation.R

echo "********* Job Ends *********"
date
