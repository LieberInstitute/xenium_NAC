#!/bin/bash
#SBATCH --job-name=add_labels
#SBATCH --mem=24G
#SBATCH --output=logs/add_labels.log
#SBATCH --error=logs/add_labels.log
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
Rscript 01_add_snRNA_labels.R

echo "********* Job Ends *********"
date
