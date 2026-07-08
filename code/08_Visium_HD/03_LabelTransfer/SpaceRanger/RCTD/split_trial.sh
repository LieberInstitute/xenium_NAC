#!/bin/bash
#SBATCH --job-name=plot
#SBATCH --mem=75G
#SBATCH --output=logs/split_trial.log
#SBATCH --error=logs/split_trial.log
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
Rscript split_trial.R

echo "********* Job Ends *********"
date
