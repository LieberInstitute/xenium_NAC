#!/bin/bash
#SBATCH --job-name=compare
#SBATCH --mem=24G
#SBATCH --output=logs/compare.log
#SBATCH --error=logs/compare.log
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
Rscript compare_RCTD_singleR.R

echo "********* Job Ends *********"
date
