#!/bin/bash
#SBATCH --job-name=addRCTD
#SBATCH --mem=150G
#SBATCH --output=logs/02_add_RCTD.log
#SBATCH --error=logs/02_add_RCTD.log
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
Rscript 02_addRCTD.R

echo "********* Job Ends *********"
date
