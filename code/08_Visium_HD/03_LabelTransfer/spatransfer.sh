#!/bin/bash
#SBATCH --job-name=spatransfer
#SBATCH --mem=700G
#SBATCH --output=logs/spatransfer.log
#SBATCH --error=logs/spatransfer.log
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
Rscript spatransfer.R

echo "********* Job Ends *********"
date
