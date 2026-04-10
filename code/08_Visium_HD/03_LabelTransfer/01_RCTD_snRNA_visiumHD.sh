#!/bin/bash
#SBATCH --job-name=rctd_array
#SBATCH --array=1-8
#SBATCH --mem=100G
#SBATCH --cpus-per-task=15
#SBATCH --output=logs/rctd_%a.out
#SBATCH --error=logs/rctd_%a.err
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
Rscript 01_RCTD_snRNA_visiumHD.R

echo "********* Job Ends *********"
date
