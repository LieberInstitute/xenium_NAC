#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=02_bin_QC
#SBATCH -t 1-0:00:00
#SBATCH -o logs/02_bin_QC.log
#SBATCH -e logs/02_bin_QC.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load the R module
module load conda_R/4.5

Rscript 02_bin_level_QC.R

echo "**** Job ends ****"
date
