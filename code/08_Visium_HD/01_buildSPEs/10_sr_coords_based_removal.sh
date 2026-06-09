#!/bin/bash
#SBATCH --mem=10G
#SBATCH --job-name=10_coords
#SBATCH -t 1-0:00:00
#SBATCH -o logs/10_coords.log
#SBATCH -e logs/10_coords.log

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

Rscript 10_sr_coords_based_removal.R

echo "**** Job ends ****"
date
