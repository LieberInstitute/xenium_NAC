#!/bin/bash
#SBATCH --mem=24G
#SBATCH --job-name=cc_cor
#SBATCH -t 1-0:00:00
#SBATCH -o logs/cc_cor.log
#SBATCH -e logs/cc_cor.log

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

Rscript cross_species_correlation.R

echo "**** Job ends ****"
date
