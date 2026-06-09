#!/bin/bash
#SBATCH --mem=64G
#SBATCH --job-name=08_mergee
#SBATCH -t 1-0:00:00
#SBATCH -o logs/08_merge.log
#SBATCH -e logs/08_merge.log

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

Rscript 08_merge_spaceranger_sfe.R

echo "**** Job ends ****"
date
