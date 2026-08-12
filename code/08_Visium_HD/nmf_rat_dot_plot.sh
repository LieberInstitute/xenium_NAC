#!/bin/bash
#SBATCH --mem=36G
#SBATCH --job-name=nmf_rat_plot
#SBATCH -t 1-0:00:00
#SBATCH -o logs/nmf_rat_plot.log
#SBATCH -e logs/nmf_rat_plot.log

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

Rscript nmf_rat_dot_plot.R

echo "**** Job ends ****"
date
