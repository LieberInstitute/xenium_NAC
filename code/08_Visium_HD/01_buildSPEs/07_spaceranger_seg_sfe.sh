#!/bin/bash
#SBATCH --mem=40G
#SBATCH --job-name=07_build
#SBATCH -t 1-0:00:00
#SBATCH -o logs/07_build_%a.log
#SBATCH -e logs/07_build_%a.log
#SBATCH --array=1-8

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
module list

Rscript 07_spaceranger_seg_sfe.R

echo "**** Job ends ****"
date
