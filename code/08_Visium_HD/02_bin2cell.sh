#!/bin/bash
#SBATCH --mem=75G
#SBATCH --job-name=02_bin2cell
#SBATCH -t 1-0:00:00
#SBATCH -o logs/02_bin2cell_%a.log
#SBATCH -e logs/02_bin2cell_%a.log
#SBATCH --array=1-4

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load visium_hd/1.0
module list

python 02_bin2cell.py

echo "**** Job ends ****"
date
