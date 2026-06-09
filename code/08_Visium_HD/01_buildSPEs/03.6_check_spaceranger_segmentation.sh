#!/bin/bash
#SBATCH --mem=44G
#SBATCH --job-name=03.6_seg_check
#SBATCH -t 1-0:00:00
#SBATCH -o logs/03.6_seg_check_%a.log
#SBATCH -e logs/03.6_seg_check_%a.log
#SBATCH --array=1-8

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

python 03.6_check_spaceranger_segmentation.py

echo "**** Job ends ****"
date
