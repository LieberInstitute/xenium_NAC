#!/bin/bash
#SBATCH --mem=128G
#SBATCH -n 4
#SBATCH --time=3-00:00:00
#SBATCH --job-name=01-NAC-AP_create_samui
#SBATCH -o logs/rerun_samui_%a.txt
#SBATCH -e logs/rerun_samui_%a.txt
#SBATCH --array=1-4


# %4

echo "**** Job starts ****"
date


echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

donor=$(awk "NR==${SLURM_ARRAY_TASK_ID}" /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/14_samui/visHD/rerun_sample_ids-visium.txt)
echo "Processing sample ${donor}"
date


module load samui/1.0.0-next.45
python 02-create_samui.py $donor

echo "**** Job ends ****"
date

