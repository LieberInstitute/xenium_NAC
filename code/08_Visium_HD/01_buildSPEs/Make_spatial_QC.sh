#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=make_qc
#SBATCH -t 1-0:00:00
#SBATCH -o /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/01_buildSPEs/logs/make_spatial_qc.log
#SBATCH -e /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/01_buildSPEs/logs/make_spatial_qc.log

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
cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
Rscript /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/01_buildSPEs/Make_spatial_QC.R

echo "**** Job ends ****"
date
