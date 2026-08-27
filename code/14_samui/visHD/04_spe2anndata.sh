#!/bin/bash
#SBATCH --mem=50G
#SBATCH --job-name=spe2anndata-xen
#SBATCH --array=1
#SBATCH -o logs/spe2anndata_%a.txt

echo -en '\n'
echo "**** Job starts ****"
date
echo -en '\n'

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.4
Rscript 03-spe_to_anndata.R

mkdir -p ./logs/

#cat ./slurm-${SLURM_JOBID}*.out >> ./logs/01_spe2anndata.txt
#echo -en '\n'
#rm ./slurm-${SLURM_JOBID}*.out

echo "**** Job ends ****"
date