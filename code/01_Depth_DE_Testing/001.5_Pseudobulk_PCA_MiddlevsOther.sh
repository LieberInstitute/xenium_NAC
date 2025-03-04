#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=001.5_celltype
#SBATCH --output=logs/MiddlevsOther.%a.log
#SBATCH --error=logs/MiddlevsOther.%a.log 
#SBATCH --mem=35G
#SBATCH --array=1-20
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load r_nac
Rscript 001.5_Pseudobulk_PCA_MiddlevsOther.R

echo "********* Job Ends *********"
date
