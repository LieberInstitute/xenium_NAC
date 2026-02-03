#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=01_snRNA_to_Xenium
#SBATCH --output=logs/01_snRNA_to_Xenium.log
#SBATCH --error=logs/01_snRNA_to_Xenium.log 
#SBATCH --mem=250G
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --time=7-00:00:00


echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
Rscript 01_snRNA_to_Xenium.R

echo "********* Job Ends *********"
date
