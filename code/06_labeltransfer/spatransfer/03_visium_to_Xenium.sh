#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=03_vis
#SBATCH --output=logs/03_vis.log
#SBATCH --error=logs/03_vis.log 
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
Rscript 03_visium_to_Xenium.R

echo "********* Job Ends *********"
date
