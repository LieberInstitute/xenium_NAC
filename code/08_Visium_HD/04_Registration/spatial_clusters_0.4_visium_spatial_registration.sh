#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=HD_spatial_visium
#SBATCH --output=logs/spatial_visium.log
#SBATCH --error=logs/spatial_visium.log 
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda_R/4.5
Rscript spatial_clusters_0.4_visium_spatial_registration.R

echo "********* Job Ends *********"
date
