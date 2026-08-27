#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=HD_spatial_snRNA
#SBATCH --output=logs/spatial_snRNA.log
#SBATCH --error=logs/spatial_snRNA.log 
#SBATCH --mem=18G
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
Rscript spatial_clusters_0.4_snRNA_spatial_registration.R

echo "********* Job Ends *********"
date
