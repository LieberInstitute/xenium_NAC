#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=xenium_visium
#SBATCH --output=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/05_clustering/Banksy/logs/xenium_visium.log
#SBATCH --error=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/05_clustering/Banksy/logs/xenium_visium.log 
#SBATCH --mem=12G
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
cd  /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/ 
Rscript /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/05_clustering/Banksy/xenium_visium_spatial_registration_all_celltypes.R

echo "********* Job Ends *********"
date
