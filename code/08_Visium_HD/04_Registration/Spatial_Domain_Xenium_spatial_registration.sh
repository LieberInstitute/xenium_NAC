#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=Xenium_spatial_domain_snRNA
#SBATCH --output=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/04_Registration/logs/Xenium_spatial_domain_snRNA.log
#SBATCH --error=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/04_Registration/logs/Xenium_spatial_domain_snRNA.log 
#SBATCH --mem=32G
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
cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/
Rscript /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/08_Visium_HD/04_Registration/Spatial_Domain_Xenium_spatial_registration.R

echo "********* Job Ends *********"
date
