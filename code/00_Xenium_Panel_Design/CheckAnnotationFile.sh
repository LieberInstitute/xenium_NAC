#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=20GB
#SBATCH --job-name=CheckAnnotationFile
#SBATCH --output=logs/CheckAnnotationFile.log
#SBATCH --error=logs/CheckAnnotationFile.log


echo "********* Job Starts *********"
date


module load r_nac
cd ~/NAc_Xenium_Panel
Rscript Code/CheckAnnotationFile.R

echo "********* Job Ends *********"
date
