#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=00_update
#SBATCH --output=logs/00_update.log
#SBATCH --error=logs/00_update.log
#SBATCH -p shared
#SBATCH --mem=24G
#SBATCH --time=10-00:00:00

echo "********* Job Starts *********"
date

#load R
module load conda_R/4.5
Rscript 00_update_VHD.R

echo "********* Job Ends *********"
date
