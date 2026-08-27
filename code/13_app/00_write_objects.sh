#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=00_write
#SBATCH --output=logs/00_write.log
#SBATCH --error=logs/00_write.log
#SBATCH -p shared
#SBATCH --mem=24G
#SBATCH --time=10-00:00:00

echo "********* Job Starts *********"
date

#load R
module load conda_R/4.5
Rscript 00_write_objects.R

echo "********* Job Ends *********"
date
