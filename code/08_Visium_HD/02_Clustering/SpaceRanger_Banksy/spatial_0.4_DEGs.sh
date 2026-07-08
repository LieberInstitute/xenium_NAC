#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=spatial_0.4_DEGs
#SBATCH --output=logs/spatial_0.4_DEGs.log
#SBATCH --error=logs/spatial_0.4_DEGs.log 
#SBATCH --mem=75G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript spatial_0.4_DEGs.R

echo "********* Job Ends *********"
date
