#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=CH_spatial
#SBATCH --output=logs/spatial_complexheatmap.log
#SBATCH --error=logs/spatial_complexheatmap.log 
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript Spatial_ComplexHeatmap.R

echo "********* Job Ends *********"
date
