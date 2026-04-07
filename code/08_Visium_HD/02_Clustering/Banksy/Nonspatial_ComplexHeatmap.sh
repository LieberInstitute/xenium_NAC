#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=CH_nonspatial
#SBATCH --output=logs/nonspatial_complexheatmap.log
#SBATCH --error=logs/nonspatial_complexheatmap.log 
#SBATCH --mem=100G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript Nonspatial_ComplexHeatmap.R

echo "********* Job Ends *********"
date
