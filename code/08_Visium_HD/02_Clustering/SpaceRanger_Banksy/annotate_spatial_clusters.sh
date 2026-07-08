#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=spatial_anno
#SBATCH --output=logs/spatial_anno.log
#SBATCH --error=logs/spatial_anno.log 
#SBATCH --mem=15G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript annotate_spatial_clusters.R

echo "********* Job Ends *********"
date
