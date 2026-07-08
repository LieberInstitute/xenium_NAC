#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=05_anno
#SBATCH --output=logs/05_anno.log
#SBATCH --error=logs/05_anno.log 
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript 05_cell_level_spatial_annotation.R

echo "********* Job Ends *********"
date
