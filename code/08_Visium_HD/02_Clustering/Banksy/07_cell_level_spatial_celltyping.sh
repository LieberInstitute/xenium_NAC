#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=07_cell
#SBATCH --output=logs/07_cell.log
#SBATCH --error=logs/07_cell.log 
#SBATCH --mem=400G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript 07_cell_level_spatial_celltyping.R

echo "********* Job Ends *********"
date
