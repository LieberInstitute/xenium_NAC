#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=07_cell_nonspatial
#SBATCH --output=logs/07_cell_nonspatial.log
#SBATCH --error=logs/07_cell_nonspatial.log 
#SBATCH --mem=75G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript 07_cell_level_nonSpatial_all.R

echo "********* Job Ends *********"
date
