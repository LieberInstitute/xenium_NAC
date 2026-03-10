#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=09_cell_nonspatial
#SBATCH --output=logs/09_cell_nonspatial_noAGF.log
#SBATCH --error=logs/09_cell_nonspatial_noAGF.log 
#SBATCH --mem=75G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript 09_cell_level_spatial_noAGF.R

echo "********* Job Ends *********"
date
