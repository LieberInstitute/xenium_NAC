#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=05_project
#SBATCH --output=logs/05_project.log
#SBATCH --error=logs/05_project.log
#SBATCH -p shared
#SBATCH --mem=120G
#SBATCH --time=10-00:00:00

echo "********* Job Starts *********"
date

#load R
module load conda_R/4.5
Rscript 05_project_rat_mouse_xenium.R

echo "********* Job Ends *********"
date
