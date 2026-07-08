#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=ct_sfe_pseudo
#SBATCH --output=logs/ct_sfe_pseudo.log
#SBATCH --error=logs/ct_sfe_pseudo.log 
#SBATCH --mem=24G
#SBATCH --cpus-per-task=1
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END
#SBATCH --mail-user=robert.phillips@libd.org

echo "********* Job Starts *********"
date

module load conda_R/4.5
Rscript make_sfe_pseudo_celltype.R

echo "********* Job Ends *********"
date
