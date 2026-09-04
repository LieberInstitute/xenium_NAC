#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=04_drug_viz
#SBATCH --output=logs/04_drug_viz.log
#SBATCH --error=logs/04_drug_viz.log
#SBATCH --mem=30G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../..
Rscript code/17_drug_analysis/04_drug_viz.R

echo "********* Job Ends *********"
date
