#!/bin/bash
#SBATCH -p rocky96
#SBATCH --job-name=03_curated_sets
#SBATCH --output=logs/03_curated_sets.log
#SBATCH --error=logs/03_curated_sets.log
#SBATCH --mem=50G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../..
Rscript code/17_drug_analysis/03_curated_drug_sets.R

echo "********* Job Ends *********"
date
