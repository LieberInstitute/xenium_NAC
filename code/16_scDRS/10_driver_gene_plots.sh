#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=10_driver_gene_plots
#SBATCH --output=logs/10_driver_gene_plots.log
#SBATCH --error=logs/10_driver_gene_plots.log
#SBATCH --mem=20G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../..
Rscript code/16_scDRS/10_driver_gene_plots.R

echo "********* Job Ends *********"
date
