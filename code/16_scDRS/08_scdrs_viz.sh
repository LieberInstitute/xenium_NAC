#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=08_scdrs_viz
#SBATCH --output=logs/08_scdrs_viz.log
#SBATCH --error=logs/08_scdrs_viz.log
#SBATCH --mem=10G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../..
Rscript code/16_scDRS/08_scdrs_viz.R

echo "********* Job Ends *********"
date
