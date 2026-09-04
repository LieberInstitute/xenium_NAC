#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_01_prep_h5ad
#SBATCH --output=logs/01_prep_h5ad.log
#SBATCH --error=logs/01_prep_h5ad.log
#SBATCH --mem=60G
#SBATCH --mail-type=END


echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../../../
Rscript code/16_scDRS/snRNA/01_prep_h5ad.R

echo "********* Job Ends *********"
date
