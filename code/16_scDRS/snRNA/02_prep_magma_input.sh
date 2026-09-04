#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_02_magma_prep
#SBATCH --output=logs/02_prep_magma_input.log
#SBATCH --error=logs/02_prep_magma_input.log
#SBATCH --mem=20G
#SBATCH --mail-type=END

set -eo pipefail

echo "********* Job Starts *********"
date

module load conda_R/4.5
cd ../../..
Rscript code/16_scDRS/snRNA/02_prep_magma_input.R

echo "********* Job Ends *********"
date
