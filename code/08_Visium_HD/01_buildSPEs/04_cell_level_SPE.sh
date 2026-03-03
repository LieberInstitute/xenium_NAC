#!/bin/bash
#SBATCH --mem=350G
#SBATCH --job-name=04_cell_spe
#SBATCH -o logs/04_cell_spe.log
#SBATCH -e logs/04_cell_spe.log

module load conda_R/4.5

Rscript 04_cell_level_SPE.R

echo "**** Job ends ****"
date

