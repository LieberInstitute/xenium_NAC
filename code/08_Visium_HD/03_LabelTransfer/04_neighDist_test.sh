#!/bin/bash
#SBATCH --job-name=crawdad_neighDist
#SBATCH --output=logs/crawdad_neighDist_%a.log
#SBATCH --error=logs/crawdad_neighDist_%a.log
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --time=12:00:00

module load conda_R/4.5

Rscript 04_neighDist_test.R
