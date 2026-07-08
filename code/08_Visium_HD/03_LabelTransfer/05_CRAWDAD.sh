#!/bin/bash
#SBATCH --job-name=crawdad
#SBATCH --output=logs/crawdad_%a.log
#SBATCH --error=logs/crawdad_%a.log
#SBATCH --array=1-8
#SBATCH --cpus-per-task=8
#SBATCH --mem=250G
#SBATCH --time=24:00:00

module load conda_R/4.5

Rscript 05_CRAWDAD.R
