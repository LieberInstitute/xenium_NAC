#!/bin/bash
#SBATCH --job-name=hvg_block
#SBATCH --output=logs/hvg_block.log
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00

module load conda_R/4.5
Rscript block_check.R
