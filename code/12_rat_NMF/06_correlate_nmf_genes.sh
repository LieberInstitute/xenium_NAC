#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=06_genes
#SBATCH --output=logs/06_genes.log
#SBATCH --error=logs/6_genes.log
#SBATCH -p shared
#SBATCH --mem=12G
#SBATCH --time=02:00:00

echo "********* Job Starts *********"
date

#load R
module load conda_R/4.5
Rscript 06_correlate_nmf_genes.R

echo "********* Job Ends *********"
date
