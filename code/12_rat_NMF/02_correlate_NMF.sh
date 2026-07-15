#!/bin/bash
#
#SBATCH --mail-type=END
#SBATCH --mail-user=Robert.Phillips@libd.org
#SBATCH --job-name=02_correlate
#SBATCH --output=logs/02_correlate.log
#SBATCH --error=logs/02_correlate.log
#SBATCH -p shared
#SBATCH --mem=16G
#SBATCH --time=10-00:00:00

echo "********* Job Starts *********"
date

#load R
module load conda_R/4.5
Rscript 02_correlate_NMF.R

echo "********* Job Ends *********"
date
