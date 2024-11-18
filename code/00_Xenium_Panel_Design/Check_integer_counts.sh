#!/bin/bash

#SBATCH -p shared
#SBATCH -c 1
#SBATCH --mem=100GB
#SBATCH --job-name=CheckIntegerCounts
#SBATCH --output=logs/CheckIntegerCounts.log
#SBATCH --error=logs/CheckIntegerCounts.log


echo "********* Job Starts *********"
date


module load r_nac
Rscript Check_integer_counts.R

echo "********* Job Ends *********"
date
