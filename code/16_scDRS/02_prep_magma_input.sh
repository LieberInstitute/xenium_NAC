#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=02_prep_magma
#SBATCH --output=logs/02_prep_magma.log
#SBATCH --error=logs/02_prep_magma.log
#SBATCH --mem=10G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date
echo "User: ${USER}; Job id: ${SLURM_JOB_ID}; Hostname: ${HOSTNAME}"

module load conda_R/4.5
cd ../..
Rscript code/16_scDRS/02_prep_magma_input.R

echo "********* Job Ends *********"
date
