#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=01_trait_manifest
#SBATCH --output=logs/01_trait_manifest.log
#SBATCH --error=logs/01_trait_manifest.log
#SBATCH --mem=4G
#SBATCH --mail-type=END

echo "********* Job Starts *********"
date
echo "User: ${USER}; Job id: ${SLURM_JOB_ID}; Hostname: ${HOSTNAME}"

module load conda_R/4.5
cd ../..
Rscript code/16_scDRS/01_trait_manifest.R

echo "********* Job Ends *********"
date
