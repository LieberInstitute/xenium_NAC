#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=02_d2c_enrich
#SBATCH --output=logs/02_d2c_enrich.log
#SBATCH --error=logs/02_d2c_enrich.log
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END

# Without this a traceback still exits 0 (the trailing echo/date succeed), so
# SLURM reports COMPLETED and afterok dependencies proceed on missing output.
set -euo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
conda activate scdrs_nac

python 02_drug2cell_enrichment.py

echo "********* Job Ends *********"
date
