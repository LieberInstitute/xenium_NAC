#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=09_risk_gene_drivers
#SBATCH --output=logs/09_risk_gene_drivers.log
#SBATCH --error=logs/09_risk_gene_drivers.log
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END

# Without this a Python traceback still exits 0 (the trailing echo/date
# succeed), so SLURM reports COMPLETED and any afterok dependency proceeds on
# missing output.
set -euo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
source activate scdrs_nac
cd ../..
python code/16_scDRS/09_risk_gene_drivers.py

echo "********* Job Ends *********"
date
