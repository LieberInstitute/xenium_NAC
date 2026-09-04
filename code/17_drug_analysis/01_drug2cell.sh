#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=01_drug2cell
#SBATCH --output=logs/01_drug2cell.log
#SBATCH --error=logs/01_drug2cell.log
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END

set -euo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
source activate scdrs_nac

python -u 01_drug2cell.py

echo "********* Job Ends *********"
date
