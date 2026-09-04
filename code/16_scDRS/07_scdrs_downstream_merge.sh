#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=07_scdrs_merge
#SBATCH --output=logs/07_scdrs_merge.log
#SBATCH --error=logs/07_scdrs_merge.log
#SBATCH --mem=30G
#SBATCH --mail-type=END

set -euo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
source activate scdrs_nac

python -u 07_scdrs_downstream_merge.py

echo "********* Job Ends *********"
date
