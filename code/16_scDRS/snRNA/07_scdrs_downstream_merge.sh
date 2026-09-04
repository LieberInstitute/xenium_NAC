#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_07_merge
#SBATCH --output=logs/07_scdrs_merge.log
#SBATCH --error=logs/07_scdrs_merge.log
#SBATCH --mem=30G
#SBATCH --mail-type=END

set -eo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
source activate scdrs_nac

cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/16_scDRS/snRNA
python -u 07_scdrs_downstream_merge.py

echo "********* Job Ends *********"
date
