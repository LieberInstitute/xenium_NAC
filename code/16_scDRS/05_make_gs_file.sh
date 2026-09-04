#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=05_make_gs
#SBATCH --output=logs/05_make_gs.log
#SBATCH --error=logs/05_make_gs.log
#SBATCH --mem=10G
#SBATCH --mail-type=END

set -euo pipefail

echo "********* Job Starts *********"
date

module load conda/3-24.3.0
source activate scdrs_nac

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
GS_DIR=${BASE}/processed-data/16_scDRS/gs

python -u 05_make_gs_file.py

scdrs munge-gs \
    --zscore-file "${GS_DIR}/magma_zscore_matrix.tsv" \
    --out-file "${GS_DIR}/nac_traits.gs" \
    --weight zscore \
    --n-max 1000

echo "********* Job Ends *********"
date
