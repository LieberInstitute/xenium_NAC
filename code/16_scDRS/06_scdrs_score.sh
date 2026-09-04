#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=06_scdrs_score
#SBATCH --output=logs/06_scdrs_score.%a.log
#SBATCH --error=logs/06_scdrs_score.%a.log
#SBATCH --mem=50G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-26%6
#SBATCH --mail-type=END

# scDRS scoring, one array task per trait (throttled to 6 concurrent jobs).
# Array range must match trait_manifest.tsv (26 traits).

set -euo pipefail

echo "********* Job Starts *********"
date
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda/3-24.3.0
source activate scdrs_nac

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
MANIFEST=${BASE}/processed-data/16_scDRS/trait_manifest.tsv

TRAIT_KEY=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR == i + 1 {print $1}' "${MANIFEST}")
echo "Trait: ${TRAIT_KEY}"

python -u 06_scdrs_score.py --trait-key "${TRAIT_KEY}"

echo "********* Job Ends *********"
date
