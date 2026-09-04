#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_06_score
#SBATCH --output=logs/06_scdrs_score.%a.log
#SBATCH --error=logs/06_scdrs_score.%a.log
#SBATCH --mem=40G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-26%6
#SBATCH --mail-type=END

set -eo pipefail

echo "********* Job Starts *********"
date
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

module load conda/3-24.3.0
source activate scdrs_nac

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
MANIFEST=${BASE}/processed-data/16_scDRS/trait_manifest.tsv

TRAIT_KEY=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR == i + 1 {print $1}' "${MANIFEST}")
echo "Trait: ${TRAIT_KEY}"

cd ${BASE}/code/16_scDRS/snRNA
python -u 06_scdrs_score.py --trait-key "${TRAIT_KEY}"

echo "********* Job Ends *********"
date
