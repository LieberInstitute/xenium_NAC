#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=07_scdrs_down
#SBATCH --output=logs/07_scdrs_downstream.%a.log
#SBATCH --error=logs/07_scdrs_downstream.%a.log
#SBATCH --array=1-26%6
#SBATCH --mem=55G
#SBATCH --cpus-per-task=4
#SBATCH --mail-type=END

set -euo pipefail

echo "********* Job Starts *********"
date
echo "User: ${USER}; Job id: ${SLURM_JOB_ID}; Array task: ${SLURM_ARRAY_TASK_ID}; Hostname: ${HOSTNAME}"

module load conda/3-24.3.0
source activate scdrs_nac

MANIFEST="../../processed-data/16_scDRS/trait_manifest.tsv"
TRAIT_KEY=$(awk -v idx="${SLURM_ARRAY_TASK_ID}" 'NR==idx+1 {print $1}' "${MANIFEST}")
echo "Trait: ${TRAIT_KEY}"

python -u 07_scdrs_downstream_one.py --trait-key "${TRAIT_KEY}"

echo "********* Job Ends *********"
date
