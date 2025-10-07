#!/bin/bash

#SBATCH -p gpu
#SBATCH --job-name=01_STAligner
#SBATCH --mem=500G
#SBATCH --time=24:00:00
#SBATCH --gres=gpu:tesh100:1
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jyao37@jhmi.edu
#SBATCH --array=0-2

set -euo pipefail

# map array idx → donor
DONORS=("Br6660" "Br6436" "All")
DONOR=${DONORS[$SLURM_ARRAY_TASK_ID]}

# fixed radius
RADIUS=0

# ensure log dir
LOG_DIR=logs
mkdir -p "${LOG_DIR}"

# now redirect all stdout/stderr to donor‑specific files
exec > "${LOG_DIR}/01_STAligner_${DONOR}.out" \
     2> "${LOG_DIR}/01_STAligner_${DONOR}.err"

echo "********* Job Starts *********"
date
echo "Donor:  ${DONOR}"
echo "Radius: ${RADIUS}"
echo "Node:   ${HOSTNAME}"
echo

module load anaconda
conda run -n STAligner --no-capture-output \
      python 01_STAligner.py "${DONOR}" "${RADIUS}"

echo
echo "********* Job Ends *********"
date
