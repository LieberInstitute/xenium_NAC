#!/bin/bash

#SBATCH -p shared
#SBATCH --job-name=01_harmony
#SBATCH --output=logs/01_harmony.out
#SBATCH --error=logs/01_harmony.err
#SBATCH --mem=200G
#SBATCH --time=24:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jyao37@jhmi.edu

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID:-N/A}"

module load anaconda
conda run -n STAligner python 01_harmony.py

echo "********* Job Ends *********"
date
