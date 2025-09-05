#!/bin/bash
#SBATCH --job-name=banksy_grid
#SBATCH --mem=250G
#SBATCH --array=0-7

##Reasoning:
#Banksy vignette says 0.2 for cell type and 0.8 for spatial information.
#We are going to run with 2 different values of lambda to see how clusters change with same res value (0.5).
#Also going to change k_geom
#R script will also save the embeddings, so we can cluster with different algorithms and compare.

# Define lambda and k_geom grids
LAMBDAS=(0.2 1.0)
KGEOMS=(10 15 25 30)

# Get this task’s index
idx=$SLURM_ARRAY_TASK_ID

# Decode it into (i,j)
i=$(( idx / ${#KGEOMS[@]} ))   # which lambda
j=$(( idx % ${#KGEOMS[@]} ))   # which k_geom

# Pick the actual values
L="${LAMBDAS[$i]}"
K="${KGEOMS[$j]}"

LOGDIR="logs"
LOGFILE="${LOGDIR}/banksy_lambda_${L}_kgeom_${K}.log"

# Redirect *all* stdout+stderr from this job into LOGFILE
exec >"$LOGFILE" 2>&1

echo "********* Job Starts *********"
date
echo "**** SLURM info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"
echo "Lambda: $L"
echo "k_geom: $K"

module load conda_R/4.5

# Run  R script
Rscript Banksy_clustering_grid.R "$L" "$K"

echo "********* Job Ends *********"
date
