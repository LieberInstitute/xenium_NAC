#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=01_build_bin_spe
#SBATCH -t 1-0:00:00
#SBATCH -o logs/01_build_bin_spe.log
#SBATCH -e logs/01_build_bin_spe.log

echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOB_ID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Node name: ${HOSTNAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

repo_dir="/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/"

#   Get spatial coordinates as a CSV (from parquet format) for each sample where
#   it doesn't exist
module load visium_hd/1.0
for sample_id in H1-XKYDCP3_A1 H1-XKYDCP3_D1 H1-M3TCP9V_A1 H1-M3TCP9V_D1; do
    spatial_dir=$repo_dir/processed-data/01_spaceranger/$sample_id/outs/binned_outputs/square_008um/spatial

    if [[ ! -f $spatial_dir/tissue_positions.csv ]]; then
        echo "Converting spatial coords to CSV for sample ${sample_id}..."
        
        parquet-tools csv $spatial_dir/tissue_positions.parquet \
            > $spatial_dir/tissue_positions.csv
    fi
done
module unload visium_hd

## Load the R module
module load conda_R/4.5

Rscript 01_build_bin_spe.R

echo "**** Job ends ****"
date
