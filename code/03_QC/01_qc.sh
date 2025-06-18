#!/bin/bash
#
# run_all_qc.sh
# Submit one SLURM job per sample in spe_raw.Rds
#

module load conda_R/4.5

# pull the list of samples
SAMPLES=( $( Rscript -e "
  spe <- readRDS('../../processed-data/02_build_spe/SPEs/spe_raw.Rds');
  cat(unique(as.character(spe\$Sample)), sep=' ')" ) )

for samp in "${SAMPLES[@]}"; do
  sbatch <<EOF
#!/bin/bash
#SBATCH --job-name=qc_${samp}
#SBATCH --output=logs/qc_${samp}.out
#SBATCH --error=logs/qc_${samp}.err
#SBATCH --time=05:00:00
#SBATCH --mem=50G

module load conda_R/4.5

echo \"Starting QC for sample: ${samp} at \$(date)\"
Rscript 01_qc.R "${samp}"
echo \"Finished QC for sample: ${samp} at \$(date)\"
EOF
done
