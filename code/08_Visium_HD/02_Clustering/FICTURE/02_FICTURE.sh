#!/bin/bash
#SBATCH --mem=150G
#SBATCH --job-name=02_FICTURE
#SBATCH -t 2-0:00:00
#SBATCH -o logs/02_FICTURE_%a.log
#SBATCH -e logs/02_FICTURE_%a.log
#SBATCH --array=1-4%4

set -e

echo "**** Job starts ****"
date

module load visium_hd/1.0

repo_dir=$(git rev-parse --show-toplevel)
all_samples=(H1-XKYDCP3_A1 H1-XKYDCP3_D1 H1-M3TCP9V_A1 H1-M3TCP9V_D1)
this_sample=${all_samples[$(($SLURM_ARRAY_TASK_ID - 1))]}

#   Path definitions
in_dir=$repo_dir/processed-data/HD_Full_Analysis/FICTURE/inputs/$this_sample

#   Loop over multiple K values
for nf in 8 12 16; do

    echo "--------------------------------------------------------------"
    echo "Running FICTURE with nFactor=${nf} for ${this_sample}"
    echo "--------------------------------------------------------------"
    date

    out_dir=$repo_dir/processed-data/HD_Full_Analysis/FICTURE/outputs/$this_sample/nF_${nf}
    mkdir -p $out_dir

    ficture run_together \
        --in-tsv $in_dir/transcripts.sorted.tsv.gz \
        --in-minmax $in_dir/minmax.tsv \
        --out-dir $out_dir \
        --mu-scale 1 \
        --major-axis X \
        --n-factor $nf \
        --fractional-count 0 \
        --all

    echo "Finished nFactor=${nf} for ${this_sample}"
    date

done

echo "**** Job ends ****"
date
