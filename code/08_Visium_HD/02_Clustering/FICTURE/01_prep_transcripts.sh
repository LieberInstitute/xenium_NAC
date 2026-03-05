#!/bin/bash
#SBATCH --mem=10G
#SBATCH --job-name=01_transcripts
#SBATCH -t 1-0:00:00
#SBATCH -o logs/01_transcripts_%a.log
#SBATCH -e logs/01_transcripts_%a.log
#SBATCH --array=1-4%2

#exit immidately if non-zero exit
set -e

echo "**** Job starts ****"
date

#Load the module
module load visium_hd/1.0

#Set variables
repo_dir=$(git rev-parse --show-toplevel)
temp_dir=$MYSCRATCH
all_samples=(H1-XKYDCP3_A1 H1-XKYDCP3_D1 H1-M3TCP9V_A1 H1-M3TCP9V_D1)

this_sample=${all_samples[$(($SLURM_ARRAY_TASK_ID - 1))]}
data_dir=$repo_dir/processed-data/01_spaceranger/$this_sample/outs/binned_outputs/square_002um #performing on 2um bins
out_dir=$repo_dir/processed-data/HD_Full_Analysis/FICTURE/inputs/$this_sample

mkdir -p $out_dir

#   Get spatial coordinates as a CSV
echo "Converting spatial coords to CSV..."
parquet-tools csv $data_dir/spatial/tissue_positions.parquet \
    | gzip -c > $temp_dir/${this_sample}_tissue_positions.csv.gz

#Pull microns per pixel
 #Split by : and pull second part; then remove(transform) and cut away the comman and space
microns_per_pixel=$(
    grep microns_per_pixel $data_dir/spatial/scalefactors_json.json \
        | cut -d ":" -f 2 \
        | tr -d ", " 
)


#Generate all input needed for FICTURE
spatula convert-sge \
    --in-sge $data_dir/raw_feature_bc_matrix \
    --pos $temp_dir/${this_sample}_tissue_positions.csv.gz \
    --units-per-um $(python -c "print(1/${microns_per_pixel})") \
    --colnames-count Count \
    --out-tsv $out_dir \
    --icols-mtx 1

## Sort by the X-coordinate
(gzip -cd $out_dir/transcripts.unsorted.tsv.gz \
    | head -1; gzip -cd $out_dir/transcripts.unsorted.tsv.gz \
    | tail -n +2 | sort -S 1G -gk1) \
    | gzip -c > $out_dir/transcripts.sorted.tsv.gz

echo "**** Job ends ****"
date
