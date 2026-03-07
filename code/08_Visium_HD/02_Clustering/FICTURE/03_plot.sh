#!/bin/bash
#SBATCH --mem=32G
#SBATCH --job-name=03_FICTURE_plot
#SBATCH -t 4:00:00
#SBATCH -o logs/03_FICTURE_plot_%a.log
#SBATCH -e logs/03_FICTURE_plot_%a.log
#SBATCH --array=1-4%4

set -e

echo "**** Job starts ****"
date

module load visium_hd/1.0

repo_dir=$(git rev-parse --show-toplevel)
all_samples=(H1-XKYDCP3_A1 H1-XKYDCP3_D1 H1-M3TCP9V_A1 H1-M3TCP9V_D1)
this_sample=${all_samples[$(($SLURM_ARRAY_TASK_ID - 1))]}

for nf in 8 12 16; do

echo "--------------------------------------------------------------"
echo "Plotting nFactor=${nf} for ${this_sample}"
echo "--------------------------------------------------------------"

out_dir1=$repo_dir/processed-data/HD_Full_Analysis/FICTURE/outputs/$this_sample/nF_${nf}/analysis/nF${nf}.d_12/
out_dir2=$repo_dir/processed-data/HD_Full_Analysis/FICTURE/outputs/$this_sample/nF_${nf}/analysis/nF${nf}.d_12/figure
plot_output=$repo_dir/plots/HD_Full_Analysis/FICTURE/$this_sample/nF_${nf}_${this_sample}_clusters.png

  # Find the pixel and color files
pixel_file=$(ls ${out_dir1}/*.pixel.sorted.tsv.gz 2>/dev/null | head -1)
color_file=$(ls ${out_dir2}/*.rgb.tsv 2>/dev/null | head -1)

if [[ -z "$pixel_file" || -z "$color_file" ]]; then
echo "WARNING: Missing pixel or color file for nF=${nf}, skipping"
echo "  pixel_file: $pixel_file"
echo "  color_file: $color_file"
fi

echo "  pixel_file: $pixel_file"
echo "  color_file: $color_file"

ficture plot_pixel_full \
--input $pixel_file \
--color_table $color_file \
--plot_um_per_pixel 2 \
--output $plot_output

echo "Plot saved: $plot_output "
date

done

echo "**** Job ends ****"
date
