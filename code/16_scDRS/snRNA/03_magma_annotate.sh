#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_03_magma_annot
#SBATCH --output=logs/03_magma_annotate.log
#SBATCH --error=logs/03_magma_annotate.log
#SBATCH --mem=10G
#SBATCH --mail-type=END

set -eo pipefail

echo "********* Job Starts *********"
date

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
MAGMA_DIR=${BASE}/processed-data/16_scDRS/snRNA/magma
MERGED=${MAGMA_DIR}/ref/1000G.EUR.QC.merged
WINDOW="35,10"

module load magma/1.10

magma \
    --annotate window=${WINDOW} \
    --snp-loc "${MERGED}.bim" \
    --gene-loc "${MAGMA_DIR}/ref/genes_hg19_symbol.gene.loc" \
    --out "${MAGMA_DIR}/annot/snrna_w${WINDOW/,/_}"

echo "********* Job Ends *********"
date
