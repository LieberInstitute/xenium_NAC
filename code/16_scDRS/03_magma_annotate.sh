#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=03_magma_annotate
#SBATCH --output=logs/03_magma_annotate.log
#SBATCH --error=logs/03_magma_annotate.log
#SBATCH --mem=10G
#SBATCH --mail-type=END

# One-time steps:
#   a) merge the per-chromosome 1000G EUR PLINK files into a single bfile
#   b) run the MAGMA annotation step (SNP -> gene assignment)
#
# Gene window: 35 kb upstream / 10 kb downstream (Bryois et al. convention,
# widely used for cell-type GWAS enrichment). NOTE: the gene loc file has no
# strand information, so MAGMA assumes + strand; for minus-strand genes the
# asymmetric window is flipped. This is a minor, commonly accepted
# approximation; use a symmetric window (e.g. window=20) if preferred.

set -euo pipefail

echo "********* Job Starts *********"
date

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
MAGMA_DIR=${BASE}/processed-data/16_scDRS/magma
REF1000G=/dcs04/lieber/shared/statsgen/LDSC/base/referencefiles/1000G_EUR_Phase3_plink
MERGED=${MAGMA_DIR}/ref/1000G.EUR.QC.merged
WINDOW="35,10"

module load plink/1.90b
module load magma/1.10

## a) Merge per-chromosome PLINK files (skip if already done)
if [ ! -f "${MERGED}.bed" ]; then
    MERGE_LIST=${MAGMA_DIR}/ref/merge_list.txt
    : > "${MERGE_LIST}"
    for chr in $(seq 2 22); do
        echo "${REF1000G}/1000G.EUR.QC.${chr}" >> "${MERGE_LIST}"
    done
    plink --bfile "${REF1000G}/1000G.EUR.QC.1" \
        --merge-list "${MERGE_LIST}" \
        --make-bed \
        --out "${MERGED}"
else
    echo "Merged bfile already exists, skipping merge."
fi

## b) MAGMA annotation
magma \
    --annotate window=${WINDOW} \
    --snp-loc "${MERGED}.bim" \
    --gene-loc "${MAGMA_DIR}/ref/genes_hg19_symbol.gene.loc" \
    --out "${MAGMA_DIR}/annot/nac_w${WINDOW/,/_}"

echo "********* Job Ends *********"
date
