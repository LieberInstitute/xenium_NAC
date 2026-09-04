#!/bin/bash
#SBATCH -p shared
#SBATCH --job-name=snRNA_04_magma_genes
#SBATCH --output=logs/04_magma_genes.%a.log
#SBATCH --error=logs/04_magma_genes.%a.log
#SBATCH --mem=5G
#SBATCH --array=1-26
#SBATCH --mail-type=END

set -eo pipefail

echo "********* Job Starts *********"
date
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

BASE=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
MAGMA_DIR=${BASE}/processed-data/16_scDRS/snRNA/magma
MANIFEST=${BASE}/processed-data/16_scDRS/trait_manifest.tsv
MERGED=${MAGMA_DIR}/ref/1000G.EUR.QC.merged
ANNOT=${MAGMA_DIR}/annot/snrna_w35_10.genes.annot

module load magma/1.10

TRAIT_KEY=$(awk -v i="${SLURM_ARRAY_TASK_ID}" 'NR == i + 1 {print $1}' "${MANIFEST}")
echo "Trait: ${TRAIT_KEY}"

# snp-wise=mean is MAGMA's default, but state it explicitly so the model is
# recorded in the flags block of every log and cannot drift if a future
# MAGMA version changes its default.
magma \
    --bfile "${MERGED}" \
    --pval "${MAGMA_DIR}/pval/${TRAIT_KEY}.pval.txt" ncol=N \
    --gene-annot "${ANNOT}" \
    --gene-model snp-wise=mean \
    --out "${MAGMA_DIR}/genes/${TRAIT_KEY}"

echo "********* Job Ends *********"
date
