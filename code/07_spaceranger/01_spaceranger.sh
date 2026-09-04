#!/bin/bash
#SBATCH --mem=80G
#SBATCH -n 8
#SBATCH --job-name=NAC-HD_spaceranger
#SBATCH --output=logs/NAC-spaceranger-2603-%a.txt
#SBATCH --array=1-4


echo "**** Job starts ****"
date

echo "**** JHPCE info ****"
echo "User: ${USER}"
echo "Job id: ${SLURM_JOBID}"
echo "Job name: ${SLURM_JOB_NAME}"
echo "Hostname: ${SLURM_NODENAME}"
echo "Task id: ${SLURM_ARRAY_TASK_ID}"

## Load spaceranger
module load spaceranger/4.0.1


## Locate file
SAMPLE=$(awk 'BEGIN {FS="\t"} {print $1}' 03-sample_ids-2603.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMAGE=$(awk 'BEGIN {FS="\t"} {print $2}' 03-sample_ids-2603.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
IMGCYT=$(awk 'BEGIN {FS="\t"} {print $3}' 03-sample_ids-2603.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
SAMPARG=$(awk 'BEGIN {FS="\t"} {print $4}' 03-sample_ids-2603.txt | awk "NR==${SLURM_ARRAY_TASK_ID}")
# SAMPLE=$(awk "NR==${SLURM_ARRAY_TASK_ID}" 02-sample_ids-2601.txt)
echo "Processing sample ${SAMPLE}"
date

## Get slide and area
SLIDE=$(echo ${SAMPLE} | cut -d "_" -f 1)
CAPTUREAREA=$(echo ${SAMPLE} | cut -d "_" -f 2)
SAM=$(paste <(echo ${SLIDE}) <(echo "-") <(echo ${CAPTUREAREA}) -d '')
echo "Slide: ${SLIDE}, capture area: ${CAPTUREAREA}"

## Find FASTQ file path
FASTQPATH=$(ls -d /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/raw-data/fastqs/${SAMPARG})

## Hank from 10x Genomics recommended setting this environment
export NUMBA_NUM_THREADS=1

spaceranger count \
    --id=${SAMPLE} \
    --sample=${SAMPARG} \
    --transcriptome=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A \
    --fastqs=${FASTQPATH} \
    --probe-set=/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/Visium_Human_Transcriptome_Probe_Set_v2.1.0_GRCh38-2024-A.csv \
    --slide=${SLIDE} \
    --area=${CAPTUREAREA} \
    --cytaimage=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/raw-data/HD/${IMGCYT}.tif \
    --image=/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/raw-data/HD/${IMAGE}.tif \
    --create-bam=false \
    --localcores=8 \
    --localmem=64 


## Move output
echo "Moving results to new location"
date
mkdir -p /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/01_spaceranger
mv ${SAMPLE} /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/01_spaceranger

echo "**** Job ends ****"
date
