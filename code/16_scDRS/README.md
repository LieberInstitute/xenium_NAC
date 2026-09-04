# 16_scDRS: mapping GWAS genetic susceptibility to Visium-HD spatial domains and cell types

Uses [scDRS](https://github.com/martinjzhang/scDRS) (Zhang et al., Nat Genet 2022)
with MAGMA gene-level statistics to score 439k SpaceRanger-segmented Visium-HD
cells (`VHD_sfe_counts.h5ad`) for 26 GWAS traits, then tests association of
spatial domains (`Spatial_Domain`, 12 domains) and cell types (`labels`,
snRNA-transferred) with each trait.

## Inputs

- GWAS: LDSC-munged sumstats at `/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain`
  (SNP, N, Z format; P derived from Z)
- Expression: `processed-data/HD_Full_Analysis/h5ad/VHD_sfe_counts.h5ad`
- MAGMA reference: 1000G EUR Phase 3 PLINK (hg19), gene coordinates from
  `gene_meta_hg19.txt` (symbol-keyed so MAGMA genes match h5ad var names)

## Run order

To run the whole pipeline unattended (SLURM dependency chain, job IDs logged
to `logs/run_all_<timestamp>.log`):

```bash
bash 00_setup_scdrs_env.sh   # one-time, login node
bash run_all.sh
```

Or step by step:

| Step | Script | Where | Notes |
|---|---|---|---|
| 0 | `00_setup_scdrs_env.sh` | login node (internet) | builds `scdrs_nac` conda env (shared with 17_drug_analysis) |
| 1–2 | `sbatch 02_prep_magma_input.sh` | SLURM | trait manifest + MAGMA pval/gene.loc files |
| 3 | `sbatch 03_magma_annotate.sh` | SLURM | merge 1000G PLINK + MAGMA annotation (window 35/10 kb) |
| 4 | `sbatch 04_magma_genes.sh` | SLURM array 1-26 | per-trait MAGMA gene analysis |
| 5 | `sbatch 05_make_gs_file.sh` | SLURM | z-score matrix → `scdrs munge-gs` (top 1000 genes, zscore weights) |
| 6 | `sbatch 06_scdrs_score.sh` | SLURM array 1-26%6 | per-trait scDRS scoring (1000 MC control sets; cov: const, log n_genes, Sample) |
| 7 | `sbatch 07_scdrs_downstream.sh` | SLURM | group association + heterogeneity for both groupings |
| 8 | `sbatch 08_scdrs_viz.sh` | SLURM | heatmaps, spatial score maps, FDR tables |

If the trait manifest changes, update the `--array` ranges in steps 4 and 6.

## Outputs

- `processed-data/16_scDRS/magma/genes/*.genes.out` — MAGMA gene z-scores
- `processed-data/16_scDRS/gs/nac_traits.gs` — scDRS gene-set file
- `processed-data/16_scDRS/scores/*.{score,full_score}.gz` — per-cell scores
- `processed-data/16_scDRS/group_analysis/scdrs_group_stats_{Spatial_Domain,labels}.tsv`
- `processed-data/16_scDRS/tables/*_fdr.csv` — with BH FDR across trait × group
- `plots/16_scDRS/` — heatmaps + spatial maps

## Interpretation notes

- `assoc_mcz`/`assoc_mcp`: Monte Carlo group-trait association (main result).
- `hetero_mcz`/`hetero_mcp`: within-group heterogeneity — significant values
  mean risk is concentrated in a subpopulation of the domain/cell type.
- Negative-control traits (HEIGHT, BMI, T2D, blood counts) should show little
  association; if they don't, treat brain-trait results with caution.
- Xenium (366-gene panel) is intentionally excluded: scDRS requires the full
  transcriptome for matched control gene sets.
