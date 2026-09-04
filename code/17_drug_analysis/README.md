# 17_drug_analysis: mapping drugs of abuse, psychiatric drugs, and other compounds to Visium-HD spatial domains and cell types

Two complementary approaches on the SpaceRanger-segmented Visium-HD cells
(439k cells; groupings: `Spatial_Domain` [12 domains], `labels`
[snRNA-transferred cell types]):

1. **drug2cell** (Kanemaru et al., Nature 2023): unbiased per-cell scoring of
   every ChEMBL drug's target genes, then Wilcoxon enrichment per domain /
   cell type. A curated panel (opioids, stimulants, nicotine, cannabinoids,
   antipsychotics, antidepressants, mood stabilizers, benzodiazepines,
   dissociatives/psychedelics) is highlighted by name matching.
2. **Curated drug-target gene sets** (R): interpretable receptor/transporter
   target sets per drug class, scored as mean logcounts (z-scaled), tested
   with Kruskal-Wallis + group-vs-rest Wilcoxon (BH FDR).

## Run order

To run everything via a SLURM dependency chain (job IDs logged to
`logs/run_all_<timestamp>.log`):

```bash
bash run_all.sh
```

Or step by step:

| Step | Script | Where | Notes |
|---|---|---|---|
| 0 | `../16_scDRS/00_setup_scdrs_env.sh` | login node | shared `scdrs_nac` conda env (includes drug2cell) |
| 1 | `sbatch 01_drug2cell.sh` | SLURM | scores ChEMBL drugs on `VHD_sfe_logcounts.h5ad` |
| 2 | `sbatch 02_drug2cell_enrichment.sh` | SLURM | Wilcoxon per domain and per cell type |
| 3 | `sbatch 03_curated_drug_sets.sh` | SLURM | curated set scoring + tests (independent of 1–2) |
| 4 | `sbatch 06_drug_viz.sh` | SLURM | combined dotplots / heatmaps / tables |

## Outputs

- `processed-data/17_drug_analysis/drug2cell/` — score h5ad + enrichment TSVs
- `processed-data/17_drug_analysis/curated_sets/` — per-cell scores + tests
- `plots/17_drug_analysis/` — figures

## Caveats

- Both approaches are *descriptive target/signature enrichment*, not
  evidence of drug effect; they indicate where the molecular substrates of a
  drug are concentrated.
- drug2cell scores depend on ChEMBL target annotation completeness; lithium
  and ethanol have poorly defined ChEMBL targets, which is why the curated
  set approach complements it.
