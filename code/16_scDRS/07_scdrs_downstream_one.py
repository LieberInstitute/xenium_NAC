"""scDRS downstream group analysis for ONE trait (array job).

Usage:
    python 07_scdrs_downstream_one.py --trait-key AD
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scdrs
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/HD_Full_Analysis/h5ad/VHD_sfe_counts.h5ad"
SCORE_DIR = BASE / "processed-data/16_scDRS/scores"
OUT_DIR = BASE / "processed-data/16_scDRS/group_analysis"
GROUP_COLS = ["Spatial_Domain", "labels"]


def main(trait_key: str):
    print(f"Loading {H5AD}", flush=True)
    adata = sc.read_h5ad(H5AD)

    # Remove cells outside the NAc (Hypo and Excitatory spatial domains)
    # and Excitatory cell-type labels from any remaining spatial domains
    sd_exclude = ["Hypo", "Excitatory"]
    label_exclude = ["Excitatory"]
    keep = ~adata.obs["Spatial_Domain"].isin(sd_exclude) & ~adata.obs["labels"].isin(label_exclude)
    print(f"Excluding {(~keep).sum()} cells (Spatial_Domain in {sd_exclude} or labels in {label_exclude})", flush=True)
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Same covariates as 06_scdrs_score.py: intercept, log(n_genes detected),
    # and Sample indicators. Without them the kNN graph below is built on data
    # still carrying donor/batch and detection-depth structure, which leaks
    # into the within-group heterogeneity statistic.
    n_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    cov = pd.DataFrame(index=adata.obs.index)
    cov["const"] = 1.0
    cov["n_genes"] = np.log1p(n_genes)
    sample_dummies = pd.get_dummies(
        adata.obs["Sample"], prefix="sample", drop_first=True, dtype=float
    )
    cov = pd.concat([cov, sample_dummies], axis=1)

    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20)
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=20, random_state=4738)

    score_file = SCORE_DIR / f"{trait_key}.full_score.gz"
    print(f"Reading scores: {score_file}", flush=True)
    df_full = pd.read_csv(score_file, sep="\t", index_col=0)
    df_full = df_full.reindex(adata.obs.index)

    print(f"Running downstream analysis for {trait_key}", flush=True)
    stats = scdrs.method.downstream_group_analysis(
        adata=adata, df_full_score=df_full, group_cols=GROUP_COLS
    )

    for g in GROUP_COLS:
        res = stats[g].copy()
        res.insert(0, "trait", trait_key)
        res.insert(1, "group", res.index)
        out_path = OUT_DIR / f"scdrs_group_stats_{g}_{trait_key}.tsv"
        res.to_csv(out_path, sep="\t", index=False)
        print(f"Wrote {out_path}", flush=True)

    # Save per-cell norm_score for this trait
    norm_path = OUT_DIR / f"norm_score_{trait_key}.tsv.gz"
    df_full[["norm_score"]].to_csv(norm_path, sep="\t", compression="gzip")
    print(f"Wrote {norm_path}", flush=True)

    print("Done.", flush=True)
    session_info.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-key", required=True)
    args = parser.parse_args()
    main(args.trait_key)
