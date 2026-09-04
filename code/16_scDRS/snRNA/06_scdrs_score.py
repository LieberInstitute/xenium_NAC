"""Compute scDRS disease scores for one GWAS trait on snRNA-seq cells.

Input: snRNA_counts.h5ad (obs carries `Sample` and `CellType.Final`).

Covariates: intercept, log(n_genes detected), and Sample indicator variables.

Usage:
    python 06_scdrs_score.py --trait-key OUD
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scdrs
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/16_scDRS/snRNA/h5ad/snRNA_counts.h5ad"
GS_FILE = BASE / "processed-data/16_scDRS/snRNA/gs/snrna_traits.gs"
OUT_DIR = BASE / "processed-data/16_scDRS/snRNA/scores"
N_CTRL = 1000


def main(trait_key: str) -> None:
    print(f"Loading {H5AD}")
    adata = sc.read_h5ad(H5AD)
    print(adata)

    # Remove Excitatory and Neuron_Ambig cell types
    ct_exclude = ["Excitatory", "Neuron_Ambig"]
    keep = ~adata.obs["CellType.Final"].isin(ct_exclude)
    print(f"Excluding {(~keep).sum()} cells (CellType.Final in {ct_exclude})")
    adata = adata[keep].copy()

    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    n_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    cov = pd.DataFrame(index=adata.obs.index)
    cov["const"] = 1.0
    cov["n_genes"] = np.log1p(n_genes)
    sample_dummies = pd.get_dummies(
        adata.obs["Sample"], prefix="sample", drop_first=True, dtype=float
    )
    cov = pd.concat([cov, sample_dummies], axis=1)

    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20)

    gs = scdrs.util.load_gs(
        str(GS_FILE),
        src_species="human",
        dst_species="human",
        to_intersect=adata.var_names,
    )
    if trait_key not in gs:
        raise KeyError(f"{trait_key} not in {GS_FILE}; available: {list(gs)}")
    gene_list, gene_weights = gs[trait_key]
    print(f"{trait_key}: {len(gene_list)} genes overlap the h5ad")

    df_score = scdrs.score_cell(
        adata,
        gene_list=gene_list,
        gene_weight=gene_weights,
        ctrl_match_key="mean_var",
        n_ctrl=N_CTRL,
        weight_opt="vs",
        return_ctrl_raw_score=False,
        return_ctrl_norm_score=True,
    )

    keep = [c for c in df_score.columns if not c.startswith("ctrl_norm_score")]
    df_score[keep].to_csv(
        OUT_DIR / f"{trait_key}.score.gz", sep="\t", compression="gzip"
    )
    df_score.to_csv(
        OUT_DIR / f"{trait_key}.full_score.gz", sep="\t", compression="gzip"
    )
    print(f"Wrote scores for {trait_key} to {OUT_DIR}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--trait-key", required=True)
    args = parser.parse_args()
    main(args.trait_key)

    session_info.show()
