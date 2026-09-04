"""Top risk-gene drivers per group (spatial domain and cell-type label).

For each trait, restrict to that trait's scDRS/MAGMA risk genes and, within each
group (a Spatial_Domain or a `labels` cell type), compute the Spearman
correlation between each risk gene's log-normalized expression and the per-cell
scDRS normalized score across the cells in that group. The top-ranked genes are
those whose expression most tracks disease risk *within that context*, i.e. the
genes driving genetic risk in that domain / cell type.

Inputs
------
- VHD_sfe_counts.h5ad         : per-cell counts (same object used for scoring)
- nac_traits.gs               : per-trait risk gene sets + MAGMA z-score weights
- norm_scores_all_traits.tsv.gz : per-cell scDRS normalized score, one col/trait

Preprocessing mirrors 06_scdrs_score.py exactly: normalize_total 1e4, log1p,
then the same covariate correction (intercept, log n_genes, Sample).

NOTE on how scDRS handles covariates: for a sparse .X, scdrs.preprocess does
NOT overwrite adata.X -- it would densify a 400k x 18k matrix. Instead it fits
the regression and stores COV_MAT / COV_BETA / COV_GENE_MEAN in
adata.uns["SCDRS_PARAM"], defining the corrected data implicitly as

    CORRECTED_X = adata.X + COV_MAT @ COV_BETA + COV_GENE_MEAN

So the correction must be applied explicitly (apply_cov_correction below); just
calling preprocess leaves adata.X untouched and the correction silently
inactive. We apply it only to each trait's ~1000 risk genes, which is cheap.

Outputs (processed-data/16_scDRS/risk_gene_drivers/)
----------------------------------------------------
- risk_gene_drivers_Spatial_Domain.tsv.gz  : all trait x domain x gene rows
- risk_gene_drivers_labels.tsv.gz          : all trait x label x gene rows
- top20_risk_gene_drivers_Spatial_Domain.tsv
- top20_risk_gene_drivers_labels.tsv
- risk_gene_expression_Spatial_Domain.tsv.gz : per trait x group x gene mean
  covariate-corrected expression, its z-score across groups, and the fraction
  of cells with a nonzero count (frac_expressing, on 0-1). This is a *level*
  summary, complementing the correlation above: it answers "where is this risk
  gene expressed", not "where does it track the score".
- risk_gene_expression_labels.tsv.gz

Usage: sbatch 09_risk_gene_drivers.sh
"""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scdrs
from scipy.stats import rankdata

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/HD_Full_Analysis/h5ad/VHD_sfe_counts.h5ad"
GS_FILE = BASE / "processed-data/16_scDRS/gs/nac_traits.gs"
SCORES = BASE / "processed-data/16_scDRS/group_analysis/norm_scores_all_traits.tsv.gz"
OUT_DIR = BASE / "processed-data/16_scDRS/risk_gene_drivers"

GROUP_COLS = ["Spatial_Domain", "labels"]
TOP_N = 20
# Skip groups with too few cells for a stable correlation.
MIN_CELLS = 30


def load_expression():
    """Load and normalize expression exactly as 06_scdrs_score.py did.

    Cell subsetting (excluded domains/labels) is handled downstream by
    intersecting with the scored cells, which are authoritative — the
    norm_scores file contains only the cells scDRS actually scored.
    """
    adata = sc.read_h5ad(H5AD)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def fit_cov_correction(adata):
    """Fit the covariate correction 06_scdrs_score.py used before scoring.

    Intercept, log(number of genes detected) and Sample indicators, with the
    same 20x20 mean/variance binning. For a sparse .X this only *stores* the
    correction in adata.uns["SCDRS_PARAM"]; apply_cov_correction() applies it.

    Must be called AFTER subsetting to the scored cells, so the regression is
    fit on the same cell set scDRS used.
    """
    n_genes = np.asarray((adata.X > 0).sum(axis=1)).ravel()
    cov = pd.DataFrame(index=adata.obs.index)
    cov["const"] = 1.0
    cov["n_genes"] = np.log1p(n_genes)
    sample_dummies = pd.get_dummies(
        adata.obs["Sample"], prefix="sample", drop_first=True, dtype=float
    )
    cov = pd.concat([cov, sample_dummies], axis=1)
    print(
        f"Fitting covariate correction on {cov.shape[1]} covariates "
        f"(const, n_genes, {sample_dummies.shape[1]} sample dummies) ...",
        flush=True,
    )
    scdrs.preprocess(adata, cov=cov, n_mean_bin=20, n_var_bin=20)

    param = adata.uns.get("SCDRS_PARAM", {})
    if "COV_MAT" not in param:
        raise SystemExit(
            "scdrs.preprocess stored no COV_MAT; the covariate correction would "
            "be silently inactive. Check that cov was passed correctly."
        )
    print(f"  stored correction for {param['COV_BETA'].shape[0]} genes", flush=True)


def apply_cov_correction(adata, expr, gene_list):
    """Return covariate-corrected expression for a subset of genes.

        CORRECTED_X = X + COV_MAT @ COV_BETA + COV_GENE_MEAN

    Orientation follows scdrs.method (COV_BETA is indexed by gene, so it is
    transposed to (n_cov, n_gene)). `expr` is the dense (n_cell, n_gene) block
    already sliced for `gene_list`, in the same gene order.
    """
    param = adata.uns["SCDRS_PARAM"]
    cov_list = list(param["COV_MAT"])
    cov_mat = param["COV_MAT"].loc[adata.obs_names, cov_list].values
    cov_beta = param["COV_BETA"].loc[gene_list, cov_list].values.T
    gene_mean = param["COV_GENE_MEAN"].loc[gene_list].values
    return expr + cov_mat @ cov_beta + gene_mean


def spearman_vec(x, Y):
    """Spearman corr between vector x (n,) and each column of matrix Y (n, g).

    Ranks once, then Pearson on ranks — vectorized across all g genes.
    """
    xr = rankdata(x)
    xr = xr - xr.mean()
    xr_norm = np.sqrt((xr ** 2).sum())

    # Rank each gene column.
    Yr = np.apply_along_axis(rankdata, 0, Y).astype(float)
    Yr -= Yr.mean(axis=0, keepdims=True)
    Yr_norm = np.sqrt((Yr ** 2).sum(axis=0))

    denom = xr_norm * Yr_norm
    with np.errstate(invalid="ignore", divide="ignore"):
        r = (xr @ Yr) / denom
    r[denom == 0] = np.nan
    return r


def main(traits=None):
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    adata = load_expression()

    scores = pd.read_csv(SCORES, sep="\t", index_col=0)
    # Align to the scored cells FIRST so the covariate regression below is fit
    # on exactly the cell set scDRS scored (as in 06_scdrs_score.py).
    common = adata.obs_names.intersection(scores.index)
    adata = adata[common].copy()
    scores = scores.loc[common]
    print(f"{len(common)} cells shared between h5ad and scores.")

    fit_cov_correction(adata)

    # var_names / gs read after preprocessing, in case it drops any genes.
    var_names = list(adata.var_names)
    gs = scdrs.util.load_gs(
        str(GS_FILE),
        src_species="human",
        dst_species="human",
        to_intersect=adata.var_names,
    )

    trait_keys = traits if traits else [t for t in gs if t in scores.columns]

    results = {g: [] for g in GROUP_COLS}
    expr_results = {g: [] for g in GROUP_COLS}

    for trait in trait_keys:
        gene_list, gene_weights = gs[trait]
        gene_list = [g for g in gene_list if g in var_names]
        if not gene_list:
            print(f"{trait}: no risk genes overlap h5ad; skipping.")
            continue
        weight_map = dict(zip(*gs[trait]))

        gidx = [var_names.index(g) for g in gene_list]
        expr = adata.X[:, gidx]
        expr = np.asarray(expr.todense()) if hasattr(expr, "todense") else np.asarray(expr)
        # Detection indicator from the *uncorrected* values: the covariate
        # correction shifts zeros off zero, so frac_expressing must be taken
        # before it is applied.
        detected = expr > 0
        # scdrs.preprocess left adata.X untouched (sparse); apply the stored
        # correction here or the covariate adjustment does nothing.
        expr = apply_cov_correction(adata, expr, gene_list)
        score_vec = scores[trait].to_numpy()

        for group_col in GROUP_COLS:
            groups = adata.obs[group_col].astype(str).to_numpy()
            for grp in pd.unique(groups):
                mask = groups == grp
                if mask.sum() < MIN_CELLS:
                    continue
                r = spearman_vec(score_vec[mask], expr[mask])
                df = pd.DataFrame(
                    {
                        "trait": trait,
                        "group": grp,
                        "gene": gene_list,
                        "spearman_r": r,
                        "magma_weight": [weight_map[g] for g in gene_list],
                        "n_cells": int(mask.sum()),
                    }
                )
                results[group_col].append(df)

                expr_results[group_col].append(
                    pd.DataFrame(
                        {
                            "trait": trait,
                            "group": grp,
                            "gene": gene_list,
                            "mean_expr": expr[mask].mean(axis=0),
                            # Fraction on 0-1, not a percentage, despite what a
                            # name like "pct" would suggest: this is the mean of
                            # a boolean.
                            "frac_expressing": detected[mask].mean(axis=0),
                            "magma_weight": [weight_map[g] for g in gene_list],
                            "n_cells": int(mask.sum()),
                        }
                    )
                )
        print(f"{trait}: done ({len(gene_list)} risk genes).")

    for group_col in GROUP_COLS:
        full = pd.concat(results[group_col], ignore_index=True)
        full = full.sort_values(
            ["trait", "group", "spearman_r"], ascending=[True, True, False]
        )
        full.to_csv(
            OUT_DIR / f"risk_gene_drivers_{group_col}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )

        top = (
            full.dropna(subset=["spearman_r"])
            .groupby(["trait", "group"], sort=False)
            .head(TOP_N)
            .reset_index(drop=True)
        )
        top.to_csv(
            OUT_DIR / f"top{TOP_N}_risk_gene_drivers_{group_col}.tsv",
            sep="\t",
            index=False,
        )
        print(f"Wrote {group_col}: {len(full)} rows, top{TOP_N} -> {len(top)} rows.")

        # Expression levels. The z-score is taken across groups within a gene,
        # so it reads as "high in this domain relative to this gene's other
        # domains" rather than "this gene is abundant" -- otherwise the map
        # would mostly recapitulate cell-type composition.
        ex = pd.concat(expr_results[group_col], ignore_index=True)
        g = ex.groupby(["trait", "gene"])["mean_expr"]
        sd = g.transform("std")
        ex["expr_z"] = (ex["mean_expr"] - g.transform("mean")) / sd
        # A gene flat across all groups has sd 0; call that zero deviation
        # rather than propagating inf/NaN into the heatmap.
        ex.loc[sd == 0, "expr_z"] = 0.0
        ex.to_csv(
            OUT_DIR / f"risk_gene_expression_{group_col}.tsv.gz",
            sep="\t",
            index=False,
            compression="gzip",
        )
        n_genes = ex.groupby("trait")["gene"].nunique()
        print(f"Wrote {group_col} expression: {len(ex)} rows.")
        print("  risk genes on panel per trait:")
        print(n_genes.to_string())


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--traits",
        nargs="*",
        default=None,
        help="Subset of trait keys (default: all in the gs file).",
    )
    args = parser.parse_args()
    main(traits=args.traits)
