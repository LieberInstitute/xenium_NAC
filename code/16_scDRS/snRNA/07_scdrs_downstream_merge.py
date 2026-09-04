"""Merge per-trait downstream results for snRNA-seq scDRS.

Run after all array tasks from 07_scdrs_downstream_one.py complete.
"""

from pathlib import Path

import pandas as pd
import scanpy as sc
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/16_scDRS/snRNA/h5ad/snRNA_counts.h5ad"
OUT_DIR = BASE / "processed-data/16_scDRS/snRNA/group_analysis"
GROUP_COLS = ["CellType.Final"]

# --- Cell metadata ---
print("Exporting cell metadata", flush=True)
adata = sc.read_h5ad(H5AD, backed="r")

# Remove Excitatory and Neuron_Ambig cell types
ct_exclude = ["Excitatory", "Neuron_Ambig"]
keep = ~adata.obs["CellType.Final"].isin(ct_exclude)
print(f"Excluding {(~keep).sum()} cells (CellType.Final in {ct_exclude})", flush=True)

meta_cols = ["Sample"] + GROUP_COLS
adata.obs.loc[keep, meta_cols].to_csv(
    OUT_DIR / "cell_metadata.tsv.gz", sep="\t", compression="gzip"
)
del adata

# --- Merge group stats ---
for g in GROUP_COLS:
    files = sorted(OUT_DIR.glob(f"scdrs_group_stats_{g}_*.tsv"))
    dfs = [pd.read_csv(f, sep="\t") for f in files]
    merged = pd.concat(dfs, ignore_index=True)
    out_path = OUT_DIR / f"scdrs_group_stats_{g}.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged {len(files)} files -> {out_path} ({merged.shape})", flush=True)

# --- Merge norm scores ---
norm_files = sorted(OUT_DIR.glob("norm_score_*.tsv.gz"))
norm_dict = {}
for f in norm_files:
    trait = f.name.replace("norm_score_", "").replace(".tsv.gz", "")
    df = pd.read_csv(f, sep="\t", index_col=0)
    norm_dict[trait] = df["norm_score"]
norm_all = pd.DataFrame(norm_dict)
norm_all.to_csv(
    OUT_DIR / "norm_scores_all_traits.tsv.gz", sep="\t", compression="gzip"
)
print(f"Merged norm scores: {norm_all.shape}", flush=True)

print("Done.", flush=True)
session_info.show()
