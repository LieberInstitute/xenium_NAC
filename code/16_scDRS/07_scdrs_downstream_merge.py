"""Merge per-trait downstream results and export cell metadata.

Run after all array tasks from 07_scdrs_downstream_one.py complete.
"""

from pathlib import Path

import pandas as pd
import scanpy as sc
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/HD_Full_Analysis/h5ad/VHD_sfe_counts.h5ad"
MANIFEST = BASE / "processed-data/16_scDRS/trait_manifest.tsv"
OUT_DIR = BASE / "processed-data/16_scDRS/group_analysis"
GROUP_COLS = ["Spatial_Domain", "labels"]

EXPECTED_TRAITS = set(pd.read_csv(MANIFEST, sep="\t")["trait_key"])


def check_traits(found: set, label: str) -> None:
    """Abort if the globbed files do not match trait_manifest.tsv exactly.

    Guards two silent failure modes. A *missing* trait means an array task
    failed and the merge would quietly emit an incomplete table. An
    *unexpected* trait means a stale file from an earlier run (with a
    different manifest) would be folded into the merged output.
    """
    missing = EXPECTED_TRAITS - found
    extra = found - EXPECTED_TRAITS
    if missing or extra:
        raise SystemExit(
            f"{label}: file set does not match trait_manifest.tsv "
            f"({len(found)} found, {len(EXPECTED_TRAITS)} expected)\n"
            f"  missing (failed task?): {sorted(missing) if missing else 'none'}\n"
            f"  unexpected (stale?)   : {sorted(extra) if extra else 'none'}"
        )

# --- Cell metadata ---
print("Exporting cell metadata", flush=True)
adata = sc.read_h5ad(H5AD, backed="r")

# Remove cells outside the NAc (Hypo and Excitatory spatial domains)
# and Excitatory cell-type labels from any remaining spatial domains
sd_exclude = ["Hypo", "Excitatory"]
label_exclude = ["Excitatory"]
keep = ~adata.obs["Spatial_Domain"].isin(sd_exclude) & ~adata.obs["labels"].isin(label_exclude)
print(f"Excluding {(~keep).sum()} cells (Spatial_Domain in {sd_exclude} or labels in {label_exclude})", flush=True)

meta_cols = ["Sample", "X_coords", "Y_coords"] + GROUP_COLS
adata.obs.loc[keep, meta_cols].to_csv(
    OUT_DIR / "cell_metadata.tsv.gz", sep="\t", compression="gzip"
)
del adata

# --- Merge group stats ---
for g in GROUP_COLS:
    prefix = f"scdrs_group_stats_{g}_"
    files = sorted(OUT_DIR.glob(f"{prefix}*.tsv"))
    check_traits(
        {f.name.removeprefix(prefix).removesuffix(".tsv") for f in files},
        f"group stats [{g}]",
    )
    dfs = [pd.read_csv(f, sep="\t") for f in files]
    merged = pd.concat(dfs, ignore_index=True)
    out_path = OUT_DIR / f"scdrs_group_stats_{g}.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged {len(files)} files -> {out_path} ({merged.shape})", flush=True)

# --- Merge norm scores ---
norm_files = sorted(OUT_DIR.glob("norm_score_*.tsv.gz"))
check_traits(
    {f.name.removeprefix("norm_score_").removesuffix(".tsv.gz") for f in norm_files},
    "norm scores",
)
norm_dict = {}
for f in norm_files:
    trait = f.name.removeprefix("norm_score_").removesuffix(".tsv.gz")
    df = pd.read_csv(f, sep="\t", index_col=0)
    norm_dict[trait] = df["norm_score"]
norm_all = pd.DataFrame(norm_dict)
norm_all.to_csv(
    OUT_DIR / "norm_scores_all_traits.tsv.gz", sep="\t", compression="gzip"
)
print(f"Merged norm scores: {norm_all.shape}", flush=True)

print("Done.", flush=True)
session_info.show()
