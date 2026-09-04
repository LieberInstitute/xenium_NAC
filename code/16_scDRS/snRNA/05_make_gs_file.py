"""Combine MAGMA gene-level results into a z-score matrix for scdrs munge-gs.

Same logic as parent pipeline, but reads from the snRNA MAGMA output.
"""

from pathlib import Path

import pandas as pd
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
MAGMA_GENES = BASE / "processed-data/16_scDRS/snRNA/magma/genes"
MANIFEST = BASE / "processed-data/16_scDRS/trait_manifest.tsv"
OUT_ZSCORE = BASE / "processed-data/16_scDRS/snRNA/gs/magma_zscore_matrix.tsv"

manifest = pd.read_csv(MANIFEST, sep="\t")

zdict = {}
for key in manifest["trait_key"]:
    f = MAGMA_GENES / f"{key}.genes.out"
    if not f.exists():
        print(f"WARNING: missing {f}, skipping {key}")
        continue
    df = pd.read_csv(f, sep=r"\s+")
    zdict[key] = df.set_index("GENE")["ZSTAT"]
    print(f"{key}: {len(df)} genes")

zmat = pd.DataFrame(zdict)
zmat.index.name = "GENE"
zmat.to_csv(OUT_ZSCORE, sep="\t")
print(f"Wrote {zmat.shape[0]} genes x {zmat.shape[1]} traits to {OUT_ZSCORE}")

session_info.show()
