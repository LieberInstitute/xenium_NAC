# srun --partition=gpu --mem=500G --gres=gpu:1 --pty /bin/bash
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate STAligner (where mclust and harmony are installed)

import subprocess
from pathlib import Path

import scanpy as sc
import anndata as ad
import harmonypy as hm
import matplotlib
matplotlib.use("Agg")             
import matplotlib.pyplot as plt


git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

adata_all = ad.read_h5ad(git_root / "processed-data" / "02_build_spe" / "h5ad" / "spe_NormCounts_nucleus_normcounts.h5ad")
gene_mask = adata_all.var["Type"] == "Gene Expression"
adata_all = adata_all[:, gene_mask].copy()

DONOR = "Br6660"
if DONOR != "All":
    adata_all = adata_all[adata_all.obs["Donor"] == DONOR].copy()

# Calculate PCA
# adata_all.layers["counts"] = adata_all.X.copy()
sc.pp.pca(adata_all, n_comps=30)

######### Run integration with harmony
# prepare metadata and PCA
meta_data = adata_all.obs
data_mat = adata_all.obsm['X_pca']
ho = hm.run_harmony(data_mat, meta_data, ['Sample'])
# mapping back the result to the adata object
adata_all.obsm['X_pca_harmony'] = ho.Z_corr.T

# Save the integrated data for
adata_all.write(git_root / "processed-data" / "05_Clustering" / f"Harmony_{DONOR}.h5ad")

sc.pp.neighbors(adata_all, use_rep='X_pca_harmony')
sc.tl.umap(adata_all)

plt.rcParams['font.sans-serif'] = "Arial"
plt.rcParams["figure.figsize"] = (3, 3)
plt.rcParams['font.size'] = 10
ax = sc.pl.umap(adata_all, color="Sample", show=False)
plt.savefig(git_root / "plots" / "05_clustering" / "Harmony" / "Harmony_umap.png", dpi=300, bbox_inches="tight")
plt.close()
