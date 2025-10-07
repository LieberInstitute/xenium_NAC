# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate STAligner

import sys
NUM_CLUST = int(sys.argv[1])

import warnings
warnings.filterwarnings("ignore")
import matplotlib.pyplot as plt
import subprocess
from pathlib import Path

# import STAligner.ST_utils
# import train_STAligner
import STAligner

# the location of R (used for the mclust clustering)
import os
os.environ['R_HOME'] = "/users/jyao/.conda/envs/proust-310/lib/R"
os.environ['R_USER'] = "/users/jyao/.conda/envs/STAligner/lib/python3.9/site-packages/rpy2"
import rpy2.robjects as robjects
import rpy2.robjects.numpy2ri

import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np
import scipy.sparse as sp
import scipy.linalg
from scipy.sparse import csr_matrix

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

# mclust
adata_concat = sc.read_h5ad(git_root / "processed-data" / "05_Clustering" / "Harmony.h5ad")
STAligner.mclust_R(adata_concat, num_cluster=NUM_CLUST, used_obsm='X_pca_harmony')
adata_concat.obs['mclust'].to_csv(
    git_root / "processed-data" / "05_Clustering" / f"Harmony_clusters.csv",
    header=True,
    index=True,
    index_label="cell"
)

# Plot individual slice
sp = adata_concat.obsm.get('spatial')
if isinstance(sp, pd.DataFrame):
    adata_concat.obsm['spatial'] = sp.values

section_ids = adata_concat.obs['Sample'].cat.remove_unused_categories().cat.categories.tolist()
spot_size = 70
for section_id in section_ids:
    print(section_id)
    ad = adata_concat[adata_concat.obs['Sample'] == section_id].copy()
    fig, ax = plt.subplots(figsize=(12, 12), dpi=300)
    sc.pl.spatial(ad, title=section_id, color="mclust", frameon=False, spot_size=spot_size, legend_fontsize=20, show=False, ax=ax)
    ax.set_title(section_id, fontsize=20)
    fig.savefig(git_root / "plots" / "05_clustering" / "Harmony" / "Pre_Annotation" / f"{section_id}_Harmony_PreAnnotation.png", bbox_inches="tight")
    plt.close(fig)