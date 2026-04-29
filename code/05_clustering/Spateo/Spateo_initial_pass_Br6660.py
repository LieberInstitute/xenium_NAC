# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/

import os
os.environ["PYVISTA_OFF_SCREEN"] = "true"   # no GUI needed
os.environ.setdefault("PYVISTA_EGL", "true")  # if VTK was built with EGL
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

import pyvista as pv
# Fallback to a virtual framebuffer if EGL isn't available
try:
    pv.start_xvfb()
except Exception:
    pass

# Quiet deprecation + numba noise, and fix Spateo_version_
import warnings
from numba.core.errors import NumbaWarning
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")
warnings.filterwarnings("ignore", category=NumbaWarning)

import importlib.metadata as md
import spateo as st
try:
    st.__version__ = md.version("spateo-release")
except md.PackageNotFoundError:
    pass

print("Spateo version:", st.__version__)
# ---------------------------------------------------------------

import torch
device = 'cuda' if torch.cuda.is_available() else 'cpu'
print("Running this notebook on: ", device)

import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import subprocess
from pathlib import Path
import pandas as pd

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

# Load gene expression data
adata_all = ad.read_h5ad(git_root / "processed-data" / "02_build_spe" / "h5ad" / "spe_NormCounts_nucleus_normcounts.h5ad")
# Load cell type annotations
CellType_df = pd.read_csv(git_root / "processed-data" / "05_Clustering" / "Banksy_cell_types.csv")
# Add cell type info to adata_all
ct_map = CellType_df.set_index("cell_id")["CellType"]
adata_all.obs["CellType"] = adata_all.obs["cell_id"].map(ct_map)

# Subset to Br6660 sample
DONOR = "Br6660"
adata_Br6660 = adata_all[adata_all.obs["Donor"] == DONOR].copy()
adata_Br6660.obs = adata_Br6660.obs.reset_index(drop=True)
adata_Br6660.obsm['spatial'] = np.asarray(getattr(adata_Br6660.obsm['spatial'], "values", adata_Br6660.obsm['spatial']))

sample_str = adata_Br6660.obs["Sample"].astype(str)
adata_Br6660.obs["z_height"] = (
    sample_str.str.extract(r'_(\d+)$', expand=False)
              .pipe(pd.to_numeric, errors="coerce")
)

s = adata_Br6660.obs["Sample"]
labels = s.cat.categories if pd.api.types.is_categorical_dtype(s) else s.dropna().unique()

for i, v in enumerate(labels, 1):
    print(f"{i:2d}. {v}")

#  1. Br6660_NAc1_580
#  2. Br6660_NAc2_1090
#  3. Br6660_NAc3_1580
#  4. Br6660_NAc4_2080
#  5. Br6660_NAc5_2580
#  6. Br6660_NAc6_3080
#  7. Br6660_NAc7_3580
#  8. Br6660_NAc8_4580
#  9. Br6660_NAc9_5080
# 10. Br6660_Nac10_4080
# 11. Br6660_Nac11_5580


Br6660_NAc1_580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc1_580'].copy()
Br6660_NAc2_1090 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc2_1090'].copy()
Br6660_NAc3_1580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc3_1580'].copy()
Br6660_NAc4_2080 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc4_2080'].copy()
Br6660_NAc5_2580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc5_2580'].copy()
Br6660_NAc6_3080 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc6_3080'].copy()
Br6660_NAc7_3580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc7_3580'].copy()
Br6660_NAc8_4580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc8_4580'].copy()
Br6660_NAc9_5080 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_NAc9_5080'].copy()
Br6660_Nac10_4080 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_Nac10_4080'].copy()
Br6660_Nac11_5580 = adata_Br6660[adata_Br6660.obs['Sample'] == 'Br6660_Nac11_5580'].copy()

slices = [Br6660_NAc1_580, Br6660_NAc2_1090, Br6660_NAc3_1580, Br6660_NAc4_2080,
          Br6660_NAc5_2580, Br6660_NAc6_3080, Br6660_NAc7_3580, Br6660_Nac10_4080, 
          Br6660_NAc8_4580, Br6660_NAc9_5080, Br6660_Nac11_5580]


def flip_y_inplace(adata):
    Y = adata.obsm['spatial'][:, 1]
    adata.obsm['spatial'][:, 1] = (Y.max() + Y.min()) - Y  # reflect around mid-Y

for adata in slices:
    flip_y_inplace(adata)

celltype_palette = {
    "Excitatory": "#DB1C5F",
    "Microglia_A": "#0D87E4",
    "D1_Island_B": "#00FF0D",
    "Astro_A": "#FD00FD",
    "Fibroblast_A": "#FFA7E2",
    "Fibroblast_B": "#2AFECA",
    "DRD2_MSN": "#7AAA16",
    "DRD1_MSN": "#9400FF",
    "Astro_B": "#823526",
    "MSN_Oligo": "#F5DEC0",
    "Inh_PVALB": "#93E5FF",
    "OPC": "#7A1699",
    "Microglia_B": "#FF0DBC",
    "Inh_SST": "#C4B3FB",
    "Astrocyte_Oligo": "#F80D2A",
    "D1_Island_A": "#224B82",
    "CHAT": "#FBA475",
    "Ependymal": "#73690D",
    "Microglia_Oligo": "#B62A7A",
    "WM": "orange",
}

# Plotting before Spateo
import matplotlib
matplotlib.use("Agg")

plt.ioff()

plot_outdir = git_root / "plots" / "05_clustering" / "Spateo"
plot_outdir.mkdir(parents=True, exist_ok=True)

st.pl.slices_2d(
    slices=slices,
    spatial_key="spatial",
    label_key="CellType",
    palette=celltype_palette,
    show_legend=True,
    point_size=0.04,
    legend_kwargs={
        "loc": "upper center",
        "bbox_to_anchor": (0.5, 0),
        "ncol": 5,
        "borderaxespad": -4,
        "frameon": False,
    },
)

plt.savefig(
    plot_outdir / "Br6660_slices_2d_before_spateo.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close("all")


# Running Spateo on the Br6660 sample
spatial_key = 'spatial'
key_added = 'align_spatial'
transformation = st.align.morpho_align_transformation(
    models=slices,
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    # sparse and chunk calculation
    sparse_calculation_mode=True,
    use_chunk=True,
    chunk_capacity=2,
    verbose=True,
    partial_robust_level=30
)

# Save the transformation matrix
import pickle
transformation_path = git_root / "processed-data" / "05_Clustering" / "Spateo" / "Spateo_transformation_Br6660.pkl"
with open(transformation_path, "wb") as f:
    pickle.dump(transformation, f, protocol=pickle.HIGHEST_PROTOCOL)
## read transformation
# with open(transformation_path, "rb") as f:
#     transformation = pickle.load(f)

# Apply the transformation to the slices
aligned_slices = st.align.morpho_align_apply_transformation(
    models=slices,
    spatial_key=spatial_key,
    key_added=key_added,
    transformation=transformation,
)

# Overlaid plot of aligned slices (first pass)
plt.ioff()
st.pl.overlay_slices_2d(
    slices=aligned_slices,
    spatial_key=key_added,
    height=2,
    overlay_type="backward",
)

plt.savefig(
    plot_outdir / "Br6660_overlay_slices_2d_initial_spateo.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close("all")

# Overlaid plot of aligned slices colored by cell type (first pass)
import re

plt.ioff()

plot_outdir = git_root / "plots" / "05_clustering" / "Spateo" / "Br6660_initial_pass_by_celltype"
plot_outdir.mkdir(parents=True, exist_ok=True)

# get cell types, then sort alphabetically (case-insensitive)
celltypes = (
    adata_all.obs["CellType"].unique().tolist()
    if "adata_all" in globals()
    else sorted(set().union(*[ad.obs["CellType"].unique().tolist() for ad in aligned_slices]))
)
celltypes = sorted(celltypes, key=lambda x: str(x).lower())

for ct in celltypes:
    aligned_slices_ct = []
    for ad in aligned_slices:
        ad_ct = ad[ad.obs["CellType"] == ct].copy()
        if ad_ct.n_obs > 0:
            aligned_slices_ct.append(ad_ct)
    if not aligned_slices_ct:
        continue
    ct_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(ct))
    st.pl.overlay_slices_2d(
        slices=aligned_slices_ct,
        spatial_key=key_added,
        height=2,
        overlay_type="backward",
    )
    plt.savefig(
        plot_outdir / f"{ct_safe}_overlay_slices_2d.png",
        dpi=300,
        bbox_inches="tight",
    )
    plt.close("all")

