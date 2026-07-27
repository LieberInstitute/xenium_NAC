# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/

import os
os.environ["PYVISTA_OFF_SCREEN"] = "true"   # no GUI needed
os.environ.setdefault("PYVISTA_EGL", "true")  # if VTK was built with EGL
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import subprocess
from pathlib import Path
import pandas as pd
import json

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

######## VisiumHD data loading ########
# Load VisiumHD gene expression data
adata_VHD = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_sfe_counts.h5ad")
adata_VHD.obsm["spatial"] = np.asarray(adata_VHD.obsm["spatial"])[:, :2]
adata_VHD.obs['sample_id']

VHD_H1_8MTH2TQ_A1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-8MTH2TQ_A1'].copy()
VHD_H1_8MTH2TQ_D1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-8MTH2TQ_D1'].copy()
VHD_H1_M3TCP9V_A1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-M3TCP9V_A1'].copy()
VHD_H1_M3TCP9V_D1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-M3TCP9V_D1'].copy()
VHD_H1_XKYDCP3_A1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-XKYDCP3_A1'].copy()
VHD_H1_XKYDCP3_D1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-XKYDCP3_D1'].copy()
VHD_H1_XNQ4F2B_A1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-XNQ4F2B_A1'].copy()
VHD_H1_XNQ4F2B_D1 = adata_VHD[adata_VHD.obs['sample_id'] == 'H1-XNQ4F2B_D1'].copy()


spaceranger_root = Path(
    "/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/01_spaceranger"
)

out_dir = git_root / "processed-data" / "HD_Full_Analysis" / "h5ad"
out_dir.mkdir(parents=True, exist_ok=True)

bin_dir = "square_008um"

vhd_samples = {
    "H1-8MTH2TQ_A1": VHD_H1_8MTH2TQ_A1,
    "H1-8MTH2TQ_D1": VHD_H1_8MTH2TQ_D1,
    "H1-M3TCP9V_A1": VHD_H1_M3TCP9V_A1,
    "H1-M3TCP9V_D1": VHD_H1_M3TCP9V_D1,
    "H1-XKYDCP3_A1": VHD_H1_XKYDCP3_A1,
    "H1-XKYDCP3_D1": VHD_H1_XKYDCP3_D1,
    "H1-XNQ4F2B_A1": VHD_H1_XNQ4F2B_A1,
    "H1-XNQ4F2B_D1": VHD_H1_XNQ4F2B_D1,
}

summary = []

for sample_id, adata in vhd_samples.items():

    sf_path = (
        spaceranger_root
        / sample_id
        / "outs"
        / "binned_outputs"
        / bin_dir
        / "spatial"
        / "scalefactors_json.json"
    )

    if not sf_path.exists():
        raise FileNotFoundError(
            f"Missing scalefactors_json.json for {sample_id}: {sf_path}"
        )

    with open(sf_path) as f:
        sf = json.load(f)

    mpp = float(sf["microns_per_pixel"])
    bin_size_um = float(sf.get("bin_size_um", np.nan))
    spot_diameter_fullres = float(sf.get("spot_diameter_fullres", np.nan))

    # Avoid double-conversion if rerunning this cell.
    # If spatial_fullres_px already exists, use it as the original pixel coordinate.
    if "spatial_fullres_px" in adata.obsm:
        spatial_px = np.asarray(adata.obsm["spatial_fullres_px"])[:, :2].astype(float)
    else:
        spatial_px = np.asarray(adata.obsm["spatial"])[:, :2].astype(float)

    adata.obsm["spatial_fullres_px"] = spatial_px.copy()

    # Convert full-resolution pixels to microns
    adata.obsm["spatial_um"] = spatial_px * mpp
    adata.obsm["spatial"] = adata.obsm["spatial_um"].copy()
    adata.obs["x_fullres_px"] = adata.obsm["spatial_fullres_px"][:, 0]
    adata.obs["y_fullres_px"] = adata.obsm["spatial_fullres_px"][:, 1]
    adata.obs["x_um"] = adata.obsm["spatial_um"][:, 0]
    adata.obs["y_um"] = adata.obsm["spatial_um"][:, 1]

    # Store metadata
    adata.uns["spatial_unit"] = "micron"
    adata.uns["spatial_source_before_conversion"] = "full-resolution image pixels"
    adata.uns["microns_per_pixel"] = mpp
    adata.uns["bin_size_um"] = bin_size_um
    adata.uns["spaceranger_scalefactors"] = sf

    spatial_min = adata.obsm["spatial_um"].min(axis=0)
    spatial_max = adata.obsm["spatial_um"].max(axis=0)
    spatial_span = spatial_max - spatial_min

    summary.append({
        "sample_id": sample_id,
        "n_obs": adata.n_obs,
        "mpp": mpp,
        "bin_size_um": bin_size_um,
        "spot_diameter_fullres": spot_diameter_fullres,
        "spot_diameter_fullres_x_mpp": spot_diameter_fullres * mpp,
        "x_um_min": spatial_min[0],
        "y_um_min": spatial_min[1],
        "x_um_max": spatial_max[0],
        "y_um_max": spatial_max[1],
        "x_um_span": spatial_span[0],
        "y_um_span": spatial_span[1],
    })

    # Save individual sample AnnData
    sample_tag = sample_id.replace("-", "_")
    out_path = out_dir / f"VHD_{sample_tag}_spaceranger_square008um_spatial_um.h5ad"

    print(f"Saving {sample_id} to:")
    print(out_path)

    adata.write_h5ad(out_path)

summary_df = pd.DataFrame(summary)

summary_path = out_dir / "VHD_spaceranger_square008um_spatial_um_conversion_summary.csv"
summary_df.to_csv(summary_path, index=False)

