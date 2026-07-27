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

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scanpy as sc
import anndata as ad
import numpy as np
import subprocess
from pathlib import Path
import pandas as pd
import scipy.sparse as sp

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

## Load VHD data
# VHD_H1_8MTH2TQ_A1
# VHD_H1_8MTH2TQ_D1
# VHD_H1_M3TCP9V_A1
# VHD_H1_M3TCP9V_D1
# VHD_H1_XKYDCP3_A1
# VHD_H1_XKYDCP3_D1
# VHD_H1_XNQ4F2B_A1
# VHD_H1_XNQ4F2B_D1

VHD_H1_8MTH2TQ_A1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_8MTH2TQ_A1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_8MTH2TQ_D1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_8MTH2TQ_D1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_M3TCP9V_A1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_M3TCP9V_A1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_M3TCP9V_D1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_M3TCP9V_D1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_XKYDCP3_A1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_XKYDCP3_A1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_XKYDCP3_D1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_XKYDCP3_D1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_XNQ4F2B_A1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_XNQ4F2B_A1_spaceranger_square008um_spatial_um.h5ad")
VHD_H1_XNQ4F2B_D1 = ad.read_h5ad(git_root / "processed-data" / "HD_Full_Analysis" / "h5ad" / "VHD_H1_XNQ4F2B_D1_spaceranger_square008um_spatial_um.h5ad")

slices = [
    VHD_H1_8MTH2TQ_A1,
    VHD_H1_8MTH2TQ_D1,
    VHD_H1_M3TCP9V_A1,
    VHD_H1_M3TCP9V_D1,
    VHD_H1_XKYDCP3_A1,
    VHD_H1_XKYDCP3_D1,
    VHD_H1_XNQ4F2B_A1,
    VHD_H1_XNQ4F2B_D1,
]

slice_names = [
    "VHD_H1_8MTH2TQ_A1",
    "VHD_H1_8MTH2TQ_D1",
    "VHD_H1_M3TCP9V_A1",
    "VHD_H1_M3TCP9V_D1",
    "VHD_H1_XKYDCP3_A1",
    "VHD_H1_XKYDCP3_D1",
    "VHD_H1_XNQ4F2B_A1",
    "VHD_H1_XNQ4F2B_D1",
]

for slice_adata in slices:
    xy = slice_adata.obsm["spatial"]  # should already be micron coordinates
    x = xy[:, 0]
    y = xy[:, 1]
    y_min, y_max = y.min(), y.max()
    xy_flip_y = xy.copy()
    xy_flip_y[:, 1] = y_min + y_max - y
    slice_adata.obsm["spatial"] = xy_flip_y


## Plot before aligning individuals to Xenium sections
plt.ioff()

plot_outdir = git_root / "plots" / "11_Xenium_VisiumHD_alignment" / "VHD_slices"
plot_outdir.mkdir(parents=True, exist_ok=True)

for slice_name, slice_adata in zip(slice_names, slices):
    for label_key in ["snRNA_label", "Spatial_Domain"]:
        st.pl.slices_2d(
            slices=[slice_adata],
            spatial_key="spatial",
            label_key=label_key,
            show_legend=True,
            point_size=0.04,
            legend_kwargs={
                "loc": "upper center",
                "bbox_to_anchor": (0.5, -0.08),
                "ncol": 5,
                "borderaxespad": 0,
                "frameon": False,
                "fontsize": 5,
                "markerscale": 0.5,
                "handletextpad": 0.2,
                "columnspacing": 0.6,
            },
        )
        fig = plt.gcf()
        for ax in fig.axes:
            ax.set_title(slice_name)
        fig.subplots_adjust(bottom=0.22)

        plt.savefig(
            plot_outdir / f"{slice_name}_{label_key}.png",
            dpi=300,
            bbox_inches="tight",
        )
        plt.close("all")


## Pairwise VHD-to-Xenium overlays before VHD alignment
# Xenium is shown in its final donor-level aligned coordinate space because it
# is the fixed reference used by Modules 04/05. VHD uses the unaligned,
# micron-scale, y-flipped coordinates prepared above.
REFERENCE_COLOR = "#4C78A8"
MOVING_COLOR = "#F58518"
PAIRWISE_ALL_POINT_SIZE = 0.20
PAIRWISE_ALL_ALPHA = 0.45
PAIRWISE_LABEL_POINT_SIZE = 0.60
PAIRWISE_XENIUM_ALPHA = 0.60
PAIRWISE_VHD_ALPHA = 0.35

pair_specs = [
    # donor, Xenium sample, VHD sample, Xenium depth, VHD depth
    ("Br6660", "Br6660_NAc4_2080", "VHD_H1_8MTH2TQ_A1", 2080, 2040),
    ("Br6660", "Br6660_NAc6_3080", "VHD_H1_M3TCP9V_A1", 3080, 3040),
    ("Br6660", "Br6660_NAc7_3580", "VHD_H1_XKYDCP3_A1", 3580, 3540),
    ("Br6660", "Br6660_NAc7_3580", "VHD_H1_8MTH2TQ_D1", 3580, 3630),
    ("Br6660", "Br6660_Nac10_4080", "VHD_H1_XKYDCP3_D1", 4080, 4120),
    ("Br6660", "Br6660_NAc8_4580", "VHD_H1_M3TCP9V_D1", 4580, 4550),
    ("Br6436", "Br6436_Nac1_650", "VHD_H1_XNQ4F2B_A1", 650, 630),
    ("Br6436", "Br6436_Nac_9_4650", "VHD_H1_XNQ4F2B_D1", 4650, 4630),
]
vhd_by_name = dict(zip(slice_names, slices))

celltype_table = pd.read_csv(
    git_root / "processed-data" / "05_Clustering" / "Banksy_cell_types.csv",
    usecols=["cell_id", "CellType"],
    dtype={"cell_id": "string"},
)
if (
    celltype_table["cell_id"].isna().any()
    or celltype_table["cell_id"].duplicated().any()
):
    raise ValueError("Banksy cell-type table contains invalid cell IDs")
celltype_map = celltype_table.set_index("cell_id")["CellType"]


def load_aligned_xenium_coordinates(donor):
    coordinate_path = (
        git_root
        / "processed-data"
        / "05_Xenium_alignment"
        / f"module02_Spateo_second_pass_{donor}"
        / f"{donor}_second_pass_coordinates.csv"
    )
    coordinates = pd.read_csv(
        coordinate_path,
        usecols=[
            "donor",
            "sample",
            "cell_id",
            "x_second_pass",
            "y_second_pass",
            "z_height",
        ],
        dtype={"donor": "string", "sample": "string", "cell_id": "string"},
    )
    coordinates = coordinates.loc[coordinates["donor"].eq(donor)].copy()
    if coordinates["cell_id"].isna().any():
        raise ValueError(f"{donor}: missing Xenium cell IDs")
    if coordinates["cell_id"].duplicated().any():
        raise ValueError(f"{donor}: duplicate Xenium cell IDs")
    coordinates["CellType"] = coordinates["cell_id"].map(celltype_map)
    return coordinates


xenium_by_donor = {
    donor: load_aligned_xenium_coordinates(donor)
    for donor in ("Br6660", "Br6436")
}


def short_xenium_name(sample, donor):
    prefix = f"{donor}_"
    return sample[len(prefix):] if sample.startswith(prefix) else sample


def plot_unaligned_pairwise_grid(cell_type, output_path):
    figure, axes = plt.subplots(
        2,
        4,
        figsize=(20.0, 9.0),
        squeeze=False,
    )
    for axis, (
        donor,
        xenium_sample,
        vhd_sample,
        xenium_depth,
        vhd_depth,
    ) in zip(axes.ravel(), pair_specs):
        xenium = xenium_by_donor[donor]
        xenium = xenium.loc[xenium["sample"].eq(xenium_sample)]
        if xenium.empty:
            raise ValueError(f"Missing Xenium sample: {xenium_sample}")
        observed_depths = xenium["z_height"].dropna().unique()
        if len(observed_depths) != 1 or not np.isclose(
            float(observed_depths[0]), xenium_depth
        ):
            raise ValueError(
                f"Unexpected Xenium depth for {xenium_sample}: "
                f"{observed_depths.tolist()}"
            )

        vhd = vhd_by_name[vhd_sample]
        vhd_xy = np.asarray(vhd.obsm["spatial"], dtype=float)[:, :2]
        vhd_labels = vhd.obs["Spatial_Domain"].astype("string")

        if cell_type is not None:
            xenium = xenium.loc[xenium["CellType"].eq(cell_type)]
            vhd_mask = (
                vhd_labels.eq(cell_type)
                .fillna(False)
                .to_numpy(dtype=bool)
            )
            vhd_xy = vhd_xy[vhd_mask]
            point_size = PAIRWISE_LABEL_POINT_SIZE
            xenium_alpha = PAIRWISE_XENIUM_ALPHA
            vhd_alpha = PAIRWISE_VHD_ALPHA
        else:
            point_size = PAIRWISE_ALL_POINT_SIZE
            xenium_alpha = PAIRWISE_ALL_ALPHA
            vhd_alpha = PAIRWISE_ALL_ALPHA

        axis.scatter(
            vhd_xy[:, 0],
            vhd_xy[:, 1],
            s=point_size,
            alpha=vhd_alpha,
            color=MOVING_COLOR,
            linewidths=0,
            label="VHD",
            zorder=1,
            rasterized=True,
        )
        axis.scatter(
            xenium["x_second_pass"],
            xenium["y_second_pass"],
            s=point_size,
            alpha=xenium_alpha,
            color=REFERENCE_COLOR,
            linewidths=0,
            label="Xenium",
            zorder=2,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (µm)")
        axis.set_ylabel("y (µm)")
        axis.set_title(
            f"{donor}: {short_xenium_name(xenium_sample, donor)} ↔ "
            f"{vhd_sample.removeprefix('VHD_')}\n"
            f"Xenium depth: {xenium_depth:g} µm | "
            f"VHD depth: {vhd_depth:g} µm",
            fontsize=8,
        )
        axis.legend(
            loc="upper right",
            fontsize=6,
            markerscale=7,
            frameon=False,
        )

    label = "all observations" if cell_type is None else cell_type
    figure.suptitle(
        f"Unaligned Xenium–VisiumHD pairwise overlays: {label}",
        y=1.002,
        fontsize=13,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


pairwise_outdir = plot_outdir / "pairwise_2d_unaligned"
pairwise_outdir.mkdir(parents=True, exist_ok=True)
plot_unaligned_pairwise_grid(
    None,
    pairwise_outdir
    / "Xenium_VHD_pairwise_unaligned_all_observations.png",
)
for label in ("D1_Island_A", "D1_Island_B", "WM"):
    plot_unaligned_pairwise_grid(
        label,
        pairwise_outdir / f"Xenium_VHD_pairwise_unaligned_{label}.png",
    )

## Pairing: VHD and Xenium sections
# Br6660
# Run 1
# VHD_H1_XKYDCP3_A1 (3540um); Br6660_NAc7_3580 (3580um)
# VHD_H1_XKYDCP3_D1 (4120um); Br6660_Nac10_4080 (4080um)

# Run 2
# VHD_H1_M3TCP9V_A1 (3040um); Br6660_NAc6_3080 (3080um)
# VHD_H1_M3TCP9V_D1 (4550um); Br6660_NAc8_4580 (4580um)

# Run 4
# VHD_H1_8MTH2TQ_A1 (2040um); Br6660_NAc4_2080 (2080um)
# VHD_H1_8MTH2TQ_D1 (3630um); Br6660_NAc7_3580 (3580um)

# Br6436
# Run 3
# VHD_H1_XNQ4F2B_A1 (630um); Br6436_Nac1_650 (650um)
# VHD_H1_XNQ4F2B_D1 (4630um); Br6436_Nac_9_4650 (4650um)
