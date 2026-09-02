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
from matplotlib.patches import Patch
import scanpy as sc
import anndata as ad
import numpy as np
import subprocess
from pathlib import Path
import pandas as pd

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

MODULE_NAME = "cell_type_grid_analysis"
OUTPUT_DIR = git_root / "processed-data" / "09_spatial_gradient" / MODULE_NAME
PLOT_DIR = git_root / "plots" / "09_spatial_gradient" / MODULE_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

# Load gene expression data
adata_all = ad.read_h5ad(git_root / "processed-data" / "02_build_spe" / "h5ad" / "spe_NormCounts_nucleus_normcounts.h5ad")
# Load cell type annotations
CellType_df = pd.read_csv(git_root / "processed-data" / "05_Clustering" / "Banksy_cell_types.csv")
# Add cell type info to adata_all
ct_map = CellType_df.set_index("cell_id")["CellType"]
adata_all.obs["CellType"] = adata_all.obs["cell_id"].map(ct_map)

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

# Subset to Br6660 sample
DONOR = "Br6660"
adata_Br6660 = adata_all[adata_all.obs["Donor"] == DONOR].copy()
adata_Br6660.obs = adata_Br6660.obs.reset_index(drop=True)
adata_Br6660.obsm['spatial'] = np.asarray(getattr(adata_Br6660.obsm['spatial'], "values", adata_Br6660.obsm['spatial']))

# Load the final (second-pass) Spateo coordinates for Br6660
csv_path = (
    git_root
    / "processed-data"
    / "05_Xenium_alignment"
    / "module02_Spateo_second_pass_Br6660"
    / "Br6660_second_pass_coordinates.csv"
)
aligned_coord_df = pd.read_csv(csv_path)
required_coordinate_columns = {
    "donor", "sample", "cell_id", "x_second_pass", "y_second_pass", "z_height"
}
missing_coordinate_columns = required_coordinate_columns - set(aligned_coord_df.columns)
if missing_coordinate_columns:
    raise ValueError(
        f"Second-pass coordinate CSV is missing columns: "
        f"{sorted(missing_coordinate_columns)}"
    )
if set(aligned_coord_df["donor"].dropna().astype(str)) != {DONOR}:
    raise ValueError("Second-pass coordinate CSV contains an unexpected donor")

adata_Br6660.obs["cell_id"] = adata_Br6660.obs["cell_id"].astype(str)
aligned_coord_df["cell_id"] = aligned_coord_df["cell_id"].astype(str)
if aligned_coord_df["cell_id"].duplicated().any():
    raise ValueError("Second-pass coordinate CSV contains duplicate cell IDs")
coord_map = aligned_coord_df.set_index("cell_id")

missing = set(adata_Br6660.obs["cell_id"]) - set(coord_map.index)
if missing:
    raise ValueError(
        f"{len(missing)} Br6660 cells were not found in the second-pass "
        "coordinate CSV"
    )
# Reorder aligned coordinates to match adata_Br6660.obs cell order
matched = coord_map.loc[adata_Br6660.obs["cell_id"]]
adata_Br6660.obsm["align_spatial"] = matched[
    ["x_second_pass", "y_second_pass"]
].to_numpy()
adata_Br6660.obs["z_height"] = matched["z_height"].to_numpy()

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

aligned_slices = [Br6660_NAc1_580, Br6660_NAc2_1090, Br6660_NAc3_1580, Br6660_NAc4_2080,
          Br6660_NAc5_2580, Br6660_NAc6_3080, Br6660_NAc7_3580, Br6660_Nac10_4080, 
          Br6660_NAc8_4580, Br6660_NAc9_5080, Br6660_Nac11_5580]
# Concatenated anndata with correct slice order and aligned coordinates
aligned_adata = ad.concat(aligned_slices)
aligned_adata.obsm['spatial_3D'] = np.concatenate([aligned_adata.obsm['align_spatial'], np.array(aligned_adata.obs['z_height'].values)[:,None]], axis=1)

# Overlaid plot of aligned consecutive slices (sanity check)
key_added = "align_spatial"
plot_outdir = PLOT_DIR


def display_sample_name(adata):
    """Return the single sample name without the exact donor prefix."""
    sample_names = adata.obs["Sample"].dropna().astype(str).unique()
    if len(sample_names) != 1:
        raise ValueError("Each plotted slice must contain exactly one Sample")
    sample_name = str(sample_names[0])
    donor_prefix = f"{DONOR}_"
    return (
        sample_name[len(donor_prefix):]
        if sample_name.startswith(donor_prefix)
        else sample_name
    )


def plot_consecutive_overlays(plotted_slices, spatial_key, output_path):
    """Plot final second-pass slice pairs in the established overlay format."""
    n_pairs = len(plotted_slices) - 1
    if n_pairs < 1:
        raise ValueError("At least two slices are required for an overlay plot")
    n_columns = min(3, n_pairs)
    n_rows = int(np.ceil(n_pairs / n_columns))
    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.ravel()
    display_names = [display_sample_name(adata) for adata in plotted_slices]

    for pair_index, axis in enumerate(axes_flat[:n_pairs]):
        reference_coordinates = np.asarray(
            plotted_slices[pair_index].obsm[spatial_key]
        )[:, :2]
        moving_coordinates = np.asarray(
            plotted_slices[pair_index + 1].obsm[spatial_key]
        )[:, :2]
        reference_name = display_names[pair_index]
        moving_name = display_names[pair_index + 1]
        axis.scatter(
            reference_coordinates[:, 0],
            reference_coordinates[:, 1],
            s=0.2,
            alpha=0.45,
            color="#4C78A8",
            linewidths=0,
            label=reference_name,
        )
        axis.scatter(
            moving_coordinates[:, 0],
            moving_coordinates[:, 1],
            s=0.2,
            alpha=0.45,
            color="#F58518",
            linewidths=0,
            label=moving_name,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(
            f"Final Spateo aligned: {reference_name} → {moving_name}",
            fontsize=9,
        )
        axis.legend(
            loc="upper right", fontsize=6, markerscale=8, frameon=False
        )

    for axis in axes_flat[n_pairs:]:
        axis.axis("off")

    figure.suptitle(
        "Final second-pass Spateo consecutive-slice overlays", y=1.002
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


plt.ioff()
plot_consecutive_overlays(
    aligned_slices,
    key_added,
    plot_outdir / "Br6660_spateo_sanity_check.png",
)

coords = aligned_adata.obsm[key_added]
plt.ioff()
plt.figure(figsize=(5, 5))
plt.scatter(coords[:, 0], coords[:, 1], s=1)
plt.gca().set_aspect("equal")
plt.xlabel("x")
plt.ylabel("y")
plt.title("Overlaid XY coordinates of Br6660 slices")
plt.savefig(
    plot_outdir / "Br6660_XY.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close("all")


############## Grid plot by slice
celltype_key  = "CellType"      
sample_key    = "Sample"            # slice/depth label
grid_size     = 5                   # 5×5 grid
pt_size       = 0.5                   # scatter size
pt_alpha      = 0.6                 # point alpha
figsize       = (18, 10)            # figure size for ~11 subplots
# -------------------------------
XY = aligned_adata.obsm[key_added]
x, y = XY[:, 0], XY[:, 1]

obs = aligned_adata.obs
samples_all = obs[sample_key].astype(str)
sample_order = (
    obs[[sample_key, "z_height"]]
    .groupby(sample_key)["z_height"].mean()
    .sort_values()
    .index.astype(str)
    .tolist()
)

# Global grid edges from ALL data (so every slice gets the same bins)
xbins = np.linspace(x.min(), x.max(), grid_size + 1)
ybins = np.linspace(y.min(), y.max(), grid_size + 1)

# Precompute per-cell colors from palette
celltypes_order = list(celltype_palette.keys())
celltypes = obs[celltype_key].astype(str)
colors = celltypes.map(celltype_palette).to_numpy()

# Plotting
n = len(sample_order)
ncols = 4
nrows = int(np.ceil(n / ncols))

plt.ioff()
fig, axes = plt.subplots(nrows, ncols, figsize=figsize, squeeze=False)
axes = axes.ravel()
for idx, sample in enumerate(sample_order):
    ax = axes[idx]
    mask = (samples_all.astype(str).values == sample)
    ax.scatter(x[mask], y[mask], s=pt_size, c=colors[mask], alpha=pt_alpha)
    ax.set_xlim(x.min(), x.max())   
    ax.set_ylim(y.min(), y.max())   
    ax.set_aspect("equal")
    ax.set_title(f"Slice {idx}: {sample}", fontsize=12)
    ax.set_xticks([])
    ax.set_yticks([])
    for xv in xbins:
        ax.axvline(x=xv, color="k", lw=0.6, alpha=0.6)
    for yv in ybins:
        ax.axhline(y=yv, color="k", lw=0.6, alpha=0.6)

for k in range(idx + 1, nrows * ncols):
    axes[k].axis("off")

# Legend
legend_handles = [Patch(facecolor=celltype_palette[ct], label=ct) for ct in celltypes_order]
fig.legend(
    legend_handles,
    celltypes_order,
    loc="center left",
    bbox_to_anchor=(0.88, 0.5),
    ncol=1,
    frameon=False,
    fontsize=9,
)

plt.subplots_adjust(right=0.86, wspace=0.15, hspace=0.25)
plt.savefig(
    plot_outdir / "Br6660_slice_grid_overlay.png",
    dpi=300,
    bbox_inches="tight",
)
plt.close("all")


############## Grid analysis and proportion plots
spatial_key   = "spatial_3D"        # (N,3): columns [x,y,z]
n_bars        = 10                  # bars per subplot for XZ/YZ (bins along remaining axis)
figsize_xy    = (16, 14)
figsize_xz_yz = (16, 14)
bars_over     = "auto"              # "auto" (recommended) or "z"
annot_thresh  = 0.05                # show prop text if > this

coords = aligned_adata.obsm[spatial_key]
x = coords[:, 0]
y = coords[:, 1]
z = coords[:, 2]

# Sample order + z_height lookup for XY bars (as before)
z_by_sample = obs[[sample_key, "z_height"]].groupby(sample_key)["z_height"].mean()
samples_xy  = z_by_sample.sort_values().index.astype(str).tolist()
z_lookup    = z_by_sample.to_dict()    # {sample: mean z}


# Helper to build stacked-bar grid for a given plane
def plot_grid_props_plane(plane="xy"):
    """
    plane: "xy", "xz", or "yz"
    - xy: grid over (x,y); bars per sample (ordered by z_height)
    - xz: grid over (x,z); bars = n_bars bins over remaining axis (y) [or z if bars_over="z"]
    - yz: grid over (y,z); bars = n_bars bins over remaining axis (x) [or z if bars_over="z"]
    """
    if plane == "xy":
        X, Y = x, y
        bar_ids = samples_all                # per-slice
        bar_labels = samples_xy
        use_samples = True
        suptitle = "Cell-type Proportions (XY grid; bars = slices; x-label = z_height)"
        figsize = figsize_xy
    elif plane == "xz":
        X, Y = x, z
        bar_axis = z if bars_over == "z" else y
        use_samples = False
        suptitle = f"Cell-type Proportions (XZ grid; bars = {n_bars} bins over {'z' if bars_over=='z' else 'y'})"
        figsize = figsize_xz_yz
    else:  # "yz"
        X, Y = y, z
        bar_axis = z if bars_over == "z" else x
        use_samples = False
        suptitle = f"Cell-type Proportions (YZ grid; bars = {n_bars} bins over {'z' if bars_over=='z' else 'x'})"
        figsize = figsize_xz_yz
    df = pd.DataFrame({
        "X": X, "Y": Y,
        "x_all": x, "y_all": y, "z_all": z,      # keep original axes for binning
        "Sample": samples_all,
        "CellType": celltypes
    })

    # Shared grid edges on the chosen plane
    Xbins = np.linspace(df["X"].min(), df["X"].max(), grid_size + 1)
    Ybins = np.linspace(df["Y"].min(), df["Y"].max(), grid_size + 1)
    df["Xbin"] = np.clip(np.digitize(df["X"], Xbins, right=False) - 1, 0, grid_size - 1)
    df["Ybin"] = np.clip(np.digitize(df["Y"], Ybins, right=False) - 1, 0, grid_size - 1)

    # Create figure & axes 
    plt.ioff()
    fig, axes = plt.subplots(
        grid_size, grid_size, figsize=figsize,
        constrained_layout=False
    )
    plt.subplots_adjust(left=0.06, right=0.86, top=0.92, bottom=0.08, hspace=0.35, wspace=0.20)
    fig.suptitle(suptitle, y=0.98, fontsize=14)

    # Legend
    legend_handles = [Patch(facecolor=celltype_palette[ct], label=ct) for ct in celltypes_order]
    fig.legend(legend_handles, celltypes_order, loc="center left", bbox_to_anchor=(0.88, 0.5),
               ncol=1, frameon=False, fontsize=8)

    # For XZ/YZ: setup bar bins & labels
    if not use_samples:
        bar_edges = np.linspace(bar_axis.min(), bar_axis.max(), n_bars + 1)
        bar_centers = 0.5 * (bar_edges[:-1] + bar_edges[1:])
        bar_labels = [f"{c:.0f}" for c in bar_centers]

    # Loop over grid cells
    for i in range(grid_size):          # columns (left to right)
        for j in range(grid_size):      # rows (bottom to top)
            ax = axes[grid_size - 1 - j, i]  # lower-left origin
            sub = df[(df["Xbin"] == i) & (df["Ybin"] == j)]
            if sub.empty:
                ax.axis("off")
                continue

            if use_samples:
                # XY: counts per (Sample × CellType) then calculate proportions within each Sample
                ct_counts = sub.groupby(["Sample", "CellType"]).size().reset_index(name="count")
                ct_props = ct_counts.groupby("Sample", group_keys=False).apply(
                    lambda d: d.assign(prop=d["count"] / d["count"].sum())
                )
                df_pivot = (
                    ct_props.pivot(index="Sample", columns="CellType", values="prop")
                    .fillna(0.0)
                    .reindex(index=samples_xy, fill_value=0.0)
                    .reindex(columns=celltypes_order, fill_value=0.0)
                )
                xlocs = np.arange(len(df_pivot))
                xticklabels = [f"{z_lookup[s]:.0f}" for s in df_pivot.index]  # show z_height
            else:
                # XZ/YZ: n_bars bins over remaining axis (or z if forced)
                if plane == "xz":
                    rem_axis = sub["z_all"] if bars_over == "z" else sub["y_all"]
                else:  # "yz"
                    rem_axis = sub["z_all"] if bars_over == "z" else sub["x_all"]

                bar_idx = np.clip(np.digitize(rem_axis.to_numpy(), bar_edges, right=False) - 1, 0, n_bars - 1)
                sub_local = sub.copy()
                sub_local["BarBin"] = bar_idx

                ct_counts = sub_local.groupby(["BarBin", "CellType"]).size().reset_index(name="count")
                ct_props = ct_counts.groupby("BarBin", group_keys=False).apply(
                    lambda d: d.assign(prop=d["count"] / d["count"].sum())
                )
                df_pivot = (
                    ct_props.pivot(index="BarBin", columns="CellType", values="prop")
                    .reindex(index=np.arange(n_bars), fill_value=0.0)
                    .reindex(columns=celltypes_order, fill_value=0.0)
                    .fillna(0.0)
                )
                xlocs = np.arange(n_bars)
                xticklabels = bar_labels

            # --- Stacked bars ---
            bottom = np.zeros(len(df_pivot), dtype=float)
            for ct in celltypes_order:
                vals = df_pivot[ct].to_numpy()
                ax.bar(xlocs, vals, bottom=bottom, color=celltype_palette[ct], label=ct)
                bottom += vals

            # annotate medium/large segments
            for k in range(len(df_pivot)):
                cum = 0.0
                for ct in celltypes_order:
                    v = float(df_pivot.iloc[k][ct])
                    if v > annot_thresh:
                        ax.text(k, cum + v/2, f"{v:.2f}", ha="center", va="center", fontsize=6)
                    cum += v

            ax.set_ylim(0, 1)
            ax.set_xticks(xlocs)
            ax.set_xticklabels(xticklabels, rotation=90, fontsize=6)
            ax.set_title(f"{plane.upper()} Grid ({i+1},{j+1}) | n={len(sub)}", fontsize=9, pad=6)
            ax.tick_params(axis="y", labelsize=7)

    plt.savefig(
        plot_outdir / f"Br6660_cell_type_proportion_grid_{plane}.png", # name according to axis
        dpi=300,
        bbox_inches="tight",
    )
    plt.close("all")

plot_grid_props_plane("xy") # along z_height (slices)
plot_grid_props_plane("xz") # along y
plot_grid_props_plane("yz") # along x


############## Grid plot for a specific slice and grid proportion plot of a specific grid
### Panel B part 1: show grid on one representative slice
sample_key = "Sample"
selected_sample = "Br6660_NAc7_3580"   # change this to the representative slice we want
pt_size = 0.25
pt_alpha = 0.7

XY = aligned_adata.obsm[key_added]
x, y = XY[:, 0], XY[:, 1]
obs = aligned_adata.obs
samples_all = obs[sample_key].astype(str)
celltypes = obs[celltype_key].astype(str)
xbins = np.linspace(x.min(), x.max(), grid_size + 1)
ybins = np.linspace(y.min(), y.max(), grid_size + 1)
colors = celltypes.map(celltype_palette).to_numpy()

# Plot one representative slice
mask = samples_all.values == selected_sample

plt.ioff()
fig, ax = plt.subplots(figsize=(5.2, 5.2))

ax.scatter(
    x[mask],
    y[mask],
    s=pt_size,
    c=colors[mask],
    alpha=pt_alpha,
    rasterized=True,
)

# Overlay grid
for xv in xbins:
    ax.axvline(x=xv, color="black", lw=0.8, alpha=0.75)
for yv in ybins:
    ax.axhline(y=yv, color="black", lw=0.8, alpha=0.75)

# # Optional: label grid cells
# for i in range(grid_size):
#     for j in range(grid_size):
#         xc = 0.5 * (xbins[i] + xbins[i + 1])
#         yc = 0.5 * (ybins[j] + ybins[j + 1])
#         ax.text(
#             xc,
#             yc,
#             f"{i + 1},{j + 1}",
#             ha="center",
#             va="center",
#             fontsize=7,
#             color="black",
#             alpha=0.8,
#         )

ax.set_xlim(x.min(), x.max())
ax.set_ylim(y.min(), y.max())
ax.set_aspect("equal")
ax.set_xticks([])
ax.set_yticks([])
ax.set_title(f"5×5 spatial grid on representative slice\n{selected_sample}", fontsize=10)

plt.savefig(
    plot_outdir / f"{selected_sample}_grid_overlay.png",
    dpi=600,
    bbox_inches="tight",
)
plt.close("all")


### Panel B part 2: stacked bar plot of cell type proportions in one grid cell across slices (ordered by z_height)
grid_i = 3   
grid_j = 2  
annot_thresh = 0.05

df = pd.DataFrame({
    "x": x,
    "y": y,
    "Sample": samples_all.values,
    "CellType": celltypes.values,
    "z_height": obs["z_height"].values,
})
df["Xbin"] = np.clip(np.digitize(df["x"], xbins, right=False) - 1, 0, grid_size - 1)
df["Ybin"] = np.clip(np.digitize(df["y"], ybins, right=False) - 1, 0, grid_size - 1)

sub = df[(df["Xbin"] == grid_i) & (df["Ybin"] == grid_j)].copy()

# Slice order by z height
z_by_sample = (
    obs[[sample_key, "z_height"]]
    .groupby(sample_key)["z_height"]
    .mean()
    .sort_values()
)
samples_xy = z_by_sample.index.astype(str).tolist()

# Count and proportion by Sample × CellType
ct_counts = sub.groupby(["Sample", "CellType"]).size().reset_index(name="count")
ct_props = ct_counts.groupby("Sample", group_keys=False).apply(
    lambda d: d.assign(prop=d["count"] / d["count"].sum())
)

celltypes_order = list(celltype_palette.keys())

df_pivot = (
    ct_props.pivot(index="Sample", columns="CellType", values="prop")
    .fillna(0.0)
    .reindex(index=samples_xy, fill_value=0.0)
    .reindex(columns=celltypes_order, fill_value=0.0)
)

# Plot stacked proportions: horizontal version
plt.ioff()
fig, ax = plt.subplots(figsize=(5.2, 4.8))

ylocs = np.arange(len(df_pivot))
left = np.zeros(len(df_pivot), dtype=float)

for ct in celltypes_order:
    vals = df_pivot[ct].to_numpy()
    ax.barh(
        ylocs,
        vals,
        left=left,
        color=celltype_palette[ct],
        height=0.8,
        label=ct,
    )

    # optional annotation
    for k, v in enumerate(vals):
        if v >= annot_thresh:
            ax.text(
                left[k] + v / 2,
                k,
                f"{v:.2f}",
                ha="center",
                va="center",
                fontsize=7.5,
            )

    left += vals

ax.set_xlim(0, 1)
ax.set_xlabel("Cell-type proportion", fontsize=14)
ax.set_ylabel("Slice number", fontsize=14)
ax.tick_params(axis="both", labelsize=12)
ax.set_yticks(ylocs)
ax.set_yticklabels(
    [str(slice_number) for slice_number in range(1, len(df_pivot) + 1)],
    fontsize=12,
)

# Make top -> bottom = small -> large
ax.invert_yaxis()

ax.set_title(
    f"Cell-type composition in grid cell ({grid_i + 1}, {grid_j + 1})",
    fontsize=10
)

legend_order = sorted(celltypes_order, key=str.casefold)
legend_handles = [
    Patch(facecolor=celltype_palette[ct], edgecolor="none", label=ct)
    for ct in legend_order
]
ax.legend(
    handles=legend_handles,
    loc="center left",
    bbox_to_anchor=(1.02, 0.5),
    frameon=False,
    fontsize=8,
)

plt.savefig(
    plot_outdir / f"Br6660_grid_{grid_i + 1}_{grid_j + 1}_proportion_by_z_horizontal.png",
    dpi=600,
    bbox_inches="tight",
)
plt.close("all")
