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

from pygam import LogisticGAM, te
import plotly.express as px
import plotly.graph_objects as go
from scipy.spatial import Delaunay


git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

MODULE_NAME = "cell_type_pseudotime"
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
# Flip z so slice order is top-to-bottom in 3D plots (plotly's default
# orientation puts higher z up, but our z_height is increasing along the 
# tissue block from anterior to posterior — flipping makes the rendered 
# stack match the anatomical viewing order)
aligned_adata.obs['z_height'] *= -1

aligned_adata.obsm['spatial_3D'] = np.concatenate([
    aligned_adata.obsm['align_spatial'],
    np.array(aligned_adata.obs['z_height'].values)[:, None]
], axis=1)

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



### Prepare data for 3D voxelization and GAM fitting
obs = aligned_adata.obs.copy()
coords = np.asarray(aligned_adata.obsm["align_spatial"])

obs["x_aligned"] = coords[:, 0]
obs["y_aligned"] = coords[:, 1]

obs = obs.loc[
    obs["CellType"].notna() &
    np.isfinite(obs["x_aligned"]) &
    np.isfinite(obs["y_aligned"])
].copy()


def make_3d_voxel_composition(
    obs,
    x_col="x_aligned",
    y_col="y_aligned",
    z_col="z_height",
    celltype_col="CellType",
    n_bins_x=30,
    n_bins_y=30,
    min_total_cells_per_voxel=20,
    min_total_cells_per_type=50,
    use_quantile_bins=True,
    verbose=True,
):
    df = obs[[x_col, y_col, z_col, celltype_col]].copy()
    n_cells_in = len(df)

    # Filter rare cell types
    ct_counts = df[celltype_col].value_counts()
    keep_ct = ct_counts[ct_counts >= min_total_cells_per_type].index
    drop_ct = ct_counts[ct_counts < min_total_cells_per_type]

    if verbose:
        print(f"=== Cell-type filter (>= {min_total_cells_per_type} cells) ===")
        print(f"  Kept:    {len(keep_ct)} types ({ct_counts[keep_ct].sum():,} cells)")
        print(f"  Dropped: {len(drop_ct)} types")
        for ct, n in drop_ct.items():
            print(f"    - {ct}: {n} cells")

    df = df.loc[df[celltype_col].isin(keep_ct)].copy()

    # Bin in x and y 
    if use_quantile_bins:
        df["x_bin"] = pd.qcut(df[x_col], q=n_bins_x, duplicates="drop")
        df["y_bin"] = pd.qcut(df[y_col], q=n_bins_y, duplicates="drop")
    else:
        df["x_bin"] = pd.cut(df[x_col], bins=n_bins_x, include_lowest=True)
        df["y_bin"] = pd.cut(df[y_col], bins=n_bins_y, include_lowest=True)

    voxel_info = (
        df.groupby(["x_bin", "y_bin", z_col], observed=True)
          .agg(
              x_min=(x_col, "min"),
              x_max=(x_col, "max"),
              y_min=(y_col, "min"),
              y_max=(y_col, "max"),
              n_total=(celltype_col, "size"),
          )
          .reset_index()
    )
    voxel_info["x_center"] = 0.5 * (voxel_info["x_min"] + voxel_info["x_max"])
    voxel_info["y_center"] = 0.5 * (voxel_info["y_min"] + voxel_info["y_max"])
    voxel_info["z_center"] = voxel_info[z_col].astype(float)
    voxel_info = voxel_info.drop(columns=["x_min", "x_max", "y_min", "y_max"])

    # Filter sparse voxels
    n_voxels_pre = len(voxel_info)
    voxel_info = voxel_info.loc[voxel_info["n_total"] >= min_total_cells_per_voxel].copy()
    n_voxels_post = len(voxel_info)

    if verbose:
        print(f"\n=== Voxel filter (>= {min_total_cells_per_voxel} cells/voxel) ===")
        print(f"  Occupied voxels: {n_voxels_pre:,}")
        print(f"  Kept:            {n_voxels_post:,}")
        print(f"  Dropped:         {n_voxels_pre - n_voxels_post:,}")
        print(f"  Cells retained:  {voxel_info['n_total'].sum():,} / {n_cells_in:,} "
              f"({100 * voxel_info['n_total'].sum() / n_cells_in:.1f}%)")

    # Cell-type counts per valid voxel
    df = df.merge(voxel_info[["x_bin", "y_bin", z_col]],
                  on=["x_bin", "y_bin", z_col], how="inner")

    count_df = (
        df.groupby(["x_bin", "y_bin", z_col, celltype_col], observed=True)
          .size()
          .rename("n_ct")
          .reset_index()
    )

    # Fill zero counts for absent cell types
    full = voxel_info[["x_bin", "y_bin", z_col]].merge(
        pd.DataFrame({celltype_col: sorted(keep_ct.tolist())}),
        how="cross",
    )
    count_df = full.merge(count_df, on=["x_bin", "y_bin", z_col, celltype_col], how="left")
    count_df["n_ct"] = count_df["n_ct"].fillna(0).astype(int)

    out = count_df.merge(
        voxel_info[["x_bin", "y_bin", z_col, "x_center", "y_center", "z_center", "n_total"]],
        on=["x_bin", "y_bin", z_col], how="left",
    )
    out["prop"] = out["n_ct"] / out["n_total"]
    return out


def fit_3d_gam_per_celltype(
    df_voxel,
    celltype_col="CellType",
    min_voxels_per_type=30,
    n_splines=(6, 6, 4),
    lam=1.0,
    verbose=True,
):
    fit_dict = {}
    pred_tables = []
    skipped = []

    all_ct = sorted(df_voxel[celltype_col].unique())

    if verbose:
        print(f"=== GAM fitting (>= {min_voxels_per_type} occupied voxels) ===")

    for ct in all_ct:
        d = df_voxel.loc[df_voxel[celltype_col] == ct].copy()
        n_occupied = int((d["n_ct"] > 0).sum())
        n_cells = int(d["n_ct"].sum())

        if n_occupied < min_voxels_per_type:
            skipped.append(ct)
            if verbose:
                print(f"  SKIP  {ct:<20s} occupied={n_occupied:>5d}  cells={n_cells:>7,d}")
            continue

        # standardize for numerical stability
        X0 = d[["x_center", "y_center", "z_center"]].to_numpy(dtype=float)
        mu, sd = X0.mean(axis=0), X0.std(axis=0)
        X = (X0 - mu) / sd

        # grouped binomial → weighted 0/1 rows
        X_fit = np.vstack([X, X])
        y_fit = np.concatenate([np.ones(len(d), dtype=int), np.zeros(len(d), dtype=int)])
        w_fit = np.concatenate([
            d["n_ct"].to_numpy(dtype=float),
            (d["n_total"] - d["n_ct"]).to_numpy(dtype=float),
        ])
        keep = w_fit > 0
        X_fit, y_fit, w_fit = X_fit[keep], y_fit[keep], w_fit[keep]

        gam = LogisticGAM(te(0, 1, 2, n_splines=list(n_splines)), lam=lam)
        gam.fit(X_fit, y_fit, weights=w_fit)
        # Sanity checks
        print(gam.terms)         # shows the term structure
        print(gam.coef_.shape)   # should be 6*6*4 + 1 = 145 if list was respected

        out = d.copy()
        out["fit"] = gam.predict_proba(X)
        pred_tables.append(out)
        fit_dict[ct] = {"gam": gam, "mean": mu, "std": sd}

        if verbose:
            print(f"  FIT   {ct:<20s} occupied={n_occupied:>5d}  cells={n_cells:>7,d}")

    if verbose:
        print(f"\nFitted: {len(fit_dict)}/{len(all_ct)} cell types")
        if skipped:
            print(f"Skipped: {skipped}")

    pred_df = pd.concat(pred_tables, ignore_index=True) if pred_tables else pd.DataFrame()
    return fit_dict, pred_df


df_voxel = make_3d_voxel_composition(
    obs,
    n_bins_x=50, n_bins_y=50,
    min_total_cells_per_voxel=20,
    min_total_cells_per_type=30,
    use_quantile_bins=True,
)

fit_3d, pred_3d = fit_3d_gam_per_celltype(
    df_voxel,
    min_voxels_per_type=30,
    n_splines=(6, 6, 4),
    lam=1.0,
)

# === Cell-type filter (>= 30 cells) ===
#   Kept:    20 types (2,425,841 cells)
#   Dropped: 0 types

# === Voxel filter (>= 20 cells/voxel) ===
#   Occupied voxels: 21,787
#   Kept:            21,612
#   Dropped:         175
#   Cells retained:  2,424,288 / 2,425,841 (99.9%)
# === GAM fitting (>= 30 occupied voxels) ===
#   FIT   Astro_A              occupied=15224  cells=169,995
#   FIT   Astro_B              occupied=17871  cells=202,868
#   FIT   Astrocyte_Oligo      occupied=19006  cells= 77,901
#   FIT   CHAT                 occupied= 5406  cells= 18,935
#   FIT   D1_Island_A          occupied=11928  cells= 60,217
#   FIT   D1_Island_B          occupied= 5501  cells= 42,045
#   FIT   DRD1_MSN             occupied=14009  cells=128,631
#   FIT   DRD2_MSN             occupied=13924  cells=130,569
#   FIT   Ependymal            occupied=  269  cells= 12,109
#   FIT   Excitatory           occupied= 4474  cells= 70,078
#   FIT   Fibroblast_A         occupied=21250  cells=165,531
#   FIT   Fibroblast_B         occupied=19822  cells= 99,425
#   FIT   Inh_PVALB            occupied=17986  cells=117,648
#   FIT   Inh_SST              occupied= 8055  cells= 18,000
#   FIT   MSN_Oligo            occupied=15524  cells= 99,423
#   FIT   Microglia_A          occupied=21535  cells=169,570
#   FIT   Microglia_B          occupied=10946  cells= 31,364
#   FIT   Microglia_Oligo      occupied=14859  cells= 33,878
#   FIT   OPC                  occupied=21319  cells=131,744
#   FIT   WM                   occupied=20684  cells=644,357

# Fitted: 20/20 cell types


### Plot 3D spatial proportion or fitted value by cell type
def plot_3d_celltype_proportion(
    pred_df,
    celltype,
    value_col="fit",                # "fit" or "prop"
    celltype_col="CellType",
    x_col="x_center",
    y_col="y_center",
    z_col="z_center",
    size_by_value=False,
    title=None,
    save_path=None,                 # if given, save instead of show
):
    """
    3D scatter:
      axes  = aligned (x, y, z)
      color = fitted/observed proportion
    """
    d = pred_df.loc[pred_df[celltype_col] == celltype].copy()
    if d.empty:
        print(f"  No data for {celltype}, skipping")
        return

    # point size
    if size_by_value:
        v = d[value_col].to_numpy()
        vmin, vmax = float(np.min(v)), float(np.max(v))
        d["_size"] = 5 + 20 * (v - vmin) / (vmax - vmin) if vmax > vmin else 8.0
    else:
        d["_size"] = 6.0

    # hover only on columns that exist
    hover_cols = {
        x_col: ":.1f", y_col: ":.1f", z_col: ":.1f",
        value_col: ":.4f",
    }
    for c in ("n_ct", "n_total"):
        if c in d.columns:
            hover_cols[c] = True

    fig = px.scatter_3d(
        d,
        x=x_col, y=y_col, z=z_col,
        color=value_col,
        size="_size",
        size_max=18,
        color_continuous_scale="Viridis",
        hover_data=hover_cols,
        title=title or f"{celltype}: 3D spatial proportion ({value_col})",
    )
    fig.update_traces(marker=dict(opacity=0.85))
    fig.update_layout(
        scene=dict(xaxis_title="Aligned x",
                   yaxis_title="Aligned y",
                   zaxis_title="Depth z"),
        width=900, height=750,
    )

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_path, include_plotlyjs="cdn")
    else:
        fig.show()


# Run for all (or top-N) cell types, saving prop and fit versions 
plot3d_outdir = plot_outdir / "3d_voxel_plots"
plot3d_outdir.mkdir(parents=True, exist_ok=True)

top_ct = sorted(fit_3d.keys())

for value_col in ("prop", "fit"):
    print(f"=== Saving {value_col} plots ===")
    for ct in sorted(top_ct):
        out_path = plot3d_outdir / f"Br6660_3d_{value_col}_{ct}.html"
        plot_3d_celltype_proportion(
            pred_3d,
            celltype=ct,
            value_col=value_col,
            size_by_value=True,
            save_path=out_path,
        )
        print(f"  saved: {out_path.name}")


# ### smooth 3D heatmap of fitted proportion
def plot_celltype_3d_heatmap(
    fit_dict,
    df_voxel,
    celltype,
    nx=60, ny=60, nz=30,
    isomin_quantile=0.5,
    isomax_quantile=1.0,
    opacity=0.1,
    surface_count=20,
    colorscale="Viridis",
    width=950, height=800,
    save_path=None,
    mask_outside_tissue=True,
    verbose=True,
):
    """
    Smooth 3D heatmap of fitted cell-type proportion.

    GAM is evaluated directly on a regular grid (no griddata interpolation).
    Grid points outside the tissue convex hull are excluded so that
    rendered values reflect the GAM inside the tissue, not extrapolation.
    """
    if celltype not in fit_dict:
        print(f"  No GAM for {celltype}, skipping")
        return

    info = fit_dict[celltype]
    gam, mu, sd = info["gam"], info["mean"], info["std"]

    # Regular grid in original coordinates
    d = df_voxel.loc[df_voxel["CellType"] == celltype]
    xg = np.linspace(d["x_center"].min(), d["x_center"].max(), nx)
    yg = np.linspace(d["y_center"].min(), d["y_center"].max(), ny)
    zg = np.linspace(d["z_center"].min(), d["z_center"].max(), nz)
    X, Y, Z = np.meshgrid(xg, yg, zg, indexing="ij")

    # Predict on grid (standardize using fit-time mu/sd)
    grid_pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    V = gam.predict_proba((grid_pts - mu) / sd).reshape(X.shape)
    v_flat = V.ravel().copy()

    # Mask grid points outside tissue convex hull (avoid extrapolation)
    inside_mask = None
    if mask_outside_tissue:
        voxel_pts = (
            d.loc[d["n_total"] > 0, ["x_center", "y_center", "z_center"]]
             .to_numpy()
        )
        try:
            hull = Delaunay(voxel_pts)
            inside_mask = hull.find_simplex(grid_pts) >= 0
            n_kept = int(inside_mask.sum())
            if verbose:
                print(f"  {celltype}: kept {n_kept:,}/{len(grid_pts):,} "
                      f"({100*n_kept/len(grid_pts):.1f}%) inside tissue")
        except Exception as e:
            print(f"  {celltype}: hull construction failed ({e}), using full grid")

    # Color range from inside-tissue points only
    valid = v_flat[inside_mask] if inside_mask is not None else v_flat
    if valid.size == 0:
        print(f"  {celltype}: no valid points after masking, skipping")
        return

    if verbose:
        print(f"    p̂ percentiles inside tissue: "
              f"50%={np.quantile(valid, 0.5):.3f}, "
              f"75%={np.quantile(valid, 0.75):.3f}, "
              f"90%={np.quantile(valid, 0.9):.3f}, "
              f"99%={np.quantile(valid, 0.99):.3f}, "
              f"max={valid.max():.3f}")

    isomin = float(np.quantile(valid, isomin_quantile))
    isomax = float(np.quantile(valid, isomax_quantile))

    # Sentinel value far below isomin so plotly skips outside-tissue points
    if inside_mask is not None:
        v_flat[~inside_mask] = -1e9

    fig = go.Figure(go.Volume(
        x=X.ravel(), y=Y.ravel(), z=Z.ravel(),
        value=v_flat,
        isomin=isomin,
        isomax=isomax,
        opacity=opacity,
        surface_count=surface_count,
        colorscale=colorscale,
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorbar=dict(title="p̂"),
    ))
    fig.update_layout(
        title=f"{celltype}: smoothed 3D fitted proportion heatmap "
              f"(p̂ inside tissue: {valid.min():.3f}–{valid.max():.3f})",
        scene=dict(
            xaxis_title="Aligned x",
            yaxis_title="Aligned y",
            zaxis_title="Depth z",
            aspectmode="data",
        ),
        width=width, height=height,
    )

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_path, include_plotlyjs="cdn")
    else:
        fig.show()


# Run for all fitted cell types
volume_outdir = plot_outdir / "3d_heatmap_plots"
volume_outdir.mkdir(parents=True, exist_ok=True)

for ct in sorted(fit_3d.keys()):
    plot_celltype_3d_heatmap(
        fit_3d, df_voxel, celltype=ct,
        nx=60, ny=60, nz=20,
        isomin_quantile=0.6,
        opacity=0.12,
        surface_count=30,
        save_path=volume_outdir / f"Br6660_3d_heatmap_{ct}.html",
    )
    print(f"  saved: {ct}\n")


### Parameter selection for Figure 2 panel C: 3D heatmap for D1_Island_A and D1_Island_B
sweep_outdir = plot_outdir / "3d_heatmap_sweep_d1_islands_Br6660"
sweep_outdir.mkdir(parents=True, exist_ok=True)

target_celltypes = ["D1_Island_A", "D1_Island_B", "WM", "Excitatory"]

# Quick nz check (1 config, 4 nz values, 1 cell type) 
print("=== Quick nz comparison ===")
for nz in [11, 15, 20, 30]:
    plot_celltype_3d_heatmap(
        fit_3d, df_voxel, celltype="D1_Island_A",
        nx=60, ny=60, nz=nz,
        isomin_quantile=0.75,
        opacity=0.12,
        surface_count=30,
        save_path=sweep_outdir / f"nz_compare__nz{nz:02d}.html",
        verbose=False,
    )
    print(f"  saved: nz_compare__nz{nz:02d}.html")

# Main sweep: isomin_quantile × opacity × surface_count
sweep_params = []
for isomin_q in [0.5, 0.6, 0.7, 0.8, 0.9]:
    for opacity in [0.12]:
        for surface_count in [30]:
            sweep_params.append({
                "isomin_quantile": isomin_q,
                "opacity": opacity,
                "surface_count": surface_count,
            })

print(f"\n=== Main sweep: {len(sweep_params)} configs × {len(target_celltypes)} cell types "
      f"= {len(sweep_params) * len(target_celltypes)} plots ===")

for ct in target_celltypes:
    print(f"\n--- {ct} ---")
    for p in sweep_params:
        tag = (f"q{int(p['isomin_quantile']*100)}"
               f"_op{int(p['opacity']*100):02d}"
               f"_sc{p['surface_count']:02d}")
        out_path = sweep_outdir / f"{ct}_{tag}.html"
        plot_celltype_3d_heatmap(
            fit_3d, df_voxel, celltype=ct,
            nx=60, ny=60, nz=20,           # pinned after nz check
            isomin_quantile=p["isomin_quantile"],
            isomax_quantile=1.0,
            opacity=p["opacity"],
            surface_count=p["surface_count"],
            save_path=out_path,
            verbose=False,
        )
        print(f"  saved: {out_path.name}\n")


def plot_celltype_3d_heatmap_clean(
    fit_dict,
    df_voxel,
    celltype,
    nx=60, ny=60, nz=20,
    isomin_quantile=0.6,
    isomax_quantile=1.0,
    opacity=0.12,
    surface_count=30,
    colorscale="Viridis",
    width=950, height=800,
    save_path=None,
    mask_outside_tissue=True,
    verbose=True,
):
    """
    Figure-quality 3D heatmap with axes / grid / background removed.
    Same algorithm as plot_celltype_3d_heatmap, only the layout differs.
    """
    if celltype not in fit_dict:
        print(f"  No GAM for {celltype}, skipping")
        return

    info = fit_dict[celltype]
    gam, mu, sd = info["gam"], info["mean"], info["std"]

    d = df_voxel.loc[df_voxel["CellType"] == celltype]
    xg = np.linspace(d["x_center"].min(), d["x_center"].max(), nx)
    yg = np.linspace(d["y_center"].min(), d["y_center"].max(), ny)
    zg = np.linspace(d["z_center"].min(), d["z_center"].max(), nz)
    X, Y, Z = np.meshgrid(xg, yg, zg, indexing="ij")

    grid_pts = np.column_stack([X.ravel(), Y.ravel(), Z.ravel()])
    V = gam.predict_proba((grid_pts - mu) / sd).reshape(X.shape)
    v_flat = V.ravel().copy()

    inside_mask = None
    if mask_outside_tissue:
        voxel_pts = (
            d.loc[d["n_total"] > 0, ["x_center", "y_center", "z_center"]]
             .to_numpy()
        )
        try:
            hull = Delaunay(voxel_pts)
            inside_mask = hull.find_simplex(grid_pts) >= 0
            if verbose:
                n_kept = int(inside_mask.sum())
                print(f"  {celltype}: kept {n_kept:,}/{len(grid_pts):,} "
                      f"({100*n_kept/len(grid_pts):.1f}%) inside tissue")
        except Exception as e:
            print(f"  {celltype}: hull construction failed ({e}), using full grid")

    valid = v_flat[inside_mask] if inside_mask is not None else v_flat
    if valid.size == 0:
        print(f"  {celltype}: no valid points after masking, skipping")
        return

    isomin = float(np.quantile(valid, isomin_quantile))
    isomax = float(np.quantile(valid, isomax_quantile))

    if inside_mask is not None:
        v_flat[~inside_mask] = -1e9

    fig = go.Figure(go.Volume(
        x=X.ravel(), y=Y.ravel(), z=Z.ravel(),
        value=v_flat,
        isomin=isomin,
        isomax=isomax,
        opacity=opacity,
        surface_count=surface_count,
        colorscale=colorscale,
        caps=dict(x_show=False, y_show=False, z_show=False),
        colorbar=dict(title="p̂"),
    ))

    # Hide tick labels & grid but keep axis lines
    hidden_axis = dict(
        showbackground=True,
        showgrid=True,
        showticklabels=False,
        zeroline=False,
        title="",
        # visible=True (default), showline=True (default) — axis lines stay
    )

    fig.update_layout(
        title=celltype,
        scene=dict(
            xaxis=hidden_axis,
            yaxis=hidden_axis,
            zaxis=hidden_axis,
            aspectmode="data",
        ),
        width=width, height=height,
        margin=dict(l=0, r=0, t=40, b=0),   # tight margins for figure use
    )

    if save_path is not None:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.write_html(save_path, include_plotlyjs="cdn")
    else:
        fig.show()


# Run for D1_Island_A and D1_Island_B with fixed figure parameters
figure_outdir = plot_outdir / "3d_heatmap_figure"
figure_outdir.mkdir(parents=True, exist_ok=True)

for ct in ["D1_Island_A", "D1_Island_B", "Excitatory", "WM"]:
    plot_celltype_3d_heatmap_clean(
        fit_3d, df_voxel, celltype=ct,
        nx=60, ny=60, nz=20,
        isomin_quantile=0.6,
        opacity=0.12,
        surface_count=30,
        save_path=figure_outdir / f"Br6660_3d_heatmap_figure_{ct}.html",
    )
    print(f"  saved figure: {ct}\n")
