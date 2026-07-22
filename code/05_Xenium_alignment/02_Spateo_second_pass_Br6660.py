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
import html
import math
import re
import scanpy as sc
import anndata as ad
import numpy as np
import subprocess
from pathlib import Path
import pandas as pd

git_root = Path(subprocess.check_output(
    ["git", "rev-parse", "--show-toplevel"], text=True).strip()
)

MODULE_NAME = "module02_Spateo_second_pass_Br6660"
OUTPUT_DIR = git_root / "processed-data" / "05_Xenium_alignment" / MODULE_NAME
PLOT_DIR = git_root / "plots" / "05_Xenium_alignment" / MODULE_NAME
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
PLOT_DIR.mkdir(parents=True, exist_ok=True)

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


def subplot_grid(n_panels):
    """Use the same compact, at-most-three-column layout as the 01 script."""
    if n_panels < 1:
        raise ValueError("At least one panel is required")
    n_columns = min(3, n_panels)
    n_rows = math.ceil(n_panels / n_columns)
    return plt.subplots(
        n_rows,
        n_columns,
        figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )


def plot_consecutive_overlays(
    plotted_slices,
    spatial_key,
    output_path,
    title_prefix,
    cell_type=None,
):
    """Plot second-pass slice pairs with the same style as the 01 script."""
    n_pairs = len(plotted_slices) - 1
    figure, axes = subplot_grid(n_pairs)
    axes_flat = axes.ravel()
    display_names = [display_sample_name(adata) for adata in plotted_slices]
    point_size = 0.6 if cell_type is not None else 0.2
    point_alpha = 0.8 if cell_type is not None else 0.45

    for pair_index, axis in enumerate(axes_flat[:n_pairs]):
        reference = plotted_slices[pair_index]
        moving = plotted_slices[pair_index + 1]
        reference_coordinates = np.asarray(reference.obsm[spatial_key])[:, :2]
        moving_coordinates = np.asarray(moving.obsm[spatial_key])[:, :2]

        if cell_type is not None:
            reference_mask = (
                reference.obs["CellType"].astype("string").eq(cell_type)
                .fillna(False).to_numpy(dtype=bool)
            )
            moving_mask = (
                moving.obs["CellType"].astype("string").eq(cell_type)
                .fillna(False).to_numpy(dtype=bool)
            )
            reference_coordinates = reference_coordinates[reference_mask]
            moving_coordinates = moving_coordinates[moving_mask]

        reference_name = display_names[pair_index]
        moving_name = display_names[pair_index + 1]
        axis.scatter(
            reference_coordinates[:, 0],
            reference_coordinates[:, 1],
            s=point_size,
            alpha=point_alpha,
            color="#4C78A8",
            linewidths=0,
            label=reference_name,
        )
        axis.scatter(
            moving_coordinates[:, 0],
            moving_coordinates[:, 1],
            s=point_size,
            alpha=point_alpha,
            color="#F58518",
            linewidths=0,
            label=moving_name,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(
            f"{title_prefix}{reference_name} → {moving_name}", fontsize=9
        )
        axis.legend(
            loc="upper right", fontsize=6, markerscale=8, frameon=False
        )

    for axis in axes_flat[n_pairs:]:
        axis.axis("off")

    figure.suptitle(
        "Second-pass Spateo consecutive-slice overlays"
        if cell_type is None
        else f"Second-pass Spateo overlays: {cell_type}",
        y=1.002,
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def add_brain_reference_axes(plotter):
    """Add a lower-left anatomical axis widget that rotates with the model."""
    plotter.add_axes(
        interactive=False,
        line_width=3,
        color="black",
        x_color="#D95F5F",
        y_color="#5FA35F",
        z_color="#5F7FD9",
        xlabel="ML (L-M)",
        ylabel="DV (D-V)",
        zlabel="AP (A-P)",
        viewport=(0.0, 0.0, 0.22, 0.22),
        label_size=(0.34, 0.12),
    )


def inject_interactive_cell_type(output_path, cell_type=None):
    """Add a fixed upper-left cell-type label to a by-cell-type HTML file."""
    if cell_type is None:
        return
    cell_type_overlay = f"""
<div style="position:fixed;left:18px;top:18px;z-index:1000;
            padding:8px 12px;background:rgba(255,255,255,0.88);
            border:1px solid #777;border-radius:4px;color:#111;
            font:600 18px Arial,sans-serif;">
  {html.escape(str(cell_type))}
</div>"""

    html_text = output_path.read_text(encoding="utf-8")
    if "</body>" not in html_text:
        raise ValueError(f"Exported HTML has no closing body tag: {output_path}")
    html_text = html_text.replace(
        "</body>",
        f"{cell_type_overlay}\n</body>",
        1,
    )
    output_path.write_text(html_text, encoding="utf-8")


def export_interactive_html(model, output_path, camera_position="iso", cell_type=None):
    """Export an interactive point cloud with fixed anatomical annotations."""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1500))
    try:
        plotter.set_background("white")
        plotter.add_mesh(
            model,
            scalars="annotation_rgba",
            rgba=True,
            style="points",
            render_points_as_spheres=True,
            point_size=3,
            ambient=0.2,
            smooth_shading=True,
            show_scalar_bar=False,
        )
        plotter.camera_position = camera_position
        add_brain_reference_axes(plotter)
        plotter.export_html(str(output_path))
        inject_interactive_cell_type(output_path, cell_type=cell_type)
    finally:
        plotter.close()

# Load the saved transformation from the initial pass
import pickle
transformation_path = (
    git_root
    / "processed-data"
    / "05_Xenium_alignment"
    / "module01_Spateo_initial_pass_Br6660"
    / "Spateo_transformation_Br6660.pkl"
)
with open(transformation_path, "rb") as f:
    transformation = pickle.load(f)

# Apply the transformation to the slices
spatial_key = 'spatial'
key_added = 'align_spatial'
aligned_slices = st.align.morpho_align_apply_transformation(
    models=slices,
    spatial_key=spatial_key,
    key_added=key_added,
    transformation=transformation,
)

# Re-align Br6660_NAc1_580 and Br6660_Nac2_1090
slice1 = slices[0].copy()
slice2 = slices[1].copy()
keep = {"Astro_A", "Astro_B"}

slice1_sub = slice1[slice1.obs["CellType"].astype(str).isin(keep)].copy()
slice2_sub = slice2[slice2.obs["CellType"].astype(str).isin(keep)].copy()

slice1_sub.obs["CellType"] = slice1_sub.obs["CellType"].astype("category")
slice2_sub.obs["CellType"] = slice2_sub.obs["CellType"].astype("category")

cur_transformation = st.align.morpho_align_transformation(
    models=[slice1_sub, slice2_sub],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    partial_robust_level=150,
    SVI_mode=False,
    rep_layer=["CellType"],
    rep_field=["obs"],
    dissimilarity=["label"]
)
transformation[0] = cur_transformation[0]


# Re-align Br6660_Nac3_1580 and Br6660_Nac4_2080
slice1 = slices[2].copy()
slice2 = slices[3].copy()
keep = {"Astro_A", "Astro_B"}

slice1_sub = slice1[slice1.obs["CellType"].astype(str).isin(keep)].copy()
slice2_sub = slice2[slice2.obs["CellType"].astype(str).isin(keep)].copy()

slice1_sub.obs["CellType"] = slice1_sub.obs["CellType"].astype("category")
slice2_sub.obs["CellType"] = slice2_sub.obs["CellType"].astype("category")

cur_transformation = st.align.morpho_align_transformation(
    models=[slice1_sub, slice2_sub],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    partial_robust_level=50,
    SVI_mode=False,
    rep_layer=["CellType"],
    rep_field=["obs"],
    dissimilarity=["label"]
)
transformation[2] = cur_transformation[0]


# Re-align Br6660_Nac4_2080 and Br6660_Nac5_2580
slice1 = slices[3].copy()
slice2 = slices[4].copy()
keep = {"Astro_A", "Astro_B"}

slice1_sub = slice1[slice1.obs["CellType"].astype(str).isin(keep)].copy()
slice2_sub = slice2[slice2.obs["CellType"].astype(str).isin(keep)].copy()

slice1_sub.obs["CellType"] = slice1_sub.obs["CellType"].astype("category")
slice2_sub.obs["CellType"] = slice2_sub.obs["CellType"].astype("category")

cur_transformation = st.align.morpho_align_transformation(
    models=[slice1_sub, slice2_sub],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    partial_robust_level=50,
    SVI_mode=False,
    rep_layer=["CellType"],
    rep_field=["obs"],
    dissimilarity=["label"]
)
transformation[3] = cur_transformation[0]


# Re-align Br6660_NAc7_3580 and Br6660_Nac10_4080
slice1 = slices[6].copy()
slice2 = slices[7].copy()
keep = {"Excitatory", "CHAT"}

slice1_sub = slice1[slice1.obs["CellType"].astype(str).isin(keep)].copy()
slice2_sub = slice2[slice2.obs["CellType"].astype(str).isin(keep)].copy()

cur_transformation = st.align.morpho_align_transformation(
    models=[slice1_sub, slice2_sub],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    partial_robust_level=50,
    SVI_mode=False,
)
transformation[6] = cur_transformation[0]


# Re-align Br6660_NAc8_4580 and Br6660_NAc9_5080
slice1 = slices[8].copy()
slice2 = slices[9].copy()
keep = {"Excitatory", "CHAT"}

slice1_sub = slice1[slice1.obs["CellType"].astype(str).isin(keep)].copy()
slice2_sub = slice2[slice2.obs["CellType"].astype(str).isin(keep)].copy()

slice1_sub.obs["CellType"] = slice1_sub.obs["CellType"].astype("category")
slice2_sub.obs["CellType"] = slice2_sub.obs["CellType"].astype("category")

cur_transformation = st.align.morpho_align_transformation(
    models=[slice1_sub, slice2_sub],
    spatial_key=spatial_key,
    key_added=key_added,
    device=device,
    verbose=True,
    partial_robust_level=50,
    SVI_mode=False,
    rep_layer=["CellType"],
    rep_field=["obs"],
    dissimilarity=["label"]
)
transformation[8] = cur_transformation[0]


# Apply the updated transformation to the slices
aligned_slices = st.align.morpho_align_apply_transformation(
    models=slices,
    spatial_key=spatial_key,
    key_added=key_added,
    transformation=transformation,
)

# Overlaid plot of aligned slices (second pass)
plt.ioff()
plot_consecutive_overlays(
    aligned_slices,
    key_added,
    PLOT_DIR / "Br6660_overlay_slices_2d_second_spateo.png",
    title_prefix="Second-pass aligned: ",
)


# Overlaid plot of aligned slices colored by cell type (second pass)
plt.ioff()

plot_outdir = PLOT_DIR / "by_celltype"
plot_outdir.mkdir(parents=True, exist_ok=True)

# get cell types, then sort alphabetically (case-insensitive)
celltypes = (
    adata_all.obs["CellType"].unique().tolist()
    if "adata_all" in globals()
    else sorted(set().union(*[ad.obs["CellType"].unique().tolist() for ad in aligned_slices]))
)
celltypes = sorted(celltypes, key=lambda x: str(x).lower())

for ct in celltypes:
    if not any(
        adata.obs["CellType"].astype("string").eq(ct).fillna(False).any()
        for adata in aligned_slices
    ):
        continue
    ct_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(ct))
    plot_consecutive_overlays(
        aligned_slices,
        key_added,
        plot_outdir / f"{ct_safe}_overlay_slices_2d.png",
        title_prefix=f"{ct}: ",
        cell_type=ct,
    )

save_path = OUTPUT_DIR / "Spateo_transformation_Br6660.pkl"
with open(save_path, "wb") as f:
    pickle.dump(transformation, f, protocol=pickle.HIGHEST_PROTOCOL)

import anndata as ad
aligned_adata = ad.concat(aligned_slices)
pd.unique(aligned_adata.obs["z_height"]) # array([ 580, 1090, 1580, 2080, 2580, 3080, 3580, 4080, 4580, 5080, 5580])
aligned_adata.obsm['spatial_3D'] = np.concatenate([aligned_adata.obsm['align_spatial'], np.array(aligned_adata.obs['z_height'].values)[:,None]], axis=1)

# Build one table for all aligned slices
aligned_coord_df = pd.DataFrame(
    {
        "donor": DONOR,
        "sample": aligned_adata.obs["Sample"].astype(str).values if "Sample" in aligned_adata.obs.columns else "",
        "cell_id": aligned_adata.obs["cell_id"].astype(str).values if "cell_id" in aligned_adata.obs.columns else aligned_adata.obs_names.astype(str),
        "cell_type": aligned_adata.obs["CellType"].astype(str).values if "CellType" in aligned_adata.obs.columns else "",
        "x_second_pass": aligned_adata.obsm['spatial_3D'][:, 0],
        "y_second_pass": aligned_adata.obsm['spatial_3D'][:, 1],
        "z_height": aligned_adata.obsm['spatial_3D'][:, 2],
    }
)

# Save updated coordinates to csv
csv_path = OUTPUT_DIR / "Br6660_second_pass_coordinates.csv"
aligned_coord_df.to_csv(csv_path, index=False)

# Build 3D reconstruction of aligned slices
point_cloud, _ = st.tdr.construct_pc(adata=aligned_adata,spatial_key="spatial_3D",groupby='CellType', key_added="annotation", colormap=celltype_palette)
export_interactive_html(
    point_cloud,
    PLOT_DIR / "Br6660_reconstruction_3d_second_spateo.html",
    camera_position="iso",
)
## Flip the z axis for better visualization
aligned_adata.obsm['spatial_3D'][:, 2] *= -1
point_cloud_flipped, _ = st.tdr.construct_pc(adata=aligned_adata,spatial_key="spatial_3D",groupby='CellType', key_added="annotation", colormap=celltype_palette)
export_interactive_html(
    point_cloud_flipped,
    PLOT_DIR / "Br6660_reconstruction_3d_second_spateo_flipped_z.html",
    camera_position="iso",
)

# Make plots of aligned slices colored by cell type using three views
########## The commented part below is using the original color code for each cell type
# cluster_key = 'CellType'
# highlight_tissues = ['D1_Island_A']
# pc_highlight = point_cloud_flipped.copy()
# pc_highlight['annotation_rgba'][:,3] = 0.003
# pc_highlight['annotation_rgba'][aligned_adata.obs[cluster_key].isin(highlight_tissues),3] = 0.2
##########

## 2d plots colored by cell type
cluster_key = "CellType"
outdir = PLOT_DIR / "by_celltype_3d_static"
outdir.mkdir(parents=True, exist_ok=True)
pv.Plotter.export_vtkjs = lambda self, filename: self.screenshot(filename)
celltypes = pd.unique(aligned_adata.obs[cluster_key])
celltypes = sorted(celltypes, key=lambda x: str(x).lower())

for ct in celltypes:
    print(f"Generating plot for: {ct}")
    pc_highlight = point_cloud_flipped.copy()
    rgba = pc_highlight["annotation_rgba"].copy()
    highlight_mask = aligned_adata.obs[cluster_key].isin([ct]).to_numpy()
    rgba[:, 0:3] = np.array([0.7, 0.7, 0.7])   # grey RGB
    rgba[:, 3] = 0.01                          # low opacity for background
    # restore original color for highlighted cells
    rgba[highlight_mask, 0:3] = point_cloud_flipped["annotation_rgba"][highlight_mask, 0:3]
    rgba[highlight_mask, 3] = 0.8              # higher opacity for highlighted cells
    pc_highlight["annotation_rgba"] = rgba
    ct_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(ct))
    st.pl.three_d_multi_plot(
        model=st.tdr.collect_models([pc_highlight, pc_highlight, pc_highlight]),
        key="annotation",
        model_style="points",
        model_size=3,
        cpo=["xy", "xz", "yz"],
        jupyter="static",   
        off_screen=True,
        window_size=(1500, 1500),
        text=[f"{ct} xy", f"{ct} xz", f"{ct} yz"],
        text_kwargs={"font_size": 12, "text_loc": "upper_edge"},
        show_legend=False,
        plotter_filename=str(outdir / f"{ct_safe}_reconstruction_3d.png"),
    )

print(f"Generating plot for: D1_Island_A and D1_Island_B")
pc_highlight = point_cloud_flipped.copy()
rgba = pc_highlight["annotation_rgba"].copy()
highlight_mask = aligned_adata.obs[cluster_key].isin(["D1_Island_A", "D1_Island_B"]).to_numpy()
rgba[:, 0:3] = np.array([0.7, 0.7, 0.7])   # grey RGB
rgba[:, 3] = 0.01                          # low opacity for background
# restore original color for highlighted cells
rgba[highlight_mask, 0:3] = point_cloud_flipped["annotation_rgba"][highlight_mask, 0:3]
rgba[highlight_mask, 3] = 0.8              # higher opacity for highlighted cells
pc_highlight["annotation_rgba"] = rgba
ct_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", "D1_Island_A_and_D1_Island_B")
st.pl.three_d_multi_plot(
    model=st.tdr.collect_models([pc_highlight, pc_highlight, pc_highlight]),
    key="annotation",
    model_style="points",
    model_size=3,
    cpo=["xy", "xz", "yz"],
    jupyter="static",   
    off_screen=True,
    window_size=(1500, 1500),
    text=[f"D1_Island_A_and_D1_Island_B xy", f"D1_Island_A_and_D1_Island_B xz", f"D1_Island_A_and_D1_Island_B yz"],
    text_kwargs={"font_size": 12, "text_loc": "upper_edge"},
    show_legend=False,
    plotter_filename=str(outdir / f"{ct_safe}_reconstruction_3d.png"),
)

## 3d plots colored by cell type and save in html
delattr(pv.Plotter, "export_vtkjs")
outdir = PLOT_DIR / "by_celltype_3d_interactive"
outdir.mkdir(parents=True, exist_ok=True)

for ct in celltypes:
    print(f"Generating interative html plot for: {ct}")
    pc_highlight = point_cloud_flipped.copy()
    rgba = pc_highlight["annotation_rgba"].copy()
    highlight_mask = aligned_adata.obs[cluster_key].isin([ct]).to_numpy()
    rgba[:, 0:3] = np.array([0.7, 0.7, 0.7])   # grey RGB
    rgba[:, 3] = 0.01                          # low opacity for background
    # restore original color for highlighted cells
    rgba[highlight_mask, 0:3] = point_cloud_flipped["annotation_rgba"][highlight_mask, 0:3]
    rgba[highlight_mask, 3] = 0.8              # higher opacity for highlighted cells
    pc_highlight["annotation_rgba"] = rgba
    ct_safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(ct))
    export_interactive_html(
        pc_highlight,
        outdir / f"{ct_safe}_reconstruction_3d.html",
        camera_position="xy",
        cell_type=ct,
    )
