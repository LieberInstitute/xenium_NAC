#   Check segmentations across a whole tissue section by tiling it into a grid.

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from skimage.segmentation import find_boundaries
import scanpy as sc
import numpy as np
import os
from pyhere import here
import pandas as pd
import bin2cell as b2c
import datetime
import json
import session_info

################################################################################
#   Pick the sample (same mechanism as the 03 script)
################################################################################

sample_info_path = here('processed-data', 'visiumHD_sample_info_NAc.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']

################################################################################
#   Paths produced by the 03 script
################################################################################

stardist_dir = here(
    'processed-data', 'HD_Full_Analysis', 'bin2cell', 'stardist'
)

pre_out_path = here(
    'processed-data', 'HD_Full_Analysis', 'bin2cell',
    f'{sample_id}_pre_bin2cell.h5ad'
)

#   Each sample gets its own subfolder inside segmentation_grid
plot_dir = here('plots', 'HD_Full_Analysis', 'segmentation_grid', sample_id)
os.makedirs(plot_dir, exist_ok=True)

#   mpp, recomputed exactly as in the 03 script (needed by b2c.get_crop)
scalefactors_path = here(
    'processed-data', '01_spaceranger', sample_id, 'outs',
    'binned_outputs', 'square_002um', 'spatial', 'scalefactors_json.json'
)
with open(scalefactors_path) as f:
    mpp = json.load(f)['microns_per_pixel']

print(f"{datetime.datetime.now()} | mpp for {sample_id}: {mpp}")

################################################################################
#   Load the pre-aggregation (bin-level) object written at the end of 03
################################################################################

print(f"{datetime.datetime.now()} | Loading {pre_out_path}")
adata = sc.read_h5ad(pre_out_path)

#   H&E image + pixel coordinates of each bin, in the space sc.pl.spatial uses
spatial_key = list(adata.uns['spatial'].keys())[0]
img_key = list(adata.uns['spatial'][spatial_key]['images'].keys())[0]
he_img = adata.uns['spatial'][spatial_key]['images'][img_key]
H, W = he_img.shape[:2]
scalef = adata.uns['spatial'][spatial_key]['scalefactors'][
    f'tissue_{img_key}_scalef'
]
coords = np.asarray(adata.obsm['spatial_cropped_150_buffer']) * scalef

#   Per-bin source + assignment
adata.obs['labels_joint_source'] = (
    adata.obs['labels_joint_source'].astype('category')
)
sources = list(adata.obs['labels_joint_source'].cat.categories)
src_arr = adata.obs['labels_joint_source'].to_numpy()
assigned_arr = (adata.obs['labels_joint'] > 0).to_numpy()


def _pick_color(cat):
    """primary -> red, secondary -> blue."""
    s = str(cat).lower()
    if 'he' in s or 'primary' in s:
        return '#d62728'   # red
    if 'gex' in s or 'secondary' in s:
        return '#1f77b4'   # blue
    return None


_fallback = ['#d62728', '#1f77b4', '#2ca02c', '#9467bd', '#ff7f0e']
palette = {}
_fi = 0
for _cat in sources:
    _col = _pick_color(_cat)
    if _col is None:
        _col = _fallback[_fi % len(_fallback)]
        _fi += 1
    palette[_cat] = _col
print(f"{datetime.datetime.now()} | source -> color: {palette}")

################################################################################
#   Find cell outlines (boundary bins) once, over the whole tissue
################################################################################

#   Rasterize the joint labels onto the regular bin grid, then mark bins that
#   sit on the edge of their cell. Scattering only those draws cell outlines.
rr = adata.obs['array_row'].to_numpy().astype(int)
cc = adata.obs['array_col'].to_numpy().astype(int)
lj = adata.obs['labels_joint'].to_numpy().astype(np.int32)
r0, c0 = rr.min(), cc.min()
label_grid = np.zeros((rr.max() - r0 + 1, cc.max() - c0 + 1), dtype=np.int32)
label_grid[rr - r0, cc - c0] = lj
boundary_grid = find_boundaries(label_grid, mode='inner', background=0)
boundary_arr = boundary_grid[rr - r0, cc - c0]

################################################################################
#   Define the grid
################################################################################

#   Tile counts along each axis. More tiles -> smaller tiles -> more zoom.
n_rows = 12   # subdivisions down the array_row axis
n_cols = 15    # subdivisions across the array_col axis

row_edges = np.linspace(
    adata.obs['array_row'].min(), adata.obs['array_row'].max() + 1, n_rows + 1
)
col_edges = np.linspace(
    adata.obs['array_col'].min(), adata.obs['array_col'].max() + 1, n_cols + 1
)


def tile_mask(r, c):
    """Boolean mask of bins falling inside tile (r, c)."""
    return (
        (adata.obs['array_row'] >= row_edges[r]) &
        (adata.obs['array_row'] < row_edges[r + 1]) &
        (adata.obs['array_col'] >= col_edges[c]) &
        (adata.obs['array_col'] < col_edges[c + 1])
    )

################################################################################
#   Overview: full H&E with the grid + tile labels drawn on top
################################################################################

print(f"{datetime.datetime.now()} | Drawing grid overview on H&E")

fig_w = 12.0
fig, ax = plt.subplots(figsize=(fig_w, fig_w * H / W))
ax.imshow(he_img)

for r in range(n_rows):
    for c in range(n_cols):
        m = tile_mask(r, c).to_numpy()
        if not m.any():
            continue
        pc = coords[m]
        x0, y0 = pc[:, 0].min(), pc[:, 1].min()
        x1, y1 = pc[:, 0].max(), pc[:, 1].max()
        ax.add_patch(
            Rectangle(
                (x0, y0), x1 - x0, y1 - y0,
                fill=False, edgecolor='cyan', linewidth=0.8
            )
        )
        ax.text(
            (x0 + x1) / 2, (y0 + y1) / 2, f'{r + 1},{c + 1}',
            color='white', fontsize=6, ha='center', va='center',
            bbox=dict(facecolor='black', alpha=0.5, pad=1, edgecolor='none')
        )

ax.axis('off')
ax.set_title(f'{sample_id}  {n_rows}x{n_cols} tile index (row,col)')
overview_path = os.path.join(plot_dir, f'{sample_id}_grid_overview.png')
fig.savefig(overview_path, dpi=200, bbox_inches='tight')
plt.close('all')
print(f"{datetime.datetime.now()} | Wrote {overview_path}")

################################################################################
#   Per tile: H&E (source-colored outlines) + GEX image (secondary segmentation)
################################################################################

print(f"{datetime.datetime.now()} | Rendering {n_rows}x{n_cols} two-panel grid")

gex_tiff = os.path.join(stardist_dir, f'gex_{sample_id}.tiff')
gex_npz = os.path.join(stardist_dir, f'gex_{sample_id}.npz')

#   Marker size for the outline dots (points^2). Bump up for thicker outlines.
marker_size = 2

for r in range(n_rows):
    for c in range(n_cols):
        mask = tile_mask(r, c)
        m = mask.to_numpy()

        out_path = os.path.join(
            plot_dir, f'{sample_id}_segmentation_r{r + 1}_c{c + 1}.png'
        )

        fig, (axL, axR) = plt.subplots(1, 2, figsize=(12, 6))

        if not m.any():
            for a in (axL, axR):
                a.text(
                    0.5, 0.5, 'No in-tissue bins in this tile',
                    transform=a.transAxes, ha='center', va='center'
                )
                a.axis('off')
            fig.suptitle(f'{sample_id}  row {r + 1}/{n_rows}, col {c + 1}/{n_cols}')
            fig.savefig(out_path, dpi=150, bbox_inches='tight')
            plt.close('all')
            continue

        #-----------------------------------------------------------------------
        #   LEFT: H&E with source-colored cell outlines
        #-----------------------------------------------------------------------
        pc_tile = coords[m]
        x0 = int(np.floor(pc_tile[:, 0].min()))
        x1 = int(np.ceil(pc_tile[:, 0].max()))
        y0 = int(np.floor(pc_tile[:, 1].min()))
        y1 = int(np.ceil(pc_tile[:, 1].max()))
        x0c, x1c = max(0, x0), min(W, x1)
        y0c, y1c = max(0, y0), min(H, y1)

        axL.imshow(he_img[y0c:y1c, x0c:x1c], extent=(x0c, x1c, y1c, y0c))

        tile_outline = m & boundary_arr & assigned_arr
        for cat in sources:
            sel = tile_outline & (src_arr == cat)
            if not sel.any():
                continue
            pc = coords[sel]
            axL.scatter(
                pc[:, 0], pc[:, 1], s=marker_size, c=palette[cat],
                label=str(cat), linewidths=0, rasterized=True
            )
        axL.set_xlim(x0c, x1c)
        axL.set_ylim(y1c, y0c)
        axL.axis('off')
        axL.set_title('H&E (outlines by source)')
        if tile_outline.any():
            axL.legend(
                loc='upper right', markerscale=6, framealpha=0.7, fontsize=7
            )

        #-----------------------------------------------------------------------
        #   RIGHT: the GEX (sum_umi) image + its StarDist (secondary) labels.
        #   Uses bin2cell on the 'array' basis, exactly like the 03 script's
        #   secondary-segmentation plot, so the coordinate mapping is handled.
        #-----------------------------------------------------------------------
        try:
            crop = b2c.get_crop(adata[mask], basis="array", mpp=mpp)
            rendered = b2c.view_labels(
                image_path=gex_tiff,
                labels_npz_path=gex_npz,
                crop=crop,
                stardist_normalize=True
            )
            axR.imshow(rendered)
        except Exception as e:
            axR.text(
                0.5, 0.5, f'GEX view failed:\n{e}',
                transform=axR.transAxes, ha='center', va='center', fontsize=7
            )
        axR.axis('off')
        axR.set_title('GEX image (secondary segmentation)')

        fig.suptitle(f'{sample_id}  row {r + 1}/{n_rows}, col {c + 1}/{n_cols}')
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close('all')

print(
    f"{datetime.datetime.now()} | Done. Wrote 1 overview + {n_rows * n_cols} "
    f"two-panel tiles to {plot_dir}"
)



session_info.show()
