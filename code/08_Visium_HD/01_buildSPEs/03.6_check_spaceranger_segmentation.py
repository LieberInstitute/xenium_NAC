#   Visual check of Space Ranger (v4.0+) cell/nucleus segmentation, tiled into
#   a grid -- the Space Ranger analog of check_segmentations.py.
#
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D
import numpy as np
import os
import json
from pyhere import here
import pandas as pd
import tifffile
import datetime
import session_info

################################################################################
#   Pick the sample (same mechanism as the 03 script)
################################################################################

sample_info_path = here('processed-data', 'visiumHD_sample_info_NAc.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']
raw_image_path = sample_info.iloc[task_id]['raw_image_path']


################################################################################
#   Paths
################################################################################

#   Space Ranger segmented outputs (v4.0+)
seg_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'segmented_outputs'
)
nucleus_geojson = os.path.join(seg_dir, 'nucleus_segmentations.geojson')
cell_geojson = os.path.join(seg_dir, 'cell_segmentations.geojson')

#   Parent dir requested: segmentation_grid_spaceranger, with a per-sample sub
plot_dir = here(
    'plots', 'HD_Full_Analysis', 'segmentation_grid_spaceranger', sample_id
)
os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Settings
################################################################################

n_rows = 12        # tiles down the image (y) axis
n_cols = 15         # tiles across the image (x) axis
downscale = 1      # integer factor to shrink the full-res image (and polygon
                   # coords) for memory/speed. 1 = full res. Raise (2, 4, ...)
                   # if the overview is slow or the image won't fit in memory.
line_width = 0.5   # outline thickness

################################################################################
#   Load the full-resolution H&E (the image Space Ranger segmented)
################################################################################

print(f"{datetime.datetime.now()} | Reading {raw_image_path}")
he = np.asarray(tifffile.imread(raw_image_path))

#   Coerce to (H, W, 3)
if he.ndim > 3:
    he = np.squeeze(he)
if he.ndim == 3:
    if he.shape[0] in (3, 4) and he.shape[2] not in (3, 4):
        he = np.moveaxis(he, 0, -1)   # (C, H, W) -> (H, W, C)
    if he.shape[2] == 4:
        he = he[..., :3]              # drop alpha
if downscale > 1:
    he = he[::downscale, ::downscale]
H, W = he.shape[:2]
print(f"{datetime.datetime.now()} | Image shape (after downscale {downscale}): {he.shape}")

################################################################################
#   Load segmentation polygons from GeoJSON
################################################################################

def load_polys(path, downscale):
    """Return a list of exterior rings (Nx2 float arrays, in image pixels)."""
    if not os.path.exists(path):
        print(f"{datetime.datetime.now()} | WARNING: missing {path}")
        return []
    with open(path) as f:
        gj = json.load(f)
    polys = []
    for feat in gj.get('features', []):
        geom = feat.get('geometry') or {}
        gtype = geom.get('type')
        if gtype == 'Polygon':
            rings = [geom['coordinates'][0]]
        elif gtype == 'MultiPolygon':
            rings = [poly[0] for poly in geom['coordinates']]
        else:
            continue
        for ring in rings:
            polys.append(np.asarray(ring, dtype=float) / downscale)
    return polys


print(f"{datetime.datetime.now()} | Loading GeoJSON polygons")
nuc_polys = load_polys(nucleus_geojson, downscale)
cell_polys = load_polys(cell_geojson, downscale)
print(
    f"{datetime.datetime.now()} | {len(nuc_polys)} nucleus, "
    f"{len(cell_polys)} cell polygons"
)


def poly_bboxes(polys):
    """Return (N, 4) array of [xmin, xmax, ymin, ymax] per polygon."""
    if not polys:
        return np.empty((0, 4))
    return np.array([
        [p[:, 0].min(), p[:, 0].max(), p[:, 1].min(), p[:, 1].max()]
        for p in polys
    ])


nuc_bb = poly_bboxes(nuc_polys)
cell_bb = poly_bboxes(cell_polys)

#   Tissue extent for the grid: union of all polygon bounding boxes
all_bb = np.vstack([bb for bb in (cell_bb, nuc_bb) if bb.shape[0] > 0])
assert all_bb.shape[0] > 0, "No polygons found in either GeoJSON"
xmin, xmax = all_bb[:, 0].min(), all_bb[:, 1].max()
ymin, ymax = all_bb[:, 2].min(), all_bb[:, 3].max()

col_edges = np.linspace(xmin, xmax, n_cols + 1)   # x -> columns
row_edges = np.linspace(ymin, ymax, n_rows + 1)   # y -> rows


def overlap(bb, tx0, tx1, ty0, ty1):
    """Boolean mask of polygons whose bbox overlaps the tile box."""
    if bb.shape[0] == 0:
        return np.zeros(0, dtype=bool)
    return (
        (bb[:, 0] <= tx1) & (bb[:, 1] >= tx0) &
        (bb[:, 2] <= ty1) & (bb[:, 3] >= ty0)
    )

################################################################################
#   Overview: full H&E with the grid + tile labels
################################################################################

print(f"{datetime.datetime.now()} | Drawing grid overview")

fig_w = 12.0
fig, ax = plt.subplots(figsize=(fig_w, fig_w * H / W))
ax.imshow(he)
for r in range(n_rows):
    for c in range(n_cols):
        tx0, tx1 = col_edges[c], col_edges[c + 1]
        ty0, ty1 = row_edges[r], row_edges[r + 1]
        ax.add_patch(
            Rectangle(
                (tx0, ty0), tx1 - tx0, ty1 - ty0,
                fill=False, edgecolor='cyan', linewidth=0.8
            )
        )
        ax.text(
            (tx0 + tx1) / 2, (ty0 + ty1) / 2, f'{r + 1},{c + 1}',
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
#   Per tile: H&E crop with nucleus (red) and cell (blue) outlines
################################################################################

print(f"{datetime.datetime.now()} | Rendering {n_rows}x{n_cols} grid")

legend_handles = [
    Line2D([0], [0], color='#d62728', label='nucleus'),
    Line2D([0], [0], color='#1f77b4', label='cell'),
]

for r in range(n_rows):
    for c in range(n_cols):
        tx0, tx1 = col_edges[c], col_edges[c + 1]
        ty0, ty1 = row_edges[r], row_edges[r + 1]

        out_path = os.path.join(
            plot_dir, f'{sample_id}_segmentation_r{r + 1}_c{c + 1}.png'
        )

        #   Clamp the crop to the image
        x0c, x1c = int(max(0, np.floor(tx0))), int(min(W, np.ceil(tx1)))
        y0c, y1c = int(max(0, np.floor(ty0))), int(min(H, np.ceil(ty1)))

        fig, ax = plt.subplots(figsize=(6, 6))
        ax.imshow(he[y0c:y1c, x0c:x1c], extent=(x0c, x1c, y1c, y0c))

        sel_nuc = overlap(nuc_bb, tx0, tx1, ty0, ty1)
        sel_cell = overlap(cell_bb, tx0, tx1, ty0, ty1)
        if sel_nuc.any():
            ax.add_collection(LineCollection(
                [nuc_polys[i] for i in np.where(sel_nuc)[0]],
                colors='#d62728', linewidths=line_width
            ))
        if sel_cell.any():
            ax.add_collection(LineCollection(
                [cell_polys[i] for i in np.where(sel_cell)[0]],
                colors='#1f77b4', linewidths=line_width
            ))

        ax.set_xlim(x0c, x1c)
        ax.set_ylim(y1c, y0c)   # inverted y to match image orientation
        ax.axis('off')
        ax.set_title(f'{sample_id}  row {r + 1}/{n_rows}, col {c + 1}/{n_cols}')
        if sel_nuc.any() or sel_cell.any():
            ax.legend(handles=legend_handles, loc='upper right',
                      framealpha=0.7, fontsize=7)
        fig.savefig(out_path, dpi=150, bbox_inches='tight')
        plt.close('all')

print(
    f"{datetime.datetime.now()} | Done. Wrote 1 overview + {n_rows * n_cols} "
    f"tiles to {plot_dir}"
)


session_info.show()
