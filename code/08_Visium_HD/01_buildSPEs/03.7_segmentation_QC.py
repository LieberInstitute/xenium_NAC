import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import os
import glob
from pyhere import here
import datetime
import session_info


BIN_AREA_UM2 = 4.0   # 2 um square bins -> 4 um^2 each


def equiv_diam_um(bins_per_cell):
    """Equivalent circle diameter (um) from a per-cell 2 um bin count."""
    area = np.asarray(bins_per_cell, dtype=float) * BIN_AREA_UM2
    return 2.0 * np.sqrt(area / np.pi)


################################################################################
#   Pick the sample (same mechanism as the 03 script)
################################################################################

sample_info_path = here('processed-data', 'visiumHD_sample_info_NAc.csv')
sample_info = pd.read_csv(sample_info_path)
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_id = sample_info.iloc[task_id]['sample_id']

################################################################################
#   Paths
################################################################################

b2c_pre = here(
    'processed-data', 'HD_Full_Analysis', 'bin2cell',
    f'{sample_id}_pre_bin2cell.h5ad'
)
b2c_final = here('processed-data', 'HD_Full_Analysis', 'bin2cell', f'{sample_id}.h5ad')

#   Space Ranger segmented outputs. barcode_mappings.parquet often sits in
#   outs/ rather than segmented_outputs/, so search the likely spots.
sr_outs = str(here('processed-data', '01_spaceranger', sample_id, 'outs'))


def _find_one(patterns):
    for pat in patterns:
        hits = sorted(glob.glob(pat, recursive=True))
        if hits:
            return hits[0]
    return None


sr_cell_h5 = _find_one([
    os.path.join(sr_outs, 'segmented_outputs', 'filtered_feature_cell_matrix.h5'),
    os.path.join(sr_outs, 'filtered_feature_cell_matrix.h5'),
    os.path.join(sr_outs, '**', 'filtered_feature_cell_matrix.h5'),
])
sr_bm = _find_one([
    os.path.join(sr_outs, 'segmented_outputs', 'barcode_mappings.parquet'),
    os.path.join(sr_outs, 'barcode_mappings.parquet'),
    os.path.join(sr_outs, '**', 'barcode_mappings.parquet'),
])
print(f"{datetime.datetime.now()} | SR cell matrix:      {sr_cell_h5}")
print(f"{datetime.datetime.now()} | SR barcode_mappings: {sr_bm}")

bins_002um_h5 = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
    'square_002um', 'filtered_feature_bc_matrix.h5'
)

out_dir = here('plots', 'HD_Full_Analysis', 'segmentation_qc_comparison', sample_id)
os.makedirs(out_dir, exist_ok=True)

################################################################################
#   Shared basis: in-tissue 2 um bins -> tissue area + per-barcode UMI
################################################################################

print(f"{datetime.datetime.now()} | Loading 2 um bins: {bins_002um_h5}")
m002 = sc.read_10x_h5(bins_002um_h5)
m002.var_names_make_unique()
n_bins_tissue = int(m002.n_obs)
tissue_area_mm2 = n_bins_tissue * BIN_AREA_UM2 / 1e6
umi002 = pd.Series(np.asarray(m002.X.sum(axis=1)).ravel(), index=m002.obs_names)
total_tissue_umi = float(umi002.sum())
print(
    f"{datetime.datetime.now()} | {n_bins_tissue} in-tissue 2um bins, "
    f"{tissue_area_mm2:.2f} mm^2"
)

#   Containers
metrics = {}        # {algo: {metric: value}}
dists = {}          # {algo: {'diam':..,'umi':..,'genes':..,'bins':..}}

################################################################################
#   bin2cell metrics
################################################################################

try:
    print(f"{datetime.datetime.now()} | bin2cell metrics")
    pre = sc.read_h5ad(b2c_pre)
    cell = sc.read_h5ad(b2c_final)
    sc.pp.calculate_qc_metrics(cell, inplace=True, percent_top=None)

    assigned = (pre.obs['labels_joint'] > 0).to_numpy()
    bin_umi = np.asarray(pre.X.sum(axis=1)).ravel()
    bins_per_cell = cell.obs['bin_count'].to_numpy().astype(float)
    umi_per_cell = cell.obs['total_counts'].to_numpy()
    genes_per_cell = cell.obs['n_genes_by_counts'].to_numpy()
    diam = equiv_diam_um(bins_per_cell)

    metrics['bin2cell'] = {
        'n_cells': int(cell.n_obs),
        'tissue_area_mm2': round(tissue_area_mm2, 3),
        'cell_density_per_mm2': round(cell.n_obs / tissue_area_mm2, 1),
        'frac_bins_assigned': round(float(assigned.mean()), 3),
        'frac_umi_assigned': round(float(bin_umi[assigned].sum() / bin_umi.sum()), 3),
        'median_umi_per_cell': float(np.median(umi_per_cell)),
        'median_genes_per_cell': float(np.median(genes_per_cell)),
        'median_bins_per_cell': float(np.median(bins_per_cell)),
        'median_equiv_diam_um': round(float(np.median(diam)), 2),
    }
    dists['bin2cell'] = {
        'diam': diam, 'umi': umi_per_cell,
        'genes': genes_per_cell, 'bins': bins_per_cell,
    }
except Exception as e:
    print(f"{datetime.datetime.now()} | bin2cell metrics FAILED: {e}")

################################################################################
#   Space Ranger metrics
################################################################################

try:
    print(f"{datetime.datetime.now()} | Space Ranger metrics")
    if sr_cell_h5 is None or sr_bm is None:
        raise FileNotFoundError(
            "Could not locate the Space Ranger cell matrix and/or "
            "barcode_mappings.parquet under outs/ or outs/segmented_outputs/"
        )
    sr_cell = sc.read_10x_h5(sr_cell_h5)
    sr_cell.var_names_make_unique()
    sc.pp.calculate_qc_metrics(sr_cell, inplace=True, percent_top=None)

    #   barcode_mappings column names vary in case (Space Ranger writes them
    #   lowercase: square_002um / cell_id / in_nucleus / in_cell). Resolve them.
    import pyarrow.parquet as pq
    bm_cols = pq.read_schema(sr_bm).names
    print(f"{datetime.datetime.now()} | barcode_mappings columns: {bm_cols}")

    def _col(*cands):
        low = {c.lower(): c for c in bm_cols}
        for cand in cands:
            if cand.lower() in low:
                return low[cand.lower()]
        raise KeyError(f"none of {cands} found in {bm_cols}")

    c_bc = _col('square_002um', 'barcode')
    c_cell = _col('cell_id')
    c_incell = _col('in_cell')

    bm = pd.read_parquet(sr_bm, columns=[c_bc, c_cell, c_incell])

    #   Coerce the in_cell flag to real booleans (handles bool / int / string)
    _ic = bm[c_incell]
    if _ic.dtype == bool:
        pass
    elif _ic.dtype.kind in 'iuf':
        bm[c_incell] = _ic.astype(bool)
    else:
        bm[c_incell] = _ic.astype(str).str.lower().isin(['true', '1', 't', 'yes'])

    #   Sanity check that 2um barcode names line up with the 2um matrix
    n_overlap = bm.set_index(c_bc).index.intersection(umi002.index).size
    print(f"{datetime.datetime.now()} | barcode overlap with 2um matrix: "
          f"{n_overlap} / {len(umi002)}")

    #   Fraction of in-tissue 2um bins / UMIs that landed in a cell
    in_cell = (
        bm.set_index(c_bc)[c_incell]
        .reindex(umi002.index).fillna(False).astype(bool)
    )
    frac_bins_assigned = float(in_cell.mean())
    frac_umi_assigned = float(umi002[in_cell.values].sum() / total_tissue_umi)

    #   Bins per cell from the barcode->cell mapping
    bins_per_cell = (
        bm.loc[bm[c_incell]].groupby(c_cell).size().to_numpy().astype(float)
    )
    umi_per_cell = sr_cell.obs['total_counts'].to_numpy()
    genes_per_cell = sr_cell.obs['n_genes_by_counts'].to_numpy()
    diam = equiv_diam_um(bins_per_cell)

    metrics['spaceranger'] = {
        'n_cells': int(sr_cell.n_obs),
        'tissue_area_mm2': round(tissue_area_mm2, 3),
        'cell_density_per_mm2': round(sr_cell.n_obs / tissue_area_mm2, 1),
        'frac_bins_assigned': round(frac_bins_assigned, 3),
        'frac_umi_assigned': round(frac_umi_assigned, 3),
        'median_umi_per_cell': float(np.median(umi_per_cell)),
        'median_genes_per_cell': float(np.median(genes_per_cell)),
        'median_bins_per_cell': float(np.median(bins_per_cell)),
        'median_equiv_diam_um': round(float(np.median(diam)), 2),
    }
    dists['spaceranger'] = {
        'diam': diam, 'umi': umi_per_cell,
        'genes': genes_per_cell, 'bins': bins_per_cell,
    }
except Exception as e:
    print(f"{datetime.datetime.now()} | Space Ranger metrics FAILED: {e}")
    print("  (check that segmented_outputs/ exists and the parquet columns "
          "are Square_002um / Cell_id / In_nucleus / In_cell)")

################################################################################
#   Comparison table
################################################################################

assert metrics, "No metrics computed for either algorithm"
table = pd.DataFrame(metrics)            # rows = metric, cols = algorithm
table.index.name = 'metric'

print(f"\n{datetime.datetime.now()} | Segmentation QC comparison for {sample_id}")
print(table.to_string())

csv_path = os.path.join(out_dir, f'{sample_id}_qc_comparison.csv')
table.to_csv(csv_path)
print(f"{datetime.datetime.now()} | Wrote {csv_path}")

################################################################################
#   Overlaid per-cell distributions
################################################################################

colors = {'bin2cell': '#d62728', 'spaceranger': '#1f77b4'}
fig, axes = plt.subplots(2, 2, figsize=(11, 8))
panels = [
    ('diam', 'equivalent diameter (um)', axes[0, 0], False),
    ('umi', 'log10(UMIs per cell + 1)', axes[0, 1], True),
    ('genes', 'genes per cell', axes[1, 0], False),
    ('bins', 'bins per cell', axes[1, 1], False),
]
for key, xlabel, ax, logx in panels:
    for algo, d in dists.items():
        vals = d[key]
        if logx:
            vals = np.log10(np.asarray(vals) + 1)
        ax.hist(vals, bins=60, density=True, histtype='step',
                linewidth=1.5, color=colors.get(algo), label=algo)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('density')
    ax.legend(fontsize=8)
fig.suptitle(f'{sample_id} segmentation QC: bin2cell vs Space Ranger')
fig.tight_layout()
png_path = os.path.join(out_dir, f'{sample_id}_qc_comparison.png')
fig.savefig(png_path, dpi=150, bbox_inches='tight')
plt.close('all')
print(f"{datetime.datetime.now()} | Wrote {png_path}")
print(f"{datetime.datetime.now()} | Done.")


session_info.show()
