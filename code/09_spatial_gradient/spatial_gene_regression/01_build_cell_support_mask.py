#!/usr/bin/env python3
"""Build gene-independent analyzed-tissue support masks from all QC-passed cells."""
from __future__ import annotations
import argparse
import gzip
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import ndimage
from skimage.morphology import disk, remove_small_holes, remove_small_objects
from rasterio.features import shapes
from affine import Affine
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--validated-grid-config', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--plot-dir', type=Path, required=True)
    p.add_argument('--resolution-um', type=float, default=25.0)
    p.add_argument('--primary-buffer-um', type=float, default=100.0)
    p.add_argument('--sensitivity-buffer-um', type=float, nargs='*', default=[75.0, 125.0])
    p.add_argument('--fill-holes-max-mm2', type=float, default=0.25)
    p.add_argument('--remove-components-min-mm2', type=float, default=0.25)
    return p.parse_args()

def sha(p):
    h = hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda : f.read(8 * 1024 * 1024), b''):
            h.update(b)
    return h.hexdigest()

def atomic_csv(d, p, compression=None):
    q = p.with_name(p.name + f".tmp.{__import__('os').getpid()}")
    d.to_csv(q, index=False, compression=compression)
    q.replace(p)

def parameter_name(buffer_um, primary):
    return 'primary_100um_buffer' if primary else f'sensitivity_{int(buffer_um)}um_buffer'

def make_mask(x, y, x0, y0, res, buffer_um, hole_px, component_px):
    pad = buffer_um + 2 * res
    xmin = x0 + math.floor((x.min() - pad - x0) / res) * res
    ymin = y0 + math.floor((y.min() - pad - y0) / res) * res
    xmax = x0 + math.ceil((x.max() + pad - x0) / res) * res
    ymax = y0 + math.ceil((y.max() + pad - y0) / res) * res
    nx = int(round((xmax - xmin) / res))
    ny = int(round((ymax - ymin) / res))
    occ = np.zeros((ny, nx), dtype=bool)
    ix = np.clip(np.floor((x - xmin) / res).astype(int), 0, nx - 1)
    iy = np.clip(np.floor((y - ymin) / res).astype(int), 0, ny - 1)
    occ[iy, ix] = True
    mask = ndimage.binary_dilation(occ, structure=disk(int(round(buffer_um / res))))
    mask = remove_small_holes(mask, area_threshold=hole_px, connectivity=2)
    mask = remove_small_objects(mask, min_size=component_px, connectivity=2)
    return (mask, xmin, ymin)

def save_geojson(mask, xmin, ymin, res, path, props):
    feats = []
    transform = Affine(res, 0, xmin, 0, res, ymin)
    for geom, val in shapes(
        mask.astype(np.uint8),
        mask=mask,
        transform=transform,
        connectivity=8,
    ):
        if val == 1:
            feats.append({'type': 'Feature', 'properties': props, 'geometry': geom})
    obj = {
        'type': 'FeatureCollection',
        'name': 'cell_coordinate_derived_analyzed_tissue_support',
        'features': feats,
    }
    q = path.with_name(path.name + f".tmp.{__import__('os').getpid()}")
    with gzip.open(q, 'wt') as f:
        json.dump(obj, f, separators=(',', ':'))
    q.replace(path)

def grid_area(mask, xmin, ymin, res, sample, variant, gx, gy, param, props, mask_hash):
    (iy, ix) = np.nonzero(mask)
    xc = xmin + (ix + 0.5) * res
    yc = ymin + (iy + 0.5) * res
    bx = np.floor((xc - gx) / 2000.0).astype(np.int64)
    by = np.floor((yc - gy) / 2000.0).astype(np.int64)
    z = (
        pd.DataFrame({'bx': bx, 'by': by})
        .value_counts(sort=False)
        .rename('n_mask_pixels')
        .reset_index()
    )
    z.insert(0, 'grid_variant', variant)
    z.insert(0, 'Sample', sample)
    z['valid_tissue_area_mm2'] = z.n_mask_pixels * res * res / 1000000.0
    z['mask_parameter_set'] = param
    z['area_source'] = 'cell-coordinate-derived analyzed-tissue support mask'
    z['area_source_sha256'] = mask_hash
    z['independent_of_gene_and_celltype'] = True
    for (k, v) in props.items():
        z[k] = v
    return z

def overlay(mask, xmin, ymin, res, x, y, sample, variant, gx, gy, path, param):
    ymax = ymin + mask.shape[0] * res
    xmax = xmin + mask.shape[1] * res
    (fig, ax) = plt.subplots(figsize=(8, 8))
    ax.imshow(
        mask,
        origin='lower',
        extent=[xmin, xmax, ymin, ymax],
        cmap='Greys',
        alpha=0.55,
        interpolation='nearest',
    )
    step = max(1, len(x) // 30000)
    ax.scatter(x[::step], y[::step], s=0.25, c='#0072B2', alpha=0.25, rasterized=True)
    xs = np.arange(gx + math.floor((xmin - gx) / 2000) * 2000, xmax + 2000, 2000)
    ys = np.arange(gy + math.floor((ymin - gy) / 2000) * 2000, ymax + 2000, 2000)
    for v in xs:
        ax.axvline(v, color='#D55E00', lw=0.35, alpha=0.7)
    for v in ys:
        ax.axhline(v, color='#D55E00', lw=0.35, alpha=0.7)
    ax.set(
        xlabel='Aligned ML (µm)',
        ylabel='Aligned DV (µm)',
        title=f'{sample} | {variant}\ncell-coordinate support ({param})',
    )
    ax.set_aspect('equal')
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)

def main():
    a = parse_args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir / 'checkpoints').mkdir(exist_ok=True)
    maskdir = a.output_dir / 'cell_support_masks'
    maskdir.mkdir(exist_ok=True)
    a.plot_dir.mkdir(parents=True, exist_ok=True)
    vg = json.loads(a.validated_grid_config.read_text())
    x0 = float(vg['x0'])
    y0 = float(vg['y0'])
    res = a.resolution_um
    if 1000 % res or 2000 % res:
        raise RuntimeError(
            'Raster resolution must divide both 1000 and 2000 um so grid '
            'intersections are exact on raster boundaries'
        )
    m = pd.read_csv(a.metadata, usecols=['cell_id', 'Sample', 'x', 'y'])
    if m.empty or m.cell_id.duplicated().any() or m[['Sample', 'x', 'y']].isna().any().any():
        raise RuntimeError('QC-passed coordinate input invalid')
    px_mm2 = res * res / 1000000.0
    hole_px = max(1, int(math.floor(a.fill_holes_max_mm2 / px_mm2)))
    component_px = max(1, int(math.ceil(a.remove_components_min_mm2 / px_mm2)))
    params = [a.primary_buffer_um] + [
        v for v in a.sensitivity_buffer_um if v != a.primary_buffer_um
    ]
    areas = []
    summaries = []
    files = {}
    for (sample, d) in m.groupby('Sample', sort=False):
        x = d.x.to_numpy(float)
        y = d.y.to_numpy(float)
        for (k, buf) in enumerate(params):
            pname = parameter_name(buf, k == 0)
            (mask, xmin, ymin) = make_mask(x, y, x0, y0, res, buf, hole_px, component_px)
            props = {
                'mask_method': (
                    'rasterized buffered union of all QC-passed cell coordinates'
                ),
                'raster_resolution_um': res,
                'cell_buffer_um': buf,
                'fill_holes_max_mm2': a.fill_holes_max_mm2,
                'remove_components_min_mm2': a.remove_components_min_mm2,
                'connectivity': 8,
                'n_qc_passed_cells': len(d),
            }
            mp = maskdir / f'{sample}__{pname}.geojson.gz'
            save_geojson(
                mask,
                xmin,
                ymin,
                res,
                mp,
                {'Sample': sample, 'mask_parameter_set': pname, **props},
            )
            mh = sha(mp)
            files[str(mp.relative_to(a.output_dir))] = mh
            grid_variants = [
                ('primary', x0, y0),
                ('shifted_xy_plus1000um', x0 + 1000, y0 + 1000),
            ]
            for variant, gx, gy in grid_variants:
                ar = grid_area(mask, xmin, ymin, res, sample, variant, gx, gy, pname, props, mh)
                areas.append(ar)
                summaries.append(
                    dict(
                        Sample=sample,
                        grid_variant=variant,
                        mask_parameter_set=pname,
                        total_valid_area_mm2=ar.valid_tissue_area_mm2.sum(),
                        n_intersected_tiles=len(ar),
                        minimum_positive_tile_area_mm2=(
                            ar.valid_tissue_area_mm2.min()
                        ),
                        maximum_tile_area_mm2=ar.valid_tissue_area_mm2.max(),
                    )
                )
                if k == 0:
                    overlay(
                        mask,
                        xmin,
                        ymin,
                        res,
                        x,
                        y,
                        sample,
                        variant,
                        gx,
                        gy,
                        a.plot_dir / f'{sample}__{variant}.png',
                        pname,
                    )
            print(f'Built {pname} support for {sample}: {mask.sum() * px_mm2:.3f} mm2', flush=True)
    area = pd.concat(areas, ignore_index=True)
    atomic_csv(area, a.output_dir / 'cell_support_mask_area_2000um.csv.gz', 'gzip')
    summary = pd.DataFrame(summaries)
    base = summary[summary.mask_parameter_set == 'primary_100um_buffer'][
        ['Sample', 'grid_variant', 'total_valid_area_mm2']
    ].rename(columns={'total_valid_area_mm2': 'primary_total_valid_area_mm2'})
    summary = summary.merge(base, on=['Sample', 'grid_variant'])
    summary['relative_total_area_change_vs_primary'] = (
        summary.total_valid_area_mm2 - summary.primary_total_valid_area_mm2
    ) / summary.primary_total_valid_area_mm2
    atomic_csv(summary, a.output_dir / 'cell_support_mask_parameter_sensitivity.csv')
    ptab = pd.DataFrame(
        [
            dict(
                mask_parameter_set=parameter_name(v, i == 0),
                method=(
                    'rasterized buffered union of all QC-passed cell coordinates'
                ),
                raster_resolution_um=res,
                cell_buffer_um=v,
                fill_holes_max_mm2=a.fill_holes_max_mm2,
                remove_components_min_mm2=a.remove_components_min_mm2,
                connectivity=8,
                role=(
                    'primary density sensitivity'
                    if i == 0
                    else 'prespecified pilot mask-parameter sensitivity'
                ),
            )
            for i, v in enumerate(params)
        ]
    )
    atomic_csv(ptab, a.output_dir / 'cell_support_mask_parameters.csv')
    manifest = {
        'analysis_support': (
            'cell-coordinate-derived analyzed-tissue support mask; not '
            'DAPI-derived and not a histological tissue mask'
        ),
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'metadata_sha256': sha(a.metadata),
        'n_qc_passed_cells': len(m),
        'parameters': ptab.to_dict('records'),
        'files': files,
    }
    (a.output_dir / 'CELL_SUPPORT_MASK_MANIFEST.json').write_text(json.dumps(manifest, indent=2))
    (a.output_dir / 'checkpoints' / 'CELL_SUPPORT_MASK_COMPLETE.txt').write_text('PASS\n')
    print('Cell-coordinate-derived analyzed-tissue support masks passed', flush=True)
if __name__ == '__main__':
    main()
