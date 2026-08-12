#!/usr/bin/env python3
"""Build audited 2000-um primary and +1000-um shifted absolute-XY caches."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy import sparse

def args():
    p = argparse.ArgumentParser()
    p.add_argument('--metadata', type=Path, required=True)
    p.add_argument('--h5ad', type=Path, required=True)
    p.add_argument('--gene-filtering', type=Path, required=True)
    p.add_argument('--validated-grid-config', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--block-size', type=int, default=100000)
    return p.parse_args()

def sha(path, block=8 * 1024 * 1024):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda : f.read(block), b''):
            h.update(b)
    return h.hexdigest()

def text(ds):
    x = ds[:]
    if x.dtype.kind == 'S':
        return x.astype(str)
    return np.asarray(
        [v.decode() if isinstance(v, bytes) else str(v) for v in x]
    )

def csr_rows(group, rows):
    if group.attrs.get('encoding-type') not in ('csr_matrix', b'csr_matrix'):
        raise RuntimeError('H5AD X is not CSR')
    shape = tuple((int(v) for v in group.attrs['shape']))
    rows = np.asarray(rows, dtype=np.int64)
    (first, last) = (int(rows[0]), int(rows[-1]))
    ip = group['indptr'][first:last + 2].astype(np.int64)
    (lo, hi) = (int(ip[0]), int(ip[-1]))
    ip = ip - lo
    span = sparse.csr_matrix(
        (group['data'][lo:hi], group['indices'][lo:hi], ip),
        shape=(last - first + 1, shape[1]),
    )
    return span[rows - first, :].tocsr()

def make_grid(m, size, x0, y0, name):
    bx = np.floor((m.x.to_numpy() - x0) / size).astype(np.int64)
    by = np.floor((m.y.to_numpy() - y0) / size).astype(np.int64)
    keys = pd.DataFrame({'Sample': m.Sample, 'CellType': m.CellType, 'bx': bx, 'by': by})
    (gid, levels) = pd.factorize(pd.MultiIndex.from_frame(keys), sort=True)
    z = keys.assign(
        tile_group_id=gid,
        AP_um=m.z.to_numpy(),
        nucleus_area_sf=m['nucleus_area.sf'].to_numpy(),
    )
    meta = (
        z.groupby('tile_group_id', sort=True, observed=True)
        .agg(
            Sample=('Sample', 'first'),
            CellType=('CellType', 'first'),
            AP_um=('AP_um', 'first'),
            bx=('bx', 'first'),
            by=('by', 'first'),
            n_cells=('Sample', 'size'),
            exposure=('nucleus_area_sf', 'sum'),
        )
        .reset_index()
    )
    meta['absolute_tile'] = meta.bx.astype(str) + ':' + meta.by.astype(str)
    meta['ML_left_um'] = x0 + meta.bx * size
    meta['ML_right_um'] = meta.ML_left_um + size
    meta['DV_bottom_um'] = y0 + meta.by * size
    meta['DV_top_um'] = meta.DV_bottom_um + size
    meta['ML_center_um'] = (meta.ML_left_um + meta.ML_right_um) / 2
    meta['DV_center_um'] = (meta.DV_bottom_um + meta.DV_top_um) / 2
    meta['grid_variant'] = name
    meta['grid_origin_x_um'] = x0
    meta['grid_origin_y_um'] = y0
    if len(meta) != len(levels) or meta.tile_group_id.duplicated().any():
        raise RuntimeError('Invalid grid factorization')
    return dict(name=name, gid=gid, meta=meta, n=len(meta), x0=x0, y0=y0, size=size)

def atomic_csv(x, path, compression=None):
    tmp = path.with_name(path.name + f".tmp.{__import__('os').getpid()}")
    x.to_csv(tmp, index=False, compression=compression)
    tmp.replace(path)

def main():
    a = args()
    a.output_dir.mkdir(parents=True, exist_ok=True)
    (a.output_dir / 'checkpoints').mkdir(exist_ok=True)
    (a.output_dir / 'logs').mkdir(exist_ok=True)
    m = pd.read_csv(
        a.metadata,
        usecols=[
            'h5ad_row',
            'cell_id',
            'Sample',
            'CellType',
            'x',
            'y',
            'z',
            'nucleus_area.sf',
        ],
    )
    if (
        m.empty
        or m.h5ad_row.duplicated().any()
        or m.cell_id.duplicated().any()
        or m.isna().any().any()
    ):
        raise RuntimeError('Eligible metadata is empty, duplicated, or incomplete')
    if (
        not m.h5ad_row.is_monotonic_increasing
        or (np.diff(m.h5ad_row.to_numpy()) <= 0).any()
        or (m['nucleus_area.sf'] <= 0).any()
    ):
        raise RuntimeError('Invalid H5AD rows or exposure')
    vg = json.loads(a.validated_grid_config.read_text())
    old = float(vg['voxel_size_um'])
    if old != 500 or not np.isfinite([vg['x0'], vg['y0']]).all():
        raise RuntimeError('Validated 500-um grid origin unavailable')
    grids = [
        make_grid(m, 2000.0, float(vg['x0']), float(vg['y0']), 'primary'),
        make_grid(
            m,
            2000.0,
            float(vg['x0']) + 1000.0,
            float(vg['y0']) + 1000.0,
            'shifted_xy_plus1000um',
        ),
    ]
    filtering = pd.read_csv(a.gene_filtering)
    with h5py.File(a.h5ad, 'r') as h:
        index_key = h['var'].attrs.get('_index', '_index')
        if isinstance(index_key, bytes):
            index_key = index_key.decode()
        genes = text(h['var'][index_key])
        X = h['X']
        if (
            len(genes) != 366
            or not np.array_equal(genes.astype(str), filtering.gene.astype(str))
            or not filtering.tested.astype(str)
            .str.lower()
            .isin(['true', 't', '1'])
            .all()
        ):
            raise RuntimeError('Frozen 366-gene universe mismatch')
        for q in grids:
            q['counts'] = sparse.csr_matrix((q['n'], len(genes)), dtype=np.float64)
        total = 0.0
        rows = m.h5ad_row.to_numpy(np.int64)
        for start in range(0, len(rows), a.block_size):
            stop = min(start + a.block_size, len(rows))
            b = csr_rows(X, rows[start:stop])
            total += float(b.sum())
            for q in grids:
                A = sparse.csr_matrix(
                    (
                        np.ones(stop - start, dtype=np.int8),
                        (q['gid'][start:stop], np.arange(stop - start)),
                    ),
                    shape=(q['n'], stop - start),
                )
                q['counts'] += A @ b
            print(f'Aggregated raw counts for {stop:,}/{len(rows):,} retained cells', flush=True)
    files = {}
    recon = []
    for q in grids:
        d = q['counts'].toarray()
        if (d < 0).any() or not np.array_equal(d, np.rint(d)):
            raise RuntimeError(f"{q['name']} counts are invalid")
        d = np.rint(d).astype(np.int64)
        meta = q['meta']
        if (
            int(meta.n_cells.sum()) != len(m)
            or not np.isclose(
                meta.exposure.sum(),
                m['nucleus_area.sf'].sum(),
                rtol=1e-12,
                atol=1e-10,
            )
            or int(d.sum()) != round(total)
        ):
            raise RuntimeError(f"{q['name']} global reconciliation failed")
        mp = a.output_dir / f"{q['name']}_metadata.csv.gz"
        cp = a.output_dir / f"{q['name']}_counts.csv.gz"
        atomic_csv(meta, mp, 'gzip')
        ct = pd.DataFrame(d, columns=genes)
        ct.insert(0, 'tile_group_id', meta.tile_group_id)
        atomic_csv(ct, cp, 'gzip')
        files[mp.name] = sha(mp)
        files[cp.name] = sha(cp)
        src = (
            m.groupby('Sample', sort=False)
            .agg(
                input_cells=('cell_id', 'size'),
                input_exposure=('nucleus_area.sf', 'sum'),
            )
            .reset_index()
        )
        dst = (
            meta.groupby('Sample', sort=False)
            .agg(
                assigned_cells=('n_cells', 'sum'),
                assigned_exposure=('exposure', 'sum'),
                n_observed_slice_celltype_tiles=('tile_group_id', 'size'),
            )
            .reset_index()
        )
        r = src.merge(dst, on='Sample', validate='one_to_one')
        r['grid_variant'] = q['name']
        r['each_cell_assigned_exactly_once'] = r.input_cells.eq(r.assigned_cells)
        r['exposure_reconciled'] = np.isclose(
            r.input_exposure,
            r.assigned_exposure,
            rtol=1e-12,
            atol=1e-10,
        )
        recon.append(r)
    definition = pd.DataFrame(
        [
            dict(
                grid_variant=q['name'],
                tile_width_um=2000,
                tile_height_um=2000,
                grid_origin_x_um=q['x0'],
                grid_origin_y_um=q['y0'],
                ML_min_edge_um=q['meta'].ML_left_um.min(),
                ML_max_edge_um=q['meta'].ML_right_um.max(),
                DV_min_edge_um=q['meta'].DV_bottom_um.min(),
                DV_max_edge_um=q['meta'].DV_top_um.max(),
                boundary_convention=(
                    '[left,right) and [bottom,top); coordinates mapped by floor'
                ),
                coordinate_system='validated aligned final second-pass ML/DV',
                n_observed_slice_celltype_tiles=q['n'],
            )
            for q in grids
        ]
    )
    atomic_csv(definition, a.output_dir / '03_global_grid_definition_2000um.csv')
    audit = pd.concat(recon, ignore_index=True)
    atomic_csv(audit, a.output_dir / '04_cell_to_tile_reconciliation_2000um.csv')
    if not audit.each_cell_assigned_exactly_once.all() or not audit.exposure_reconciled.all():
        raise RuntimeError('Cell-to-tile reconciliation failed')
    atomic_csv(pd.DataFrame({'Gene': genes}), a.output_dir / 'genes.tsv')
    atomic_csv(filtering, a.output_dir / 'gene_filtering.csv')
    manifest = {
        'analysis_id': 'Br6660_absolute_xy_2000um_density_v1',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'n_cells': len(m),
        'n_genes': len(genes),
        'raw_count_total': round(total),
        'input_hashes': {
            'metadata': sha(a.metadata),
            'h5ad': sha(a.h5ad),
            'gene_filtering': sha(a.gene_filtering),
            'validated_grid_config': sha(a.validated_grid_config),
        },
        'grid_origin_inherited': {'x0': vg['x0'], 'y0': vg['y0']},
        'files': files,
        'software': {
            'python': platform.python_version(),
            'numpy': np.__version__,
            'pandas': pd.__version__,
            'h5py': h5py.__version__,
        },
    }
    (a.output_dir / 'INPUT_CACHE_MANIFEST.json').write_text(json.dumps(manifest, indent=2))
    (a.output_dir / 'checkpoints' / 'GRID_INPUT_COMPLETE.txt').write_text('PASS\n')
    print('2000-um primary and shifted absolute-grid input audit passed', flush=True)
if __name__ == '__main__':
    main()
