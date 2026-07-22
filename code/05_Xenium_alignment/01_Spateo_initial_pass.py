#!/usr/bin/env python3
"""Run the initial Spateo alignment pass for an ordered set of donor slices.

The alignment uses the complete expression matrix in ``adata.X`` together
with the configured spatial coordinates. Cell-type annotations are attached
from a CSV for QC and plotting, but they do not define the alignment features.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata as md
import json
import logging
import math
import os
import pickle
import platform
import re
import shutil
import sys
import tempfile
import time
import warnings
from pathlib import Path
from typing import Any, Sequence

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ.setdefault("NUMBA_CACHE_DIR", f"/tmp/numba_spateo_{os.getuid()}")

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
from numba.core.errors import NumbaWarning


LOGGER = logging.getLogger("module02_spateo_initial_pass")
SCRIPT_VERSION = "1.0"

CELLTYPE_PALETTE = {
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


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path)
    parser.add_argument("--celltype-csv", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument("--donor")
    parser.add_argument("--donor-column", default="Donor")
    parser.add_argument("--sample-column", default="Sample")
    parser.add_argument("--sample-order", nargs="+")
    parser.add_argument("--adata-cell-id-column", default="cell_id")
    parser.add_argument("--celltype-csv-id-column", default="cell_id")
    parser.add_argument("--celltype-csv-column", default="CellType")
    parser.add_argument("--celltype-column", default="CellType")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--key-added", default="align_spatial")
    parser.add_argument("--device", choices=("cuda", "cpu"), default="cuda")
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--dtype", choices=("float32", "float64"), default="float32")
    parser.add_argument("--partial-robust-level", type=float, default=30.0)
    parser.add_argument("--chunk-capacity", type=float, default=2.0)
    parser.add_argument(
        "--sparse-calculation-mode",
        action=argparse.BooleanOptionalAction,
        default=True,
    )
    parser.add_argument(
        "--use-chunk", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument(
        "--flip-y", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument(
        "--verbose", action=argparse.BooleanOptionalAction, default=True
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--run-self-tests", action="store_true")
    args = parser.parse_args(argv)

    required = (
        "input_h5ad", "celltype_csv", "output_dir", "plot_dir", "donor",
        "sample_order",
    )
    if not args.run_self_tests and any(getattr(args, name) is None for name in required):
        parser.error(
            "--input-h5ad, --celltype-csv, --output-dir, --plot-dir, --donor, "
            "and --sample-order are required unless --run-self-tests is used"
        )
    if args.sample_order and len(args.sample_order) != len(set(args.sample_order)):
        parser.error("--sample-order contains duplicate names")
    if args.sample_order and len(args.sample_order) < 2:
        parser.error("--sample-order requires at least two ordered slices")
    if args.max_iter < 1:
        parser.error("--max-iter must be positive")
    if args.partial_robust_level <= 0:
        parser.error("--partial-robust-level must be positive")
    if args.chunk_capacity <= 0:
        parser.error("--chunk-capacity must be positive")
    if args.resume and args.overwrite:
        parser.error("--resume and --overwrite are mutually exclusive")
    if (
        not args.run_self_tests
        and args.output_dir.resolve() == args.plot_dir.resolve()
    ):
        parser.error("--output-dir and --plot-dir must be different directories")
    return args


def safe_filename(value: str) -> str:
    """Convert an annotation value into a filesystem-safe component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "unnamed"


def file_sha256(path: Path) -> str:
    """Return a stable content fingerprint for a file."""
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def transformation_path(output_dir: Path, donor: str) -> Path:
    """Return the donor-specific transformation output path."""
    return output_dir / f"Spateo_transformation_{safe_filename(donor)}.pkl"


def first_pass_coordinates_path(output_dir: Path, donor: str) -> Path:
    """Return the donor-specific first-pass coordinate-table path."""
    return output_dir / f"{safe_filename(donor)}_first_pass_coordinates.csv"


def metadata_path(output_dir: Path) -> Path:
    """Return the run-metadata output path."""
    return output_dir / "module02_initial_pass_metadata.json"


def validate_device(device: str, cuda_available: bool) -> None:
    """Reject an unavailable explicitly requested CUDA device."""
    if device == "cuda" and not cuda_available:
        raise RuntimeError(
            "--device cuda was requested, but torch.cuda.is_available() is False. "
            "Run on a GPU node with a CUDA-enabled environment or use --device cpu."
        )


def initialize_spateo_runtime() -> tuple[Any, Any, str]:
    """Initialize headless rendering and import the heavy runtime libraries."""
    import pyvista as pv

    try:
        pv.start_xvfb()
    except Exception:
        pass

    warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")
    warnings.filterwarnings("ignore", category=NumbaWarning)

    import spateo as st
    import torch

    try:
        st.__version__ = md.version("spateo-release")
    except md.PackageNotFoundError:
        pass
    return st, torch, str(getattr(st, "__version__", "unknown"))


def flip_y_inplace(adata: ad.AnnData, spatial_key: str) -> float:
    """Reflect one slice around its own mid-Y and return that midpoint."""
    coordinates = np.asarray(adata.obsm[spatial_key], dtype=float)
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError(
            f"adata.obsm[{spatial_key!r}] must have shape (n_cells, >=2)"
        )
    if not np.isfinite(coordinates[:, :2]).all():
        raise ValueError(f"adata.obsm[{spatial_key!r}] contains nonfinite coordinates")
    midpoint = 0.5 * (coordinates[:, 1].min() + coordinates[:, 1].max())
    coordinates = coordinates.copy()
    coordinates[:, 1] = 2.0 * midpoint - coordinates[:, 1]
    adata.obsm[spatial_key] = coordinates
    return float(midpoint)


def prepare_donor_slices(
    adata: ad.AnnData,
    annotations: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[list[ad.AnnData], pd.DataFrame, pd.DataFrame]:
    """Subset one donor, attach cell types, and build anatomically ordered slices."""
    required_obs = [args.donor_column, args.sample_column, args.adata_cell_id_column]
    missing_obs = [column for column in required_obs if column not in adata.obs]
    if missing_obs:
        raise ValueError(f"Missing required adata.obs column(s): {missing_obs}")
    if args.spatial_key not in adata.obsm:
        raise ValueError(f"Missing adata.obsm[{args.spatial_key!r}]")

    required_csv = [args.celltype_csv_id_column, args.celltype_csv_column]
    missing_csv = [column for column in required_csv if column not in annotations]
    if missing_csv:
        raise ValueError(f"Missing required cell-type CSV column(s): {missing_csv}")

    donor_mask = (
        adata.obs[args.donor_column].astype("string").eq(str(args.donor))
        .fillna(False).to_numpy(dtype=bool)
    )
    if not donor_mask.any():
        observed = sorted(
            adata.obs[args.donor_column].dropna().astype(str).unique().tolist()
        )
        raise ValueError(f"Donor {args.donor!r} is absent; observed donors: {observed}")

    donor_adata = adata[donor_mask].copy()
    donor_ids = donor_adata.obs[args.adata_cell_id_column].astype("string")
    if donor_ids.isna().any() or donor_ids.duplicated().any():
        raise ValueError(
            f"adata.obs[{args.adata_cell_id_column!r}] must be nonmissing and "
            "unique within the selected donor"
        )

    annotation_ids = annotations[args.celltype_csv_id_column].astype("string")
    if annotation_ids.isna().any() or annotation_ids.duplicated().any():
        raise ValueError(
            f"CSV column {args.celltype_csv_id_column!r} must be nonmissing and unique"
        )
    annotation_values = annotations[args.celltype_csv_column].astype("string")
    mapping = pd.Series(annotation_values.to_numpy(), index=annotation_ids.to_numpy())
    mapped = donor_ids.map(mapping)
    missing_labels = mapped.isna()
    annotation_qc = pd.DataFrame.from_records(
        [{
            "donor": str(args.donor),
            "n_cells_selected_donor": donor_adata.n_obs,
            "n_annotation_rows": len(annotations),
            "n_cells_with_celltype": int((~missing_labels).sum()),
            "n_cells_missing_celltype": int(missing_labels.sum()),
            "fraction_cells_with_celltype": float((~missing_labels).mean()),
        }]
    )
    if missing_labels.any():
        examples = donor_ids[missing_labels].astype(str).head(10).tolist()
        raise ValueError(
            f"Cell-type CSV has no annotation for {int(missing_labels.sum())} "
            f"selected donor cells; example IDs: {examples}"
        )

    donor_adata.obs[args.celltype_column] = pd.Categorical(mapped.astype(str))
    donor_adata.obs_names = donor_ids.astype(str).to_numpy()
    donor_adata.obsm[args.spatial_key] = np.asarray(
        getattr(donor_adata.obsm[args.spatial_key], "values", donor_adata.obsm[args.spatial_key]),
        dtype=float,
    )

    sample_values = donor_adata.obs[args.sample_column].astype("string")
    observed_samples = set(sample_values.dropna().astype(str))
    missing_samples = [sample for sample in args.sample_order if sample not in observed_samples]
    if missing_samples:
        raise ValueError(
            f"Configured samples absent from donor {args.donor}: {missing_samples}"
        )

    slices: list[ad.AnnData] = []
    slice_qc_rows: list[dict[str, Any]] = []
    for sample in args.sample_order:
        sample_mask = sample_values.eq(sample).fillna(False).to_numpy(dtype=bool)
        slice_adata = donor_adata[sample_mask].copy()
        z_match = re.search(r"_(\d+)$", sample)
        if z_match is None:
            raise ValueError(
                f"Sample {sample!r} does not end with an integer z height"
            )
        z_height = int(z_match.group(1))
        slice_adata.obs["z_height"] = z_height
        midpoint = np.nan
        if args.flip_y:
            midpoint = flip_y_inplace(slice_adata, args.spatial_key)
        else:
            coordinates = np.asarray(slice_adata.obsm[args.spatial_key], dtype=float)
            if not np.isfinite(coordinates[:, :2]).all():
                raise ValueError(
                    f"Slice {sample} contains nonfinite spatial coordinates"
                )
        coordinates = np.asarray(slice_adata.obsm[args.spatial_key])[:, :2]
        slices.append(slice_adata)
        slice_qc_rows.append(
            {
                "slice": sample,
                "z_height": z_height,
                "n_cells": slice_adata.n_obs,
                "n_celltypes": slice_adata.obs[args.celltype_column].nunique(),
                "x_min": float(coordinates[:, 0].min()),
                "x_max": float(coordinates[:, 0].max()),
                "y_min": float(coordinates[:, 1].min()),
                "y_max": float(coordinates[:, 1].max()),
                "y_axis_flipped": bool(args.flip_y),
                "y_flip_reflection_midpoint": midpoint,
            }
        )
    return slices, annotation_qc, pd.DataFrame.from_records(slice_qc_rows)


def build_palette(celltypes: Sequence[str]) -> dict[str, Any]:
    """Extend the established palette deterministically for unseen annotations."""
    palette = dict(CELLTYPE_PALETTE)
    missing = sorted(set(map(str, celltypes)) - set(palette), key=str.lower)
    if missing:
        cmap = plt.get_cmap("tab20", len(missing))
        for index, cell_type in enumerate(missing):
            palette[cell_type] = matplotlib.colors.to_hex(cmap(index))
    return palette


def alignment_metadata(args: argparse.Namespace) -> dict[str, Any]:
    """Return parameters that must match before a transformation is reused."""
    input_stat = args.input_h5ad.stat()
    return {
        "script_version": SCRIPT_VERSION,
        "input_h5ad": str(args.input_h5ad.resolve()),
        "input_h5ad_size": input_stat.st_size,
        "input_h5ad_mtime_ns": input_stat.st_mtime_ns,
        "celltype_csv": str(args.celltype_csv.resolve()),
        "celltype_csv_sha256": file_sha256(args.celltype_csv),
        "donor": str(args.donor),
        "donor_column": args.donor_column,
        "sample_column": args.sample_column,
        "sample_order": list(args.sample_order),
        "adata_cell_id_column": args.adata_cell_id_column,
        "celltype_csv_id_column": args.celltype_csv_id_column,
        "celltype_csv_column": args.celltype_csv_column,
        "celltype_column": args.celltype_column,
        "spatial_key": args.spatial_key,
        "key_added": args.key_added,
        "flip_y": args.flip_y,
        "expression_representation": "all_genes_from_adata_X",
        "device": args.device,
        "max_iter": args.max_iter,
        "dtype": args.dtype,
        "partial_robust_level": args.partial_robust_level,
        "chunk_capacity": args.chunk_capacity,
        "sparse_calculation_mode": args.sparse_calculation_mode,
        "use_chunk": args.use_chunk,
    }


def metadata_mismatches(saved: dict[str, Any], expected: dict[str, Any]) -> list[str]:
    """Describe all transformation metadata differences."""
    return [
        f"{key}: saved={saved.get(key)!r}, current={value!r}"
        for key, value in expected.items()
        if saved.get(key) != value
    ]


def prepare_output_dirs(args: argparse.Namespace) -> None:
    """Create separate output directories and enforce overwrite protection."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    expected_nonplots = [
        transformation_path(args.output_dir, args.donor),
        metadata_path(args.output_dir),
        first_pass_coordinates_path(args.output_dir, args.donor),
        args.output_dir / "celltype_annotation_qc.csv",
        args.output_dir / "slice_qc.csv",
        args.output_dir / "session_info.txt",
        args.output_dir / "runtime_summary.csv",
    ]
    plot_files = list(args.plot_dir.rglob("*.png"))
    if args.overwrite:
        for path in expected_nonplots:
            if path.is_file():
                path.unlink()
        for path in plot_files:
            path.unlink()
        celltype_plot_dir = args.plot_dir / "by_celltype"
        if celltype_plot_dir.exists() and not any(celltype_plot_dir.iterdir()):
            celltype_plot_dir.rmdir()
        return
    if args.resume:
        required_resume = expected_nonplots[:2]
        missing = [str(path) for path in required_resume if not path.is_file()]
        if missing:
            raise FileNotFoundError(
                "--resume requires an existing transformation and metadata; missing: "
                + ", ".join(missing)
            )
        return
    conflicts = [str(path) for path in expected_nonplots if path.exists()]
    conflicts.extend(str(path) for path in plot_files)
    if conflicts:
        raise FileExistsError(
            "Module 2 outputs already exist. Use --overwrite for a new run or "
            f"--resume to reuse a compatible transformation: {conflicts}"
        )


def atomic_pickle_dump(value: Any, path: Path) -> None:
    """Atomically persist a pickle payload."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("wb") as handle:
        pickle.dump(value, handle, protocol=pickle.HIGHEST_PROTOCOL)
    temporary.replace(path)


def atomic_json_dump(value: dict[str, Any], path: Path) -> None:
    """Atomically persist a JSON payload."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def build_first_pass_coordinates(
    aligned_slices: Sequence[ad.AnnData], args: argparse.Namespace
) -> pd.DataFrame:
    """Build a cell-ID-keyed table of first-pass aligned coordinates."""
    if len(aligned_slices) != len(args.sample_order):
        raise ValueError(
            "The number of aligned slices does not match --sample-order"
        )
    records: list[pd.DataFrame] = []
    for expected_sample, aligned_slice in zip(args.sample_order, aligned_slices):
        if args.key_added not in aligned_slice.obsm:
            raise KeyError(
                f"Aligned slice {expected_sample!r} is missing "
                f"obsm[{args.key_added!r}]"
            )
        coordinates = np.asarray(aligned_slice.obsm[args.key_added], dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[0] != aligned_slice.n_obs:
            raise ValueError(
                f"obsm[{args.key_added!r}] for {expected_sample!r} has "
                f"unexpected shape {coordinates.shape}"
            )
        if coordinates.shape[1] < 2 or not np.isfinite(coordinates[:, :2]).all():
            raise ValueError(
                f"First-pass coordinates for {expected_sample!r} must contain "
                "at least two finite columns"
            )
        if args.adata_cell_id_column not in aligned_slice.obs:
            raise KeyError(
                f"Aligned slice {expected_sample!r} is missing obs column "
                f"{args.adata_cell_id_column!r}"
            )
        cell_ids = aligned_slice.obs[args.adata_cell_id_column].astype("string")
        if cell_ids.isna().any():
            raise ValueError(
                f"Aligned slice {expected_sample!r} contains missing cell IDs"
            )
        observed_samples = aligned_slice.obs[args.sample_column].astype("string")
        if not observed_samples.eq(expected_sample).fillna(False).all():
            raise ValueError(
                f"Aligned slice order/sample mismatch at {expected_sample!r}"
            )
        records.append(
            pd.DataFrame(
                {
                    "donor": str(args.donor),
                    "sample": expected_sample,
                    "cell_id": cell_ids.astype(str).to_numpy(),
                    "x_first_pass": coordinates[:, 0],
                    "y_first_pass": coordinates[:, 1],
                }
            )
        )
    result = pd.concat(records, ignore_index=True)
    if result["cell_id"].duplicated().any():
        examples = result.loc[
            result["cell_id"].duplicated(keep=False), "cell_id"
        ].head(10).tolist()
        raise ValueError(
            "First-pass coordinate output contains duplicate cell IDs; "
            f"examples: {examples}"
        )
    return result


def atomic_csv_dump(frame: pd.DataFrame, path: Path) -> None:
    """Atomically persist a CSV table."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def slice_names_from_models(
    models: Sequence[ad.AnnData], sample_column: str
) -> list[str]:
    """Extract the single sample name represented by each slice model."""
    names: list[str] = []
    for model in models:
        values = model.obs[sample_column].dropna().astype(str).unique()
        if len(values) != 1:
            raise ValueError(
                f"Each plotted slice must contain one {sample_column!r} value"
            )
        names.append(str(values[0]))
    return names


def display_sample_name(sample: str, donor: str) -> str:
    """Remove an exact donor prefix from a sample name for plot display only."""
    prefix = f"{donor}_"
    return str(sample)[len(prefix):] if str(sample).startswith(prefix) else str(sample)


def display_slice_names(
    models: Sequence[ad.AnnData], sample_column: str, donor: str
) -> list[str]:
    """Return slice names shortened only for subplot display."""
    return [
        display_sample_name(sample, donor)
        for sample in slice_names_from_models(models, sample_column)
    ]


def subplot_grid(n_panels: int) -> tuple[Any, np.ndarray]:
    """Create the same compact, at-most-three-column layout as Module 00."""
    if n_panels < 1:
        raise ValueError("At least one panel is required")
    n_columns = min(3, n_panels)
    n_rows = math.ceil(n_panels / n_columns)
    return plt.subplots(
        n_rows, n_columns, figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )


def plot_slices_by_celltype(
    slices: Sequence[ad.AnnData],
    args: argparse.Namespace,
    palette: dict[str, Any],
    output_path: Path,
) -> None:
    """Plot each pre-alignment slice using the common Module 00 panel style."""
    figure, axes = subplot_grid(len(slices))
    axes_flat = axes.ravel()
    display_names = display_slice_names(slices, args.sample_column, args.donor)
    for axis, ad_slice, display_name in zip(axes_flat, slices, display_names):
        coordinates = np.asarray(ad_slice.obsm[args.spatial_key], dtype=float)[:, :2]
        labels = ad_slice.obs[args.celltype_column].astype("string")
        for cell_type in palette:
            mask = labels.eq(cell_type).fillna(False).to_numpy(dtype=bool)
            if mask.any():
                axis.scatter(
                    coordinates[mask, 0], coordinates[mask, 1], s=0.2,
                    alpha=0.65, color=palette[cell_type], linewidths=0,
                )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(f"Before Spateo: {display_name}", fontsize=9)
    for axis in axes_flat[len(slices):]:
        axis.axis("off")
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markersize=4,
            markerfacecolor=color, markeredgewidth=0, label=cell_type,
        )
        for cell_type, color in palette.items()
    ]
    if handles:
        figure.legend(
            handles=handles, loc="upper right", ncol=min(5, len(handles)),
            fontsize=7, frameon=False, bbox_to_anchor=(0.99, 0.99),
        )
    figure.suptitle("Slices before initial Spateo alignment", y=1.002)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_consecutive_overlays(
    slices: Sequence[ad.AnnData],
    args: argparse.Namespace,
    spatial_key: str,
    output_path: Path,
    title_prefix: str,
    cell_type: str | None = None,
) -> None:
    """Plot ordered slice pairs using the same blue/orange style as Module 00."""
    n_pairs = len(slices) - 1
    figure, axes = subplot_grid(n_pairs)
    axes_flat = axes.ravel()
    display_names = display_slice_names(slices, args.sample_column, args.donor)
    point_size = 0.6 if cell_type is not None else 0.2
    point_alpha = 0.8 if cell_type is not None else 0.45
    for pair_index, axis in enumerate(axes_flat[:n_pairs]):
        reference = slices[pair_index]
        moving = slices[pair_index + 1]
        reference_coordinates = np.asarray(reference.obsm[spatial_key], dtype=float)[:, :2]
        moving_coordinates = np.asarray(moving.obsm[spatial_key], dtype=float)[:, :2]
        if cell_type is not None:
            reference_mask = (
                reference.obs[args.celltype_column].astype("string").eq(cell_type)
                .fillna(False).to_numpy(dtype=bool)
            )
            moving_mask = (
                moving.obs[args.celltype_column].astype("string").eq(cell_type)
                .fillna(False).to_numpy(dtype=bool)
            )
            reference_coordinates = reference_coordinates[reference_mask]
            moving_coordinates = moving_coordinates[moving_mask]
        reference_name = display_names[pair_index]
        moving_name = display_names[pair_index + 1]
        axis.scatter(
            reference_coordinates[:, 0], reference_coordinates[:, 1],
            s=point_size, alpha=point_alpha, color="#4C78A8", linewidths=0,
            label=reference_name,
        )
        axis.scatter(
            moving_coordinates[:, 0], moving_coordinates[:, 1],
            s=point_size, alpha=point_alpha, color="#F58518", linewidths=0,
            label=moving_name,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(f"{title_prefix}{reference_name} → {moving_name}", fontsize=9)
        axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    for axis in axes_flat[n_pairs:]:
        axis.axis("off")
    figure.suptitle(
        "Initial Spateo consecutive-slice overlays"
        if cell_type is None else f"Initial Spateo overlays: {cell_type}",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_initial_results(
    slices: list[ad.AnnData],
    aligned_slices: list[ad.AnnData],
    args: argparse.Namespace,
) -> None:
    """Write all Module 2 plots exclusively beneath ``args.plot_dir``."""
    celltypes = sorted(
        set().union(
            *[
                set(ad_slice.obs[args.celltype_column].dropna().astype(str))
                for ad_slice in aligned_slices
            ]
        ),
        key=str.lower,
    )
    palette = build_palette(celltypes)
    donor_safe = safe_filename(args.donor)
    plt.ioff()

    plot_slices_by_celltype(
        slices, args, palette,
        args.plot_dir / f"{donor_safe}_slices_2d_before_spateo.png",
    )
    plot_consecutive_overlays(
        aligned_slices, args, args.key_added,
        args.plot_dir / f"{donor_safe}_overlay_slices_2d_initial_spateo.png",
        title_prefix="Initial aligned: ",
    )

    celltype_plot_dir = args.plot_dir / "by_celltype"
    celltype_plot_dir.mkdir(parents=True, exist_ok=True)
    for cell_type in celltypes:
        plot_consecutive_overlays(
            aligned_slices, args, args.key_added,
            celltype_plot_dir
            / f"{safe_filename(cell_type)}_overlay_slices_2d.png",
            title_prefix=f"{cell_type}: ", cell_type=cell_type,
        )


def write_nonplot_outputs(
    args: argparse.Namespace,
    annotation_qc: pd.DataFrame,
    slice_qc: pd.DataFrame,
    metadata: dict[str, Any],
    runtime_rows: list[dict[str, Any]],
    software_versions: dict[str, str],
) -> None:
    """Write every non-plot deliverable beneath ``args.output_dir``."""
    annotation_qc.to_csv(args.output_dir / "celltype_annotation_qc.csv", index=False)
    slice_qc.to_csv(args.output_dir / "slice_qc.csv", index=False)
    pd.DataFrame.from_records(runtime_rows).to_csv(
        args.output_dir / "runtime_summary.csv", index=False
    )
    session_lines = [f"{key}: {value}" for key, value in software_versions.items()]
    (args.output_dir / "session_info.txt").write_text("\n".join(session_lines) + "\n")
    atomic_json_dump(
        {
            "alignment_metadata": metadata,
            "output_dir": str(args.output_dir),
            "plot_dir": str(args.plot_dir),
            "first_pass_coordinates_csv": str(
                first_pass_coordinates_path(args.output_dir, args.donor)
            ),
            "software_versions": software_versions,
            "status": "complete",
        },
        metadata_path(args.output_dir),
    )


def run_self_tests() -> None:
    """Run lightweight tests without importing Spateo or requiring a GPU."""
    assert display_sample_name("Br6660_NAc1_580", "Br6660") == "NAc1_580"
    assert display_sample_name("NAc1_580", "Br6660") == "NAc1_580"
    defaults = parse_args(["--run-self-tests"])
    assert defaults.device == "cuda"
    assert defaults.key_added == "align_spatial"
    assert defaults.sparse_calculation_mode and defaults.use_chunk
    try:
        validate_device("cuda", False)
    except RuntimeError:
        pass
    else:
        raise AssertionError("Unavailable CUDA device was not rejected")
    validate_device("cpu", False)
    assert transformation_path(Path("results"), "D 1").name == "Spateo_transformation_D_1.pkl"
    assert (
        first_pass_coordinates_path(Path("results"), "D 1").name
        == "D_1_first_pass_coordinates.csv"
    )

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary = Path(temporary_directory)
        annotation_path = temporary / "celltypes.csv"
        pd.DataFrame(
            {"cell_id": ["a", "b", "c"], "CellType": ["A", "B", "C"]}
        ).to_csv(annotation_path, index=False)
        annotations = pd.read_csv(annotation_path, dtype={"cell_id": "string"})
        test_adata = ad.AnnData(
            X=np.arange(6, dtype=float).reshape(3, 2),
            obs=pd.DataFrame(
                {
                    "Donor": ["D1", "D1", "D2"],
                    "Sample": ["S_100", "S_200", "S_300"],
                    "cell_id": ["a", "b", "c"],
                },
                index=["obs_a", "obs_b", "obs_c"],
            ),
            obsm={"spatial": np.array([[1.0, 10.0], [2.0, 20.0], [3.0, 30.0]])},
        )
        args = argparse.Namespace(
            donor="D1", donor_column="Donor", sample_column="Sample",
            sample_order=["S_200", "S_100"], adata_cell_id_column="cell_id",
            celltype_csv_id_column="cell_id", celltype_csv_column="CellType",
            celltype_column="CellType", spatial_key="spatial",
            key_added="align_spatial", flip_y=True,
        )
        slices, annotation_qc, slice_qc = prepare_donor_slices(
            test_adata, annotations, args
        )
        assert [value.obs["Sample"].iloc[0] for value in slices] == ["S_200", "S_100"]
        assert [value.obs["CellType"].astype(str).iloc[0] for value in slices] == ["B", "A"]
        assert annotation_qc.loc[0, "n_cells_missing_celltype"] == 0
        assert slice_qc["y_axis_flipped"].all()
        assert slice_names_from_models(slices, "Sample") == ["S_200", "S_100"]
        figure, axes = subplot_grid(2)
        assert axes.size == 2
        plt.close(figure)
        aligned_slices = [value.copy() for value in slices]
        for index, value in enumerate(aligned_slices):
            value.obsm["align_spatial"] = (
                np.asarray(value.obsm["spatial"])[:, :2] + index
            )
        coordinate_table = build_first_pass_coordinates(aligned_slices, args)
        assert coordinate_table.columns.tolist() == [
            "donor", "sample", "cell_id", "x_first_pass", "y_first_pass"
        ]
        assert coordinate_table["cell_id"].tolist() == ["b", "a"]
        assert coordinate_table["sample"].tolist() == ["S_200", "S_100"]
        coordinate_test_path = temporary / "coordinates.csv"
        atomic_csv_dump(coordinate_table, coordinate_test_path)
        roundtrip = pd.read_csv(coordinate_test_path, dtype={"cell_id": "string"})
        assert roundtrip["cell_id"].tolist() == ["b", "a"]

        output_dir = temporary / "results"
        plot_dir = temporary / "plots"
        output_dir.mkdir()
        plot_dir.mkdir()
        stale_transformation = transformation_path(output_dir, "D1")
        stale_transformation.write_text("stale")
        stale_coordinates = first_pass_coordinates_path(output_dir, "D1")
        stale_coordinates.write_text("stale")
        stale_plot = plot_dir / "stale.png"
        stale_plot.write_text("stale")
        output_args = argparse.Namespace(
            output_dir=output_dir, plot_dir=plot_dir, donor="D1",
            overwrite=True, resume=False,
        )
        prepare_output_dirs(output_args)
        assert not stale_transformation.exists()
        assert not stale_coordinates.exists()
        assert not stale_plot.exists()
        args.plot_dir = plot_dir
        plot_initial_results(slices, aligned_slices, args)
        assert (plot_dir / "D1_slices_2d_before_spateo.png").is_file()
        assert (plot_dir / "D1_overlay_slices_2d_initial_spateo.png").is_file()
        assert (plot_dir / "by_celltype" / "A_overlay_slices_2d.png").is_file()
        assert (plot_dir / "by_celltype" / "B_overlay_slices_2d.png").is_file()

        assert metadata_mismatches(
            {"key_added": "old"}, {"key_added": "align_spatial"}
        )
    LOGGER.info("All self-tests passed")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the Module 2 initial-pass alignment workflow."""
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
    )
    if args.run_self_tests:
        run_self_tests()
        return 0

    if not args.input_h5ad.is_file():
        raise FileNotFoundError(f"Input AnnData file does not exist: {args.input_h5ad}")
    if not args.celltype_csv.is_file():
        raise FileNotFoundError(f"Cell-type CSV does not exist: {args.celltype_csv}")

    prepare_output_dirs(args)
    expected_metadata = alignment_metadata(args)
    transform_file = transformation_path(args.output_dir, args.donor)
    metadata_file = metadata_path(args.output_dir)

    transformation = None
    if args.resume:
        saved_payload = json.loads(metadata_file.read_text())
        mismatches = metadata_mismatches(
            saved_payload.get("alignment_metadata", {}), expected_metadata
        )
        if mismatches:
            raise ValueError(
                "Cannot resume because transformation metadata differs:\n"
                + "\n".join(mismatches)
            )
        with transform_file.open("rb") as handle:
            transformation = pickle.load(handle)
        LOGGER.info("Reusing compatible transformation: %s", transform_file)

    st, torch, spateo_version = initialize_spateo_runtime()
    validate_device(args.device, bool(torch.cuda.is_available()))
    LOGGER.info("Spateo version: %s", spateo_version)
    LOGGER.info("Alignment device: %s", args.device)

    runtime_rows: list[dict[str, Any]] = []
    total_started = time.perf_counter()
    stage_started = time.perf_counter()
    adata = ad.read_h5ad(args.input_h5ad)
    annotations = pd.read_csv(
        args.celltype_csv, dtype={args.celltype_csv_id_column: "string"}
    )
    slices, annotation_qc, slice_qc = prepare_donor_slices(
        adata, annotations, args
    )
    del adata
    runtime_rows.append(
        {
            "stage": "load_and_prepare_slices",
            "elapsed_seconds": time.perf_counter() - stage_started,
            "n_slices": len(slices),
            "n_cells": sum(value.n_obs for value in slices),
        }
    )

    if transformation is None:
        stage_started = time.perf_counter()
        transformation = st.align.morpho_align_transformation(
            models=slices,
            spatial_key=args.spatial_key,
            key_added=args.key_added,
            device=args.device,
            max_iter=args.max_iter,
            dtype=args.dtype,
            sparse_calculation_mode=args.sparse_calculation_mode,
            use_chunk=args.use_chunk,
            chunk_capacity=args.chunk_capacity,
            verbose=args.verbose,
            partial_robust_level=args.partial_robust_level,
        )
        atomic_pickle_dump(transformation, transform_file)
        atomic_json_dump(
            {
                "alignment_metadata": expected_metadata,
                "output_dir": str(args.output_dir),
                "plot_dir": str(args.plot_dir),
                "status": "transformation_complete",
            },
            metadata_file,
        )
        runtime_rows.append(
            {
                "stage": "morpho_align_transformation",
                "elapsed_seconds": time.perf_counter() - stage_started,
                "n_slices": len(slices),
                "n_cells": sum(value.n_obs for value in slices),
            }
        )

    stage_started = time.perf_counter()
    aligned_slices = st.align.morpho_align_apply_transformation(
        models=slices,
        spatial_key=args.spatial_key,
        key_added=args.key_added,
        transformation=transformation,
        verbose=args.verbose,
    )
    runtime_rows.append(
        {
            "stage": "apply_transformation",
            "elapsed_seconds": time.perf_counter() - stage_started,
            "n_slices": len(aligned_slices),
            "n_cells": sum(value.n_obs for value in aligned_slices),
        }
    )

    stage_started = time.perf_counter()
    first_pass_coordinates = build_first_pass_coordinates(aligned_slices, args)
    coordinates_file = first_pass_coordinates_path(args.output_dir, args.donor)
    atomic_csv_dump(first_pass_coordinates, coordinates_file)
    runtime_rows.append(
        {
            "stage": "write_first_pass_coordinates",
            "elapsed_seconds": time.perf_counter() - stage_started,
            "n_slices": len(aligned_slices),
            "n_cells": len(first_pass_coordinates),
        }
    )
    LOGGER.info("Saved first-pass coordinates: %s", coordinates_file)

    stage_started = time.perf_counter()
    plot_initial_results(slices, aligned_slices, args)
    runtime_rows.append(
        {
            "stage": "plotting",
            "elapsed_seconds": time.perf_counter() - stage_started,
            "n_slices": len(aligned_slices),
            "n_cells": sum(value.n_obs for value in aligned_slices),
        }
    )
    runtime_rows.append(
        {
            "stage": "total",
            "elapsed_seconds": time.perf_counter() - total_started,
            "n_slices": len(aligned_slices),
            "n_cells": sum(value.n_obs for value in aligned_slices),
        }
    )

    software_versions = {
        "python": platform.python_version(),
        "anndata": ad.__version__,
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "matplotlib": matplotlib.__version__,
        "spateo": spateo_version,
        "torch": torch.__version__,
    }
    write_nonplot_outputs(
        args, annotation_qc, slice_qc, expected_metadata, runtime_rows,
        software_versions,
    )
    LOGGER.info(
        "Module 2 complete. Non-plot outputs: %s; plots: %s",
        args.output_dir, args.plot_dir,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())


# Example:
# python 02_Spateo_initial_pass.py \
#     --input-h5ad /path/to/spe_NormCounts_nucleus_normcounts.h5ad \
#     --celltype-csv /path/to/Banksy_cell_types.csv \
#     --output-dir /path/to/module02_initial_pass_Br6660 \
#     --plot-dir /path/to/module02_initial_pass_Br6660_plots \
#     --donor Br6660 \
#     --donor-column Donor \
#     --sample-column Sample \
#     --sample-order \
#         Br6660_NAc1_580 \
#         Br6660_NAc2_1090 \
#         Br6660_NAc3_1580 \
#         Br6660_NAc4_2080 \
#         Br6660_NAc5_2580 \
#         Br6660_NAc6_3080 \
#         Br6660_NAc7_3580 \
#         Br6660_Nac10_4080 \
#         Br6660_NAc8_4580 \
#         Br6660_NAc9_5080 \
#         Br6660_Nac11_5580 \
#     --adata-cell-id-column cell_id \
#     --celltype-csv-id-column cell_id \
#     --celltype-csv-column CellType \
#     --celltype-column CellType \
#     --spatial-key spatial \
#     --key-added align_spatial \
#     --device cuda \
#     --max-iter 200 \
#     --partial-robust-level 30 \
#     --chunk-capacity 2 \
#     --sparse-calculation-mode \
#     --use-chunk \
#     --flip-y \
#     --overwrite
