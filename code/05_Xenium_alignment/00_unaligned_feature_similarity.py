#!/usr/bin/env python3
"""Evaluate unaligned consecutive Xenium slices using all genes in ``adata.X``.

Coordinates are converted explicitly to micrometers and optionally flipped
along Y within each slice to match the orientation used by the Spateo pipeline.
The script evaluates ordered consecutive pairs without modifying or aligning
their coordinates. Quantitative results and plots are written separately.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import math
import os
import platform
import sys
import tempfile
import time
from pathlib import Path
from typing import Any, Sequence

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from feature_similarity import (
    DEFAULT_GRID_SIZES_UM,
    SUPPORTED_METRICS,
    compute_pair_gene_similarity,
    summarize_pair_gene_similarity,
)


LOGGER = logging.getLogger("module00_unaligned_feature_similarity")
SCRIPT_VERSION = "1.0"
PLOT_METRICS = (
    ("pcc_all", "Pearson correlation"),
    ("cos_sim_all", "Cosine similarity"),
    ("ssim_all", "SSIM"),
    ("mi_all", "Mutual information"),
)
PAIR_SUMMARY_METRICS = (
    "pcc_all",
    "cos_sim_all",
    "ssim_all",
    "mi_all",
)


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
    parser.add_argument("--coordinate-unit", choices=("um", "mm"))
    parser.add_argument(
        "--flip-y", action=argparse.BooleanOptionalAction, default=True,
        help="Reflect Y independently within each slice before evaluation.",
    )
    parser.add_argument(
        "--grid-sizes-um", type=float, nargs="+",
        default=list(DEFAULT_GRID_SIZES_UM),
    )
    parser.add_argument("--min-cells-reference", type=int, default=5)
    parser.add_argument("--min-cells-moving", type=int, default=5)
    parser.add_argument("--min-overlap-grids", type=int, default=10)
    parser.add_argument("--min-gene-nz-grids", type=int, default=5)
    parser.add_argument("--min-gene-valid-offsets", type=int, default=1)
    parser.add_argument("--min-gene-valid-scales", type=int, default=1)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument(
        "--n-jobs", type=int, default=4,
        help="Consecutive pairs evaluated concurrently; use -1 for all CPUs.",
    )
    parser.add_argument("--histogram-bins", type=int, default=30)
    parser.add_argument("--plot-dpi", type=int, default=250)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--run-self-tests", action="store_true")
    args = parser.parse_args(argv)

    required = (
        "input_h5ad", "celltype_csv", "output_dir", "plot_dir", "donor",
        "sample_order", "coordinate_unit",
    )
    if not args.run_self_tests and any(getattr(args, name) is None for name in required):
        parser.error(
            "--input-h5ad, --celltype-csv, --output-dir, --plot-dir, --donor, "
            "--sample-order, and --coordinate-unit are required unless "
            "--run-self-tests is used"
        )
    if args.sample_order and len(args.sample_order) < 2:
        parser.error("--sample-order requires at least two ordered slices")
    if args.sample_order and len(args.sample_order) != len(set(args.sample_order)):
        parser.error("--sample-order contains duplicate names")
    if args.grid_sizes_um and (
        any(not np.isfinite(value) or value <= 0 for value in args.grid_sizes_um)
        or len(args.grid_sizes_um) != len(set(args.grid_sizes_um))
    ):
        parser.error("--grid-sizes-um must contain unique positive finite values")
    for name in (
        "min_cells_reference", "min_cells_moving", "min_overlap_grids",
        "min_gene_valid_offsets", "min_gene_valid_scales", "histogram_bins",
        "plot_dpi",
    ):
        if getattr(args, name) < 1:
            parser.error(f"--{name.replace('_', '-')} must be positive")
    if args.min_gene_nz_grids < 0:
        parser.error("--min-gene-nz-grids must be nonnegative")
    if args.n_jobs == 0 or args.n_jobs < -1:
        parser.error("--n-jobs must be -1 or a positive integer")
    if (
        not args.run_self_tests
        and args.output_dir.resolve() == args.plot_dir.resolve()
    ):
        parser.error("--output-dir and --plot-dir must be different directories")
    return args


def safe_filename(value: str) -> str:
    """Convert a value into a conservative filename component."""
    import re

    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "unnamed"


def display_sample_name(sample: str, donor: str) -> str:
    """Remove an exact donor prefix from a sample name for plot display only."""
    prefix = f"{donor}_"
    return str(sample)[len(prefix):] if str(sample).startswith(prefix) else str(sample)


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    """Return every generated output path."""
    donor = safe_filename(args.donor)
    return {
        "multiscale": args.output_dir / f"{donor}_unaligned_feature_similarity_multiscale_by_gene.csv",
        "scale": args.output_dir / f"{donor}_unaligned_feature_similarity_by_scale_gene.csv",
        "offset": args.output_dir / f"{donor}_unaligned_feature_similarity_by_offset_gene.csv",
        "pair_summary": args.output_dir / f"{donor}_unaligned_feature_similarity_pair_summary.csv",
        "metadata": args.output_dir / f"{donor}_unaligned_feature_similarity_metadata.json",
        "runtime": args.output_dir / f"{donor}_unaligned_feature_similarity_runtime.csv",
        "plot": args.plot_dir / f"{donor}_unaligned_feature_similarity_distributions.png",
        "overlay_plot": args.plot_dir / f"{donor}_overlay_slices_2d_unaligned.png",
    }


def prepare_output_dirs(args: argparse.Namespace) -> dict[str, Path]:
    """Create output directories and protect existing deliverables."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(args)
    conflicts = [path for path in paths.values() if path.exists()]
    celltype_plot_dir = args.plot_dir / "by_celltype"
    celltype_conflicts = (
        list(celltype_plot_dir.glob("*.png"))
        if celltype_plot_dir.exists() else []
    )
    conflicts.extend(celltype_conflicts)
    if conflicts and not args.overwrite:
        raise FileExistsError(
            "Outputs already exist; use --overwrite: "
            + ", ".join(map(str, conflicts))
        )
    if args.overwrite:
        for path in conflicts:
            if path.is_file():
                path.unlink()
    return paths


def atomic_csv_dump(frame: pd.DataFrame, path: Path) -> None:
    """Atomically write a CSV table."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def atomic_json_dump(payload: dict[str, Any], path: Path) -> None:
    """Atomically write JSON metadata."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def prepare_ordered_slices(
    adata: ad.AnnData,
    annotations: pd.DataFrame,
    args: argparse.Namespace,
) -> tuple[list[ad.AnnData], pd.DataFrame]:
    """Subset one donor and return validated slices in anatomical order."""
    required_obs = [
        args.donor_column, args.sample_column, args.adata_cell_id_column,
    ]
    missing_obs = [column for column in required_obs if column not in adata.obs]
    if missing_obs:
        raise ValueError(f"Missing required adata.obs columns: {missing_obs}")
    required_annotation_columns = [
        args.celltype_csv_id_column, args.celltype_csv_column,
    ]
    missing_annotation_columns = [
        column for column in required_annotation_columns
        if column not in annotations.columns
    ]
    if missing_annotation_columns:
        raise ValueError(
            "Missing required cell-type CSV columns: "
            f"{missing_annotation_columns}"
        )
    if args.spatial_key not in adata.obsm:
        raise ValueError(f"Missing adata.obsm[{args.spatial_key!r}]")
    donor_mask = (
        adata.obs[args.donor_column].astype("string").eq(str(args.donor))
        .fillna(False).to_numpy(dtype=bool)
    )
    if not donor_mask.any():
        observed = sorted(adata.obs[args.donor_column].dropna().astype(str).unique())
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
            f"CSV column {args.celltype_csv_id_column!r} must be nonmissing "
            "and unique"
        )
    annotation_values = annotations[args.celltype_csv_column].astype("string")
    celltype_mapping = pd.Series(
        annotation_values.to_numpy(), index=annotation_ids.to_numpy()
    )
    mapped_celltypes = donor_ids.map(celltype_mapping)
    if mapped_celltypes.isna().any():
        examples = donor_ids[mapped_celltypes.isna()].astype(str).head(10).tolist()
        raise ValueError(
            "Cell-type annotations are missing for selected donor cells; "
            f"example IDs: {examples}"
        )
    donor_adata.obs[args.celltype_column] = pd.Categorical(
        mapped_celltypes.astype(str)
    )
    samples = donor_adata.obs[args.sample_column].astype("string")
    observed_samples = set(samples.dropna().astype(str))
    absent = [sample for sample in args.sample_order if sample not in observed_samples]
    if absent:
        raise ValueError(f"Requested samples are absent from donor {args.donor}: {absent}")

    conversion = 1000.0 if args.coordinate_unit == "mm" else 1.0
    slices: list[ad.AnnData] = []
    qc_rows: list[dict[str, Any]] = []
    for sample in args.sample_order:
        mask = samples.eq(sample).fillna(False).to_numpy(dtype=bool)
        ad_slice = donor_adata[mask].copy()
        coordinates = np.asarray(ad_slice.obsm[args.spatial_key], dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[0] != ad_slice.n_obs or coordinates.shape[1] < 2:
            raise ValueError(
                f"Coordinates for {sample!r} have unexpected shape {coordinates.shape}"
            )
        coordinates = coordinates[:, :2].copy() * conversion
        if not np.isfinite(coordinates).all():
            raise ValueError(f"Coordinates for {sample!r} contain nonfinite values")
        midpoint = np.nan
        if args.flip_y:
            midpoint = float(coordinates[:, 1].min() + coordinates[:, 1].max())
            coordinates[:, 1] = midpoint - coordinates[:, 1]
        ad_slice.obsm[args.spatial_key] = coordinates
        slices.append(ad_slice)
        qc_rows.append(
            {
                "slice": sample,
                "n_cells": ad_slice.n_obs,
                "n_genes": ad_slice.n_vars,
                "x_min_um": float(coordinates[:, 0].min()),
                "x_max_um": float(coordinates[:, 0].max()),
                "y_min_um": float(coordinates[:, 1].min()),
                "y_max_um": float(coordinates[:, 1].max()),
                "y_axis_flipped": bool(args.flip_y),
                "y_flip_reflection_sum_um": midpoint,
            }
        )
    return slices, pd.DataFrame.from_records(qc_rows)


def _pair_columns(
    frame: pd.DataFrame, pair_number: int, reference: str, moving: str
) -> pd.DataFrame:
    """Prepend stable ordered-pair identifiers to a result table."""
    result = frame.copy()
    result.insert(0, "moving_slice", moving)
    result.insert(0, "reference_slice", reference)
    result.insert(0, "pair_id", f"pair_{pair_number:02d}_{reference}__{moving}")
    result.insert(0, "pair_number", pair_number)
    return result


def evaluate_consecutive_pairs(
    slices: Sequence[ad.AnnData], args: argparse.Namespace
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    """Evaluate every ordered anatomical neighbor pair, optionally in threads."""
    pair_inputs = [
        (pair_number, reference_slice, moving_slice)
        for pair_number, (reference_slice, moving_slice) in enumerate(
            zip(slices[:-1], slices[1:]), 1
        )
    ]

    def evaluate_one_pair(
        pair_number: int,
        reference_slice: ad.AnnData,
        moving_slice: ad.AnnData,
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
        reference_name = args.sample_order[pair_number - 1]
        moving_name = args.sample_order[pair_number]
        LOGGER.info(
            "Evaluating unaligned pair %d: %s -> %s",
            pair_number, reference_name, moving_name,
        )
        started = time.perf_counter()
        multiscale, scale, offset = compute_pair_gene_similarity(
            reference_slice,
            moving_slice,
            reference_coord_key=args.spatial_key,
            moving_coord_key=args.spatial_key,
            gene_list=None,
            expression_layer=None,
            grid_sizes_um=args.grid_sizes_um,
            grid_offsets_um=None,
            min_cells_reference=args.min_cells_reference,
            min_cells_moving=args.min_cells_moving,
            min_overlap_grids=args.min_overlap_grids,
            metrics=SUPPORTED_METRICS,
            random_state=args.random_state,
        )
        metric_summaries = {
            primary_metric: summarize_pair_gene_similarity(
                scale,
                primary_metric=primary_metric,
                min_gene_nz_grids=args.min_gene_nz_grids,
                min_gene_valid_offsets=args.min_gene_valid_offsets,
                min_gene_valid_scales=args.min_gene_valid_scales,
            )
            for primary_metric in PAIR_SUMMARY_METRICS
        }
        pair_id = f"pair_{pair_number:02d}_{reference_name}__{moving_name}"
        summary_row: dict[str, Any] = {
            "pair_number": pair_number,
            "pair_id": pair_id,
            "reference_slice": reference_name,
            "moving_slice": moving_name,
        }
        for metric, metric_summary in metric_summaries.items():
            summary_row[f"{metric}_multiscale_median"] = (
                metric_summary["multiscale_median_primary_score"]
            )
            summary_row[f"{metric}_n_genes_finite"] = (
                metric_summary["n_genes_primary_metric_finite"]
            )
            summary_row[f"{metric}_n_genes_eligible"] = (
                metric_summary["n_genes_primary_metric_eligible"]
            )
            summary_row[f"{metric}_finite_gene_fraction"] = (
                metric_summary["finite_gene_fraction"]
            )
            summary_row[f"{metric}_eligible_gene_fraction"] = (
                metric_summary["eligible_gene_fraction"]
            )
            summary_row[f"{metric}_per_scale_median"] = json.dumps(
                metric_summary["per_scale_median_primary_score"], sort_keys=True
            )
            summary_row[f"{metric}_per_scale_n_genes_eligible"] = json.dumps(
                metric_summary["per_scale_n_genes_eligible"], sort_keys=True
            )

            values = multiscale[metric].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if values.size:
                q25 = float(np.quantile(values, 0.25))
                median = float(np.median(values))
                q75 = float(np.quantile(values, 0.75))
                summary_row[f"{metric}_gene_median"] = median
                summary_row[f"{metric}_gene_q25"] = q25
                summary_row[f"{metric}_gene_q75"] = q75
                summary_row[f"{metric}_gene_iqr"] = q75 - q25
                summary_row[f"{metric}_n_genes_distribution"] = int(values.size)
            else:
                summary_row[f"{metric}_gene_median"] = np.nan
                summary_row[f"{metric}_gene_q25"] = np.nan
                summary_row[f"{metric}_gene_q75"] = np.nan
                summary_row[f"{metric}_gene_iqr"] = np.nan
                summary_row[f"{metric}_n_genes_distribution"] = 0
        return (
            _pair_columns(multiscale, pair_number, reference_name, moving_name),
            _pair_columns(scale, pair_number, reference_name, moving_name),
            _pair_columns(offset, pair_number, reference_name, moving_name),
            summary_row,
            {
                "pair_number": pair_number,
                "reference_slice": reference_name,
                "moving_slice": moving_name,
                "elapsed_seconds": time.perf_counter() - started,
                "n_reference_cells": reference_slice.n_obs,
                "n_moving_cells": moving_slice.n_obs,
                "n_genes_scored": multiscale["gene"].nunique(),
            },
        )

    requested_workers = (os.cpu_count() or 1) if args.n_jobs == -1 else args.n_jobs
    n_workers = min(len(pair_inputs), requested_workers)
    LOGGER.info(
        "Evaluating %d consecutive pairs with %d CPU worker thread(s)",
        len(pair_inputs), n_workers,
    )
    if n_workers == 1:
        completed = [evaluate_one_pair(*pair_input) for pair_input in pair_inputs]
    else:
        # AnnData objects remain shared in memory. Limiting native BLAS pools
        # prevents each pair worker from creating another full CPU thread pool.
        with threadpool_limits(limits=1):
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                completed = list(executor.map(lambda values: evaluate_one_pair(*values), pair_inputs))
    completed.sort(key=lambda value: int(value[3]["pair_number"]))
    multiscale_tables = [value[0] for value in completed]
    scale_tables = [value[1] for value in completed]
    offset_tables = [value[2] for value in completed]
    summary_rows = [value[3] for value in completed]
    runtime_rows = [value[4] for value in completed]
    return (
        pd.concat(multiscale_tables, ignore_index=True),
        pd.concat(scale_tables, ignore_index=True),
        pd.concat(offset_tables, ignore_index=True),
        pd.DataFrame.from_records(summary_rows),
        runtime_rows,
    )


def plot_metric_distributions(
    multiscale_result: pd.DataFrame,
    output_path: Path,
    dpi: int,
    bins: int,
    donor: str,
) -> None:
    """Plot all-gene multiscale metric distributions as pair rows and metric columns."""
    pairs = (
        multiscale_result[
            ["pair_number", "pair_id", "reference_slice", "moving_slice"]
        ]
        .drop_duplicates()
        .sort_values("pair_number")
    )
    figure, axes = plt.subplots(
        len(pairs), len(PLOT_METRICS),
        figsize=(4.2 * len(PLOT_METRICS), 2.7 * len(pairs)),
        squeeze=False,
    )
    for row_index, pair in enumerate(pairs.itertuples(index=False)):
        pair_data = multiscale_result[multiscale_result["pair_id"] == pair.pair_id]
        for column_index, (metric, title) in enumerate(PLOT_METRICS):
            axis = axes[row_index, column_index]
            values = pair_data[metric].to_numpy(float)
            values = values[np.isfinite(values)]
            if values.size:
                axis.hist(values, bins=bins, color="#4C78A8", alpha=0.85)
                median = float(np.median(values))
                axis.axvline(median, color="#D62728", linestyle="--", linewidth=1)
                axis.text(
                    0.98, 0.95, f"n={values.size}\nmedian={median:.3g}",
                    transform=axis.transAxes, ha="right", va="top", fontsize=8,
                )
            else:
                axis.text(0.5, 0.5, "No finite values", ha="center", va="center")
            axis.set_title(
                f"{title}\n"
                f"{display_sample_name(pair.reference_slice, donor)} → "
                f"{display_sample_name(pair.moving_slice, donor)}",
                fontsize=9,
            )
            if column_index == 0:
                axis.set_ylabel(
                    f"{display_sample_name(pair.reference_slice, donor)}\n→ "
                    f"{display_sample_name(pair.moving_slice, donor)}\nGene count"
                )
            else:
                axis.set_ylabel("Gene count")
            if row_index == len(pairs) - 1:
                axis.set_xlabel(metric)
    figure.suptitle("Unaligned multiscale feature similarity across all genes", y=1.002)
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_unaligned_pair_overlays(
    slices: Sequence[ad.AnnData],
    sample_order: Sequence[str],
    spatial_key: str,
    output_path: Path,
    dpi: int,
    donor: str,
    celltype_column: str = "CellType",
    cell_type: str | None = None,
) -> None:
    """Overlay every unaligned consecutive pair with informative slice titles."""
    n_pairs = len(slices) - 1
    if n_pairs < 1:
        raise ValueError("At least two slices are required for an overlay plot")
    n_columns = min(3, n_pairs)
    n_rows = math.ceil(n_pairs / n_columns)
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )
    axes_flat = axes.ravel()
    point_size = 0.6 if cell_type is not None else 0.2
    point_alpha = 0.8 if cell_type is not None else 0.45
    for pair_index, axis in enumerate(axes_flat[:n_pairs]):
        reference_coordinates = np.asarray(
            slices[pair_index].obsm[spatial_key], dtype=float
        )[:, :2]
        moving_coordinates = np.asarray(
            slices[pair_index + 1].obsm[spatial_key], dtype=float
        )[:, :2]
        if cell_type is not None:
            reference_mask = (
                slices[pair_index].obs[celltype_column].astype("string")
                .eq(cell_type).fillna(False).to_numpy(dtype=bool)
            )
            moving_mask = (
                slices[pair_index + 1].obs[celltype_column].astype("string")
                .eq(cell_type).fillna(False).to_numpy(dtype=bool)
            )
            reference_coordinates = reference_coordinates[reference_mask]
            moving_coordinates = moving_coordinates[moving_mask]
        reference_name = sample_order[pair_index]
        moving_name = sample_order[pair_index + 1]
        reference_display = display_sample_name(reference_name, donor)
        moving_display = display_sample_name(moving_name, donor)
        axis.scatter(
            reference_coordinates[:, 0], reference_coordinates[:, 1],
            s=point_size, alpha=point_alpha, color="#4C78A8", linewidths=0,
            label=reference_display,
        )
        axis.scatter(
            moving_coordinates[:, 0], moving_coordinates[:, 1],
            s=point_size, alpha=point_alpha, color="#F58518", linewidths=0,
            label=moving_display,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(
            (
                f"Unaligned: {reference_display} → {moving_display}"
                if cell_type is None
                else f"{cell_type}: {reference_display} → {moving_display}"
            ),
            fontsize=9,
        )
        axis.legend(
            loc="upper right", fontsize=6, markerscale=8, frameon=False,
        )
    for axis in axes_flat[n_pairs:]:
        axis.axis("off")
    figure.suptitle(
        "Unaligned consecutive-slice overlays"
        if cell_type is None else f"Unaligned overlays: {cell_type}",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def run_self_tests() -> None:
    """Run lightweight integration tests for slice preparation and pair evaluation."""
    assert display_sample_name("Br6660_NAc1_580", "Br6660") == "NAc1_580"
    assert display_sample_name("NAc1_580", "Br6660") == "NAc1_580"
    coordinates = np.array(
        [(x, y) for y in range(4) for x in range(4)], dtype=float
    )
    expression = np.column_stack(
        [coordinates[:, 0] + coordinates[:, 1] + 1.0, coordinates[:, 0] + 1.0]
    )
    full_expression = np.vstack([expression, expression, expression])
    test = ad.AnnData(
        X=full_expression,
        obs=pd.DataFrame(
            {
                "Donor": ["D1"] * 48,
                "Sample": np.repeat(["S1", "S2", "S3"], 16),
                "cell_id": [f"cell_{index}" for index in range(48)],
            }
        ),
        obsm={"spatial": np.vstack([coordinates, coordinates, coordinates])},
    )
    test.var_names = ["g1", "g2"]
    annotations = pd.DataFrame(
        {
            "cell_id": [f"cell_{index}" for index in range(48)],
            "CellType": np.where(np.arange(48) % 2 == 0, "A", "B"),
        }
    )
    args = argparse.Namespace(
        donor="D1", donor_column="Donor", sample_column="Sample",
        sample_order=["S1", "S2", "S3"], spatial_key="spatial",
        adata_cell_id_column="cell_id", celltype_csv_id_column="cell_id",
        celltype_csv_column="CellType", celltype_column="CellType",
        coordinate_unit="um", flip_y=True, grid_sizes_um=[1.0, 2.0],
        min_cells_reference=1, min_cells_moving=1, min_overlap_grids=2,
        min_gene_nz_grids=1, min_gene_valid_offsets=1,
        min_gene_valid_scales=1, random_state=0, n_jobs=2,
    )
    slices, qc = prepare_ordered_slices(test, annotations, args)
    assert len(slices) == 3 and qc["y_axis_flipped"].all()
    multiscale, scale, offset, summary, runtime = evaluate_consecutive_pairs(
        slices, args
    )
    assert multiscale["pair_id"].nunique() == 2
    assert set(metric for metric, _ in PLOT_METRICS) <= set(multiscale.columns)
    assert scale["pair_id"].nunique() == 2
    assert offset["pair_id"].nunique() == 2
    assert len(summary) == 2 and len(runtime) == 2
    assert summary["pair_id"].nunique() == 2
    for metric in PAIR_SUMMARY_METRICS:
        for suffix in (
            "multiscale_median", "n_genes_eligible",
            "eligible_gene_fraction", "per_scale_median", "gene_median",
            "gene_q25", "gene_q75", "gene_iqr",
        ):
            assert f"{metric}_{suffix}" in summary.columns
    sequential_args = argparse.Namespace(**vars(args))
    sequential_args.n_jobs = 1
    sequential_multiscale, _, _, sequential_summary, _ = evaluate_consecutive_pairs(
        slices, sequential_args
    )
    pd.testing.assert_frame_equal(multiscale, sequential_multiscale)
    pd.testing.assert_frame_equal(summary, sequential_summary)
    with tempfile.TemporaryDirectory() as temporary_directory:
        plot_path = Path(temporary_directory) / "distribution.png"
        plot_metric_distributions(
            multiscale, plot_path, dpi=80, bins=5, donor=args.donor
        )
        assert plot_path.is_file()
        overlay_path = Path(temporary_directory) / "overlay.png"
        plot_unaligned_pair_overlays(
            slices, args.sample_order, args.spatial_key, overlay_path, dpi=80,
            donor=args.donor,
        )
        assert overlay_path.is_file()
        celltype_overlay_path = Path(temporary_directory) / "celltype_overlay.png"
        plot_unaligned_pair_overlays(
            slices, args.sample_order, args.spatial_key,
            celltype_overlay_path, dpi=80, donor=args.donor,
            celltype_column=args.celltype_column, cell_type="A",
        )
        assert celltype_overlay_path.is_file()
    LOGGER.info("All self-tests passed")


def main(argv: Sequence[str] | None = None) -> int:
    """Run unaligned consecutive-pair feature-similarity evaluation."""
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
        raise FileNotFoundError(
            f"Cell-type CSV does not exist: {args.celltype_csv}"
        )
    paths = prepare_output_dirs(args)
    total_started = time.perf_counter()
    LOGGER.info("Loading AnnData: %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    annotations = pd.read_csv(
        args.celltype_csv, dtype={args.celltype_csv_id_column: "string"}
    )
    slices, slice_qc = prepare_ordered_slices(adata, annotations, args)
    del adata
    multiscale, scale, offset, pair_summary, runtime_rows = (
        evaluate_consecutive_pairs(slices, args)
    )
    atomic_csv_dump(multiscale, paths["multiscale"])
    atomic_csv_dump(scale, paths["scale"])
    atomic_csv_dump(offset, paths["offset"])
    atomic_csv_dump(pair_summary, paths["pair_summary"])
    plot_metric_distributions(
        multiscale, paths["plot"], dpi=args.plot_dpi, bins=args.histogram_bins,
        donor=args.donor,
    )
    plot_unaligned_pair_overlays(
        slices, args.sample_order, args.spatial_key,
        paths["overlay_plot"], dpi=args.plot_dpi, donor=args.donor,
    )
    celltype_plot_dir = args.plot_dir / "by_celltype"
    celltype_plot_dir.mkdir(parents=True, exist_ok=True)
    celltypes = sorted(
        set().union(
            *[
                set(ad_slice.obs[args.celltype_column].dropna().astype(str))
                for ad_slice in slices
            ]
        ),
        key=str.lower,
    )
    for cell_type in celltypes:
        plot_unaligned_pair_overlays(
            slices, args.sample_order, args.spatial_key,
            celltype_plot_dir
            / f"{safe_filename(cell_type)}_overlay_slices_2d.png",
            dpi=args.plot_dpi, donor=args.donor,
            celltype_column=args.celltype_column, cell_type=cell_type,
        )
    runtime_rows.append(
        {
            "pair_number": np.nan,
            "reference_slice": "ALL",
            "moving_slice": "ALL",
            "elapsed_seconds": time.perf_counter() - total_started,
            "n_reference_cells": np.nan,
            "n_moving_cells": np.nan,
            "n_genes_scored": multiscale["gene"].nunique(),
        }
    )
    atomic_csv_dump(pd.DataFrame.from_records(runtime_rows), paths["runtime"])
    input_stat = args.input_h5ad.stat()
    metadata = {
        "script_version": SCRIPT_VERSION,
        "input_h5ad": str(args.input_h5ad.resolve()),
        "input_h5ad_size": input_stat.st_size,
        "input_h5ad_mtime_ns": input_stat.st_mtime_ns,
        "celltype_csv": str(args.celltype_csv.resolve()),
        "celltype_csv_id_column": args.celltype_csv_id_column,
        "celltype_csv_column": args.celltype_csv_column,
        "celltype_column": args.celltype_column,
        "donor": str(args.donor),
        "donor_column": args.donor_column,
        "sample_column": args.sample_column,
        "sample_order": list(args.sample_order),
        "spatial_key": args.spatial_key,
        "coordinate_unit_input": args.coordinate_unit,
        "coordinate_unit_analysis": "um",
        "flip_y": args.flip_y,
        "expression_source": "adata.X",
        "gene_selection": "all_genes_shared_by_each_pair",
        "metrics": list(SUPPORTED_METRICS),
        "distribution_metric_columns": [metric for metric, _ in PLOT_METRICS],
        "pair_summary_metrics": list(PAIR_SUMMARY_METRICS),
        "pair_summary_aggregation": (
            "median across eligible genes within scale, followed by median "
            "across finite scales"
        ),
        "gene_distribution_statistics": ["median", "q25", "q75", "iqr"],
        "grid_sizes_um": list(map(float, args.grid_sizes_um)),
        "grid_offsets": "feature_similarity_default_four_offsets_per_scale",
        "min_cells_reference": args.min_cells_reference,
        "min_cells_moving": args.min_cells_moving,
        "min_overlap_grids": args.min_overlap_grids,
        "min_gene_nz_grids": args.min_gene_nz_grids,
        "min_gene_valid_offsets": args.min_gene_valid_offsets,
        "min_gene_valid_scales": args.min_gene_valid_scales,
        "random_state": args.random_state,
        "n_jobs_requested": args.n_jobs,
        "parallelization": "shared-memory consecutive-pair threads",
        "n_pairs": int(multiscale["pair_id"].nunique()),
        "n_genes_union_in_results": int(multiscale["gene"].nunique()),
        "slice_qc": slice_qc.to_dict(orient="records"),
        "outputs": {name: str(path) for name, path in paths.items()},
        "by_celltype_plot_dir": str(celltype_plot_dir),
        "software_versions": {
            "python": platform.python_version(),
            "anndata": ad.__version__,
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    atomic_json_dump(metadata, paths["metadata"])
    LOGGER.info("Module 00 complete. Metrics: %s; plot: %s", args.output_dir, paths["plot"])
    return 0


if __name__ == "__main__":
    sys.exit(main())


# Example usage:
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
# python 00_unaligned_feature_similarity.py \
#     --input-h5ad /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/02_build_spe/h5ad/spe_NormCounts_nucleus_normcounts.h5ad \
#     --output-dir /path/to/unaligned_feature_similarity \
#     --plot-dir /path/to/unaligned_feature_similarity_plots \
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
#     --spatial-key spatial \
#     --coordinate-unit um \
#     --grid-sizes-um 200 300 400 \
#     --min-cells-reference 5 \
#     --min-cells-moving 5 \
#     --min-overlap-grids 10 \
#     --min-gene-nz-grids 5 \
#     --min-gene-valid-offsets 1 \
#     --min-gene-valid-scales 2 \
#     --n-jobs 4 \
#     --random-state 0 \
#     --overwrite
