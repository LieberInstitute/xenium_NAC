#!/usr/bin/env python3
"""Evaluate aligned Spateo coordinates across consecutive Xenium slices.

Expression and annotations come from the original AnnData object, while the
evaluation coordinates come exclusively from either an aligned-coordinate
CSV or a saved Spateo transformation applied again to reconstructed slices.
The alignment pass controls coordinate-column selection and output naming.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import logging
import math
import os
import pickle
import platform
import re
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


LOGGER = logging.getLogger("module03_spateo_feature_similarity")
SCRIPT_VERSION = "2.0"
PLOT_METRICS = (
    ("pcc_all", "Pearson correlation"),
    ("cos_sim_all", "Cosine similarity"),
    ("ssim_all", "SSIM"),
    ("mi_all", "Mutual information"),
)
PAIR_SUMMARY_METRICS = tuple(metric for metric, _ in PLOT_METRICS)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-h5ad", type=Path)
    parser.add_argument(
        "--alignment-pass", choices=("first", "second"),
        help="Alignment pass represented by the supplied coordinates.",
    )
    coordinate_source = parser.add_mutually_exclusive_group()
    coordinate_source.add_argument("--aligned-coordinate-csv", type=Path)
    coordinate_source.add_argument("--transformation-pkl", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument("--donor")
    parser.add_argument("--donor-column", default="Donor")
    parser.add_argument("--sample-column", default="Sample")
    parser.add_argument("--sample-order", nargs="+")
    parser.add_argument("--adata-cell-id-column", default="cell_id")
    parser.add_argument("--spatial-key", default="spatial")
    parser.add_argument("--key-added", default="align_spatial")
    parser.add_argument("--coordinate-unit", choices=("um", "mm"))
    parser.add_argument(
        "--flip-y", action=argparse.BooleanOptionalAction, default=True,
        help="PKL mode only: reproduce the per-slice Y reflection used in alignment.",
    )
    parser.add_argument(
        "--verbose", action=argparse.BooleanOptionalAction, default=True,
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
    parser.add_argument("--n-jobs", type=int, default=4)
    parser.add_argument("--histogram-bins", type=int, default=30)
    parser.add_argument("--plot-dpi", type=int, default=250)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--run-self-tests", action="store_true")
    args = parser.parse_args(argv)

    required = (
        "input_h5ad", "alignment_pass", "output_dir", "plot_dir", "donor",
        "sample_order", "coordinate_unit",
    )
    if not args.run_self_tests and any(getattr(args, name) is None for name in required):
        parser.error(
            "--input-h5ad, --alignment-pass, --output-dir, --plot-dir, --donor, "
            "--sample-order, and --coordinate-unit are required unless "
            "--run-self-tests is used"
        )
    if not args.run_self_tests and not (
        args.aligned_coordinate_csv or args.transformation_pkl
    ):
        parser.error(
            "Exactly one of --aligned-coordinate-csv or --transformation-pkl "
            "is required"
        )
    if args.sample_order and len(args.sample_order) < 2:
        parser.error("--sample-order requires at least two slices")
    if args.sample_order and len(args.sample_order) != len(set(args.sample_order)):
        parser.error("--sample-order contains duplicates")
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
    if not args.run_self_tests and args.output_dir.resolve() == args.plot_dir.resolve():
        parser.error("--output-dir and --plot-dir must differ")
    return args


def safe_filename(value: str) -> str:
    """Convert a value into a conservative filename component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "unnamed"


def alignment_pass_name(args: argparse.Namespace) -> str:
    """Return the filename-safe alignment-pass name."""
    return f"{args.alignment_pass}_pass"


def alignment_pass_label(args: argparse.Namespace) -> str:
    """Return the display label for the selected alignment pass."""
    return f"{args.alignment_pass.capitalize()}-pass"


def coordinate_columns(args: argparse.Namespace) -> tuple[str, str]:
    """Return the coordinate columns required for the selected pass."""
    pass_name = alignment_pass_name(args)
    return f"x_{pass_name}", f"y_{pass_name}"


def display_sample_name(sample: str, donor: str) -> str:
    """Remove an exact donor prefix for plot display only."""
    prefix = f"{donor}_"
    return str(sample)[len(prefix):] if str(sample).startswith(prefix) else str(sample)


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    """Return all generated output paths."""
    donor = safe_filename(args.donor)
    pass_name = alignment_pass_name(args)
    return {
        "multiscale": args.output_dir / f"{donor}_{pass_name}_feature_similarity_multiscale_by_gene.csv",
        "scale": args.output_dir / f"{donor}_{pass_name}_feature_similarity_by_scale_gene.csv",
        "offset": args.output_dir / f"{donor}_{pass_name}_feature_similarity_by_offset_gene.csv",
        "pair_summary": args.output_dir / f"{donor}_{pass_name}_feature_similarity_pair_summary.csv",
        "metadata": args.output_dir / f"{donor}_{pass_name}_feature_similarity_metadata.json",
        "runtime": args.output_dir / f"{donor}_{pass_name}_feature_similarity_runtime.csv",
        "plot": args.plot_dir / f"{donor}_{pass_name}_feature_similarity_distributions.png",
        "overlay_plot": args.plot_dir / f"{donor}_{pass_name}_aligned_overlays.png",
    }


def prepare_output_dirs(args: argparse.Namespace) -> dict[str, Path]:
    """Create output directories and enforce overwrite protection."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(args)
    conflicts = [path for path in paths.values() if path.exists()]
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


def prepare_base_slices(
    adata: ad.AnnData, args: argparse.Namespace
) -> list[ad.AnnData]:
    """Subset one donor and create cell-ID-preserving slices in requested order."""
    required_obs = [
        args.donor_column, args.sample_column, args.adata_cell_id_column,
    ]
    missing = [column for column in required_obs if column not in adata.obs]
    if missing:
        raise ValueError(f"Missing required adata.obs columns: {missing}")
    donor_mask = (
        adata.obs[args.donor_column].astype("string").eq(str(args.donor))
        .fillna(False).to_numpy(dtype=bool)
    )
    if not donor_mask.any():
        raise ValueError(f"Donor {args.donor!r} is absent")
    donor_adata = adata[donor_mask].copy()
    donor_ids = donor_adata.obs[args.adata_cell_id_column].astype("string")
    if donor_ids.isna().any() or donor_ids.duplicated().any():
        raise ValueError("Cell IDs must be nonmissing and unique within the donor")
    donor_adata.obs_names = donor_ids.astype(str).to_numpy()
    sample_values = donor_adata.obs[args.sample_column].astype("string")
    observed = set(sample_values.dropna().astype(str))
    absent = [sample for sample in args.sample_order if sample not in observed]
    if absent:
        raise ValueError(f"Requested samples are absent: {absent}")
    return [
        donor_adata[
            sample_values.eq(sample).fillna(False).to_numpy(dtype=bool)
        ].copy()
        for sample in args.sample_order
    ]


def attach_coordinates_from_csv(
    slices: Sequence[ad.AnnData], args: argparse.Namespace
) -> list[ad.AnnData]:
    """Attach aligned coordinates after strict sample/cell-ID validation."""
    table = pd.read_csv(
        args.aligned_coordinate_csv, dtype={"cell_id": "string"}
    )
    x_column, y_column = coordinate_columns(args)
    required = {"donor", "sample", "cell_id", x_column, y_column}
    missing = required - set(table.columns)
    if missing:
        raise ValueError(f"Aligned-coordinate CSV is missing columns: {sorted(missing)}")
    table = table[table["donor"].astype(str).eq(str(args.donor))].copy()
    if table.empty:
        raise ValueError(f"Coordinate CSV contains no rows for donor {args.donor!r}")
    if table["cell_id"].isna().any() or table["cell_id"].duplicated().any():
        raise ValueError("Coordinate CSV cell_id values must be nonmissing and unique")
    coordinates = table.set_index("cell_id")
    conversion = 1000.0 if args.coordinate_unit == "mm" else 1.0
    results: list[ad.AnnData] = []
    for sample, ad_slice in zip(args.sample_order, slices):
        cell_ids = ad_slice.obs[args.adata_cell_id_column].astype("string")
        absent = cell_ids[~cell_ids.isin(coordinates.index)].astype(str).head(10).tolist()
        if absent:
            raise ValueError(
                f"Coordinate CSV is missing cells from {sample!r}; examples: {absent}"
            )
        rows = coordinates.loc[cell_ids.astype(str).to_numpy()]
        if not rows["sample"].astype(str).eq(sample).all():
            raise ValueError(f"Coordinate CSV sample labels disagree for {sample!r}")
        xy = rows[[x_column, y_column]].to_numpy(dtype=float) * conversion
        if not np.isfinite(xy).all():
            raise ValueError(f"Coordinate CSV contains nonfinite values for {sample!r}")
        result = ad_slice.copy()
        result.obsm[args.key_added] = xy
        results.append(result)
    return results


def flip_y_inplace(adata: ad.AnnData, spatial_key: str) -> None:
    """Reproduce the per-slice Y reflection before applying a transformation."""
    coordinates = np.asarray(adata.obsm[spatial_key], dtype=float).copy()
    if coordinates.ndim != 2 or coordinates.shape[1] < 2:
        raise ValueError(f"obsm[{spatial_key!r}] has unexpected shape {coordinates.shape}")
    if not np.isfinite(coordinates[:, :2]).all():
        raise ValueError(f"obsm[{spatial_key!r}] contains nonfinite coordinates")
    y = coordinates[:, 1]
    coordinates[:, 1] = y.max() + y.min() - y
    adata.obsm[spatial_key] = coordinates


def attach_coordinates_from_transformation(
    slices: Sequence[ad.AnnData],
    args: argparse.Namespace,
    spateo_module: Any | None = None,
) -> list[ad.AnnData]:
    """Apply a saved Spateo transformation and retain its aligned coordinates."""
    prepared = [ad_slice.copy() for ad_slice in slices]
    for sample, ad_slice in zip(args.sample_order, prepared):
        if args.spatial_key not in ad_slice.obsm:
            raise KeyError(f"Slice {sample!r} is missing obsm[{args.spatial_key!r}]")
        if args.flip_y:
            flip_y_inplace(ad_slice, args.spatial_key)
    with args.transformation_pkl.open("rb") as handle:
        transformation = pickle.load(handle)
    if spateo_module is None:
        import spateo as spateo_module
    aligned = spateo_module.align.morpho_align_apply_transformation(
        models=prepared,
        spatial_key=args.spatial_key,
        key_added=args.key_added,
        transformation=transformation,
        verbose=args.verbose,
    )
    if len(aligned) != len(prepared):
        raise ValueError("Applied transformation returned an unexpected number of slices")
    conversion = 1000.0 if args.coordinate_unit == "mm" else 1.0
    for sample, ad_slice in zip(args.sample_order, aligned):
        if args.key_added not in ad_slice.obsm:
            raise KeyError(f"Transformed slice {sample!r} lacks obsm[{args.key_added!r}]")
        xy = np.asarray(ad_slice.obsm[args.key_added], dtype=float)
        if xy.ndim != 2 or xy.shape[0] != ad_slice.n_obs or xy.shape[1] < 2:
            raise ValueError(f"Transformed coordinates for {sample!r} have shape {xy.shape}")
        xy = xy[:, :2].copy() * conversion
        if not np.isfinite(xy).all():
            raise ValueError(f"Transformed coordinates for {sample!r} are nonfinite")
        ad_slice.obsm[args.key_added] = xy
    return list(aligned)


def pair_columns(
    frame: pd.DataFrame, pair_number: int, reference: str, moving: str
) -> pd.DataFrame:
    """Prepend stable ordered-pair identifiers."""
    result = frame.copy()
    result.insert(0, "moving_slice", moving)
    result.insert(0, "reference_slice", reference)
    result.insert(0, "pair_id", f"pair_{pair_number:02d}_{reference}__{moving}")
    result.insert(0, "pair_number", pair_number)
    return result


def evaluate_consecutive_pairs(
    slices: Sequence[ad.AnnData], args: argparse.Namespace
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, list[dict[str, Any]]]:
    """Evaluate selected-pass coordinates for every consecutive pair."""
    pair_inputs = [
        (number, reference, moving)
        for number, (reference, moving) in enumerate(zip(slices[:-1], slices[1:]), 1)
    ]

    def evaluate_one_pair(
        pair_number: int, reference: ad.AnnData, moving: ad.AnnData
    ) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, Any], dict[str, Any]]:
        reference_name = args.sample_order[pair_number - 1]
        moving_name = args.sample_order[pair_number]
        LOGGER.info(
            "Evaluating %s pair %s -> %s",
            alignment_pass_label(args), reference_name, moving_name,
        )
        started = time.perf_counter()
        multiscale, scale, offset = compute_pair_gene_similarity(
            reference,
            moving,
            reference_coord_key=args.key_added,
            moving_coord_key=args.key_added,
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
        summaries = {
            metric: summarize_pair_gene_similarity(
                scale,
                primary_metric=metric,
                min_gene_nz_grids=args.min_gene_nz_grids,
                min_gene_valid_offsets=args.min_gene_valid_offsets,
                min_gene_valid_scales=args.min_gene_valid_scales,
            )
            for metric in PAIR_SUMMARY_METRICS
        }
        pair_id = f"pair_{pair_number:02d}_{reference_name}__{moving_name}"
        summary_row: dict[str, Any] = {
            "pair_number": pair_number, "pair_id": pair_id,
            "reference_slice": reference_name, "moving_slice": moving_name,
        }
        for metric, summary in summaries.items():
            for output_suffix, summary_key in (
                ("multiscale_median", "multiscale_median_primary_score"),
                ("n_genes_finite", "n_genes_primary_metric_finite"),
                ("n_genes_eligible", "n_genes_primary_metric_eligible"),
                ("finite_gene_fraction", "finite_gene_fraction"),
                ("eligible_gene_fraction", "eligible_gene_fraction"),
            ):
                summary_row[f"{metric}_{output_suffix}"] = summary[summary_key]
            summary_row[f"{metric}_per_scale_median"] = json.dumps(
                summary["per_scale_median_primary_score"], sort_keys=True
            )
            summary_row[f"{metric}_per_scale_n_genes_eligible"] = json.dumps(
                summary["per_scale_n_genes_eligible"], sort_keys=True
            )
            values = multiscale[metric].to_numpy(float)
            values = values[np.isfinite(values)]
            if values.size:
                q25, median, q75 = np.quantile(values, [0.25, 0.5, 0.75])
                summary_row[f"{metric}_gene_median"] = float(median)
                summary_row[f"{metric}_gene_q25"] = float(q25)
                summary_row[f"{metric}_gene_q75"] = float(q75)
                summary_row[f"{metric}_gene_iqr"] = float(q75 - q25)
                summary_row[f"{metric}_n_genes_distribution"] = int(values.size)
            else:
                for suffix in ("gene_median", "gene_q25", "gene_q75", "gene_iqr"):
                    summary_row[f"{metric}_{suffix}"] = np.nan
                summary_row[f"{metric}_n_genes_distribution"] = 0
        return (
            pair_columns(multiscale, pair_number, reference_name, moving_name),
            pair_columns(scale, pair_number, reference_name, moving_name),
            pair_columns(offset, pair_number, reference_name, moving_name),
            summary_row,
            {
                "pair_number": pair_number,
                "reference_slice": reference_name,
                "moving_slice": moving_name,
                "elapsed_seconds": time.perf_counter() - started,
                "n_reference_cells": reference.n_obs,
                "n_moving_cells": moving.n_obs,
                "n_genes_scored": multiscale["gene"].nunique(),
            },
        )

    requested_workers = (os.cpu_count() or 1) if args.n_jobs == -1 else args.n_jobs
    n_workers = min(len(pair_inputs), requested_workers)
    if n_workers == 1:
        completed = [evaluate_one_pair(*values) for values in pair_inputs]
    else:
        with threadpool_limits(limits=1):
            with concurrent.futures.ThreadPoolExecutor(max_workers=n_workers) as executor:
                completed = list(executor.map(lambda values: evaluate_one_pair(*values), pair_inputs))
    completed.sort(key=lambda value: int(value[3]["pair_number"]))
    return (
        pd.concat([value[0] for value in completed], ignore_index=True),
        pd.concat([value[1] for value in completed], ignore_index=True),
        pd.concat([value[2] for value in completed], ignore_index=True),
        pd.DataFrame.from_records([value[3] for value in completed]),
        [value[4] for value in completed],
    )


def plot_metric_distributions(
    result: pd.DataFrame,
    output_path: Path,
    donor: str,
    alignment_label: str,
    dpi: int,
    bins: int,
) -> None:
    """Plot pair rows by four metric columns using multiscale gene scores."""
    pairs = result[
        ["pair_number", "pair_id", "reference_slice", "moving_slice"]
    ].drop_duplicates().sort_values("pair_number")
    figure, axes = plt.subplots(
        len(pairs), len(PLOT_METRICS),
        figsize=(4.2 * len(PLOT_METRICS), 2.7 * len(pairs)), squeeze=False,
    )
    for row_index, pair in enumerate(pairs.itertuples(index=False)):
        pair_data = result[result["pair_id"] == pair.pair_id]
        reference = display_sample_name(pair.reference_slice, donor)
        moving = display_sample_name(pair.moving_slice, donor)
        for column_index, (metric, label) in enumerate(PLOT_METRICS):
            axis = axes[row_index, column_index]
            values = pair_data[metric].to_numpy(float)
            values = values[np.isfinite(values)]
            if values.size:
                median = float(np.median(values))
                axis.hist(values, bins=bins, color="#4C78A8", alpha=0.85)
                axis.axvline(median, color="#D62728", linestyle="--", linewidth=1)
                axis.text(
                    0.98, 0.95, f"n={values.size}\nmedian={median:.3g}",
                    transform=axis.transAxes, ha="right", va="top", fontsize=8,
                )
            else:
                axis.text(0.5, 0.5, "No finite values", ha="center", va="center")
            axis.set_title(f"{label}\n{reference} → {moving}", fontsize=9)
            axis.set_ylabel("Gene count")
            if row_index == len(pairs) - 1:
                axis.set_xlabel(metric)
    figure.suptitle(
        f"{alignment_label} multiscale feature similarity across all genes",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_aligned_overlays(
    slices: Sequence[ad.AnnData], args: argparse.Namespace, output_path: Path
) -> None:
    """Plot aligned consecutive pairs in the common blue/orange style."""
    n_pairs = len(slices) - 1
    n_columns = min(3, n_pairs)
    n_rows = math.ceil(n_pairs / n_columns)
    figure, axes = plt.subplots(
        n_rows, n_columns, figsize=(5.0 * n_columns, 4.5 * n_rows), squeeze=False,
    )
    axes_flat = axes.ravel()
    for index, axis in enumerate(axes_flat[:n_pairs]):
        reference_xy = np.asarray(slices[index].obsm[args.key_added])[:, :2]
        moving_xy = np.asarray(slices[index + 1].obsm[args.key_added])[:, :2]
        reference = display_sample_name(args.sample_order[index], args.donor)
        moving = display_sample_name(args.sample_order[index + 1], args.donor)
        axis.scatter(
            reference_xy[:, 0], reference_xy[:, 1], s=0.2, alpha=0.45,
            color="#4C78A8", linewidths=0, label=reference,
        )
        axis.scatter(
            moving_xy[:, 0], moving_xy[:, 1], s=0.2, alpha=0.45,
            color="#F58518", linewidths=0, label=moving,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(
            f"{alignment_pass_label(args)} aligned: {reference} → {moving}",
            fontsize=9,
        )
        axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    for axis in axes_flat[n_pairs:]:
        axis.axis("off")
    figure.suptitle(
        f"{alignment_pass_label(args)} aligned consecutive-slice overlays",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=args.plot_dpi, bbox_inches="tight")
    plt.close(figure)


def run_self_tests() -> None:
    """Test both coordinate sources and the complete evaluation integration."""
    coordinates = np.array([(x, y) for y in range(4) for x in range(4)], float)
    expression = np.column_stack(
        [coordinates[:, 0] + coordinates[:, 1] + 1, coordinates[:, 0] + 1]
    )
    adata = ad.AnnData(
        X=np.vstack([expression, expression, expression]),
        obs=pd.DataFrame(
            {
                "Donor": ["D1"] * 48,
                "Sample": np.repeat(["D1_S1", "D1_S2", "D1_S3"], 16),
                "cell_id": [f"c{index}" for index in range(48)],
            }
        ),
        obsm={"spatial": np.vstack([coordinates, coordinates, coordinates])},
    )
    adata.var_names = ["g1", "g2"]
    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary = Path(temporary_directory)
        csv_path = temporary / "aligned.csv"
        csv_rows = []
        for sample_index, sample in enumerate(("D1_S1", "D1_S2", "D1_S3")):
            for within_index, xy in enumerate(coordinates):
                csv_rows.append(
                    {
                        "donor": "D1", "sample": sample,
                        "cell_id": f"c{sample_index * 16 + within_index}",
                        "x_first_pass": xy[0] + sample_index,
                        "y_first_pass": xy[1],
                        "x_second_pass": xy[0] + sample_index + 100,
                        "y_second_pass": xy[1] + 100,
                    }
                )
        pd.DataFrame(csv_rows[::-1]).to_csv(csv_path, index=False)
        pkl_path = temporary / "transformation.pkl"
        with pkl_path.open("wb") as handle:
            pickle.dump({"test": True}, handle)
        args = argparse.Namespace(
            alignment_pass="first", donor="D1", donor_column="Donor",
            sample_column="Sample",
            sample_order=["D1_S1", "D1_S2", "D1_S3"],
            adata_cell_id_column="cell_id", spatial_key="spatial",
            key_added="align_spatial", coordinate_unit="um", flip_y=True,
            verbose=False, aligned_coordinate_csv=csv_path,
            transformation_pkl=None, grid_sizes_um=[1.0, 2.0],
            min_cells_reference=1, min_cells_moving=1, min_overlap_grids=2,
            min_gene_nz_grids=1, min_gene_valid_offsets=1,
            min_gene_valid_scales=1, random_state=0, n_jobs=2,
            plot_dpi=80, output_dir=temporary / "results",
            plot_dir=temporary / "plots",
        )
        base = prepare_base_slices(adata, args)
        csv_slices = attach_coordinates_from_csv(base, args)
        assert np.allclose(csv_slices[1].obsm["align_spatial"][:, 0], coordinates[:, 0] + 1)
        assert coordinate_columns(args) == ("x_first_pass", "y_first_pass")
        assert "first_pass" in output_paths(args)["multiscale"].name

        second_args = argparse.Namespace(**vars(args))
        second_args.alignment_pass = "second"
        second_slices = attach_coordinates_from_csv(base, second_args)
        assert np.allclose(
            second_slices[1].obsm["align_spatial"][:, 0],
            coordinates[:, 0] + 101,
        )
        assert coordinate_columns(second_args) == (
            "x_second_pass", "y_second_pass",
        )
        assert "second_pass" in output_paths(second_args)["multiscale"].name

        class FakeAlign:
            @staticmethod
            def morpho_align_apply_transformation(**kwargs: Any) -> list[ad.AnnData]:
                output = [value.copy() for value in kwargs["models"]]
                for index, value in enumerate(output):
                    value.obsm[kwargs["key_added"]] = (
                        np.asarray(value.obsm[kwargs["spatial_key"]])[:, :2] + index
                    )
                return output

        class FakeSpateo:
            align = FakeAlign()

        args.aligned_coordinate_csv = None
        args.transformation_pkl = pkl_path
        pkl_slices = attach_coordinates_from_transformation(base, args, FakeSpateo())
        assert all("align_spatial" in value.obsm for value in pkl_slices)
        args.aligned_coordinate_csv = csv_path
        args.transformation_pkl = None
        multiscale, scale, offset, summary, _ = evaluate_consecutive_pairs(
            csv_slices, args
        )
        assert multiscale["pair_id"].nunique() == 2
        assert scale["pair_id"].nunique() == 2
        assert offset["pair_id"].nunique() == 2
        assert len(summary) == 2
        for metric in PAIR_SUMMARY_METRICS:
            assert f"{metric}_multiscale_median" in summary
        sequential_args = argparse.Namespace(**vars(args))
        sequential_args.n_jobs = 1
        sequential_multiscale, _, _, sequential_summary, _ = (
            evaluate_consecutive_pairs(csv_slices, sequential_args)
        )
        pd.testing.assert_frame_equal(multiscale, sequential_multiscale)
        pd.testing.assert_frame_equal(summary, sequential_summary)
        distribution_path = temporary / "distribution.png"
        overlay_path = temporary / "overlay.png"
        plot_metric_distributions(
            multiscale, distribution_path, "D1",
            alignment_pass_label(args), 80, 5,
        )
        plot_aligned_overlays(csv_slices, args, overlay_path)
        assert distribution_path.is_file() and overlay_path.is_file()
    LOGGER.info("All self-tests passed")


def main(argv: Sequence[str] | None = None) -> int:
    """Run Spateo feature-similarity evaluation for one alignment pass."""
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
    )
    if args.run_self_tests:
        run_self_tests()
        return 0
    for path, label in (
        (args.input_h5ad, "Input AnnData"),
        (args.aligned_coordinate_csv, "Aligned-coordinate CSV"),
        (args.transformation_pkl, "Transformation PKL"),
    ):
        if path is not None and not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    paths = prepare_output_dirs(args)
    total_started = time.perf_counter()
    LOGGER.info("Loading original expression AnnData: %s", args.input_h5ad)
    adata = ad.read_h5ad(args.input_h5ad)
    base_slices = prepare_base_slices(adata, args)
    del adata
    if args.aligned_coordinate_csv is not None:
        coordinate_source = "aligned_coordinate_csv"
        aligned_slices = attach_coordinates_from_csv(base_slices, args)
    else:
        coordinate_source = "reapplied_spateo_transformation_pkl"
        aligned_slices = attach_coordinates_from_transformation(base_slices, args)
    multiscale, scale, offset, summary, runtime = evaluate_consecutive_pairs(
        aligned_slices, args
    )
    atomic_csv_dump(multiscale, paths["multiscale"])
    atomic_csv_dump(scale, paths["scale"])
    atomic_csv_dump(offset, paths["offset"])
    atomic_csv_dump(summary, paths["pair_summary"])
    plot_metric_distributions(
        multiscale, paths["plot"], args.donor, alignment_pass_label(args),
        args.plot_dpi, args.histogram_bins,
    )
    plot_aligned_overlays(aligned_slices, args, paths["overlay_plot"])
    runtime.append(
        {
            "pair_number": np.nan, "reference_slice": "ALL", "moving_slice": "ALL",
            "elapsed_seconds": time.perf_counter() - total_started,
            "n_reference_cells": np.nan, "n_moving_cells": np.nan,
            "n_genes_scored": multiscale["gene"].nunique(),
        }
    )
    atomic_csv_dump(pd.DataFrame.from_records(runtime), paths["runtime"])
    metadata = {
        "script_version": SCRIPT_VERSION,
        "input_h5ad": str(args.input_h5ad.resolve()),
        "alignment_pass": args.alignment_pass,
        "coordinate_columns": list(coordinate_columns(args)),
        "coordinate_source": coordinate_source,
        "aligned_coordinate_csv": (
            str(args.aligned_coordinate_csv.resolve())
            if args.aligned_coordinate_csv is not None else None
        ),
        "transformation_pkl": (
            str(args.transformation_pkl.resolve())
            if args.transformation_pkl is not None else None
        ),
        "donor": str(args.donor),
        "donor_column": args.donor_column,
        "sample_column": args.sample_column,
        "sample_order": list(args.sample_order),
        "adata_cell_id_column": args.adata_cell_id_column,
        "original_spatial_key_used_for_evaluation": False,
        "spatial_key_for_transformation_reapplication": (
            args.spatial_key if args.transformation_pkl is not None else None
        ),
        "aligned_evaluation_key": args.key_added,
        "coordinate_unit_input": args.coordinate_unit,
        "coordinate_unit_analysis": "um",
        "flip_y_before_transformation_reapplication": (
            args.flip_y if args.transformation_pkl is not None else None
        ),
        "expression_source": "original_adata.X",
        "gene_selection": "all_genes_shared_by_each_pair",
        "metrics": list(SUPPORTED_METRICS),
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
        "outputs": {name: str(path) for name, path in paths.items()},
        "software_versions": {
            "python": platform.python_version(), "anndata": ad.__version__,
            "numpy": np.__version__, "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    atomic_json_dump(metadata, paths["metadata"])
    LOGGER.info("Module 03 complete. Outputs: %s", args.output_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())


# Example using an aligned-coordinate CSV:
# python 03_Spateo_feature_similarity.py \
#     --input-h5ad /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/processed-data/02_build_spe/h5ad/spe_NormCounts_nucleus_normcounts.h5ad \
#     --alignment-pass first \
#     --aligned-coordinate-csv /path/to/Br6660_first_pass_coordinates.csv \
#     --output-dir /path/to/module03_Spateo_feature_similarity_first_pass_Br6660 \
#     --plot-dir /path/to/module03_Spateo_feature_similarity_first_pass_Br6660 \
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
#     --spatial-key spatial \
#     --key-added align_spatial \
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
#
# Alternative PKL mode: replace --aligned-coordinate-csv above with:
#     --transformation-pkl /path/to/Spateo_transformation_Br6660.pkl \
#     --flip-y
