#!/usr/bin/env python3
"""Permutation-free diagnosis of CRAWDAD neighborhood distance for D1 islands.

This module measures directional cell-level proximity, neighborhood coverage,
and grid-defined island geometry.  It deliberately does not call CRAWDAD or
perform permutations; null-shuffle scale is outside its scope.
"""

from __future__ import annotations

import argparse
import copy
import concurrent.futures
import fcntl
import hashlib
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile
from collections.abc import Iterable, Sequence
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-distance-qc")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Circle
from scipy import ndimage
from scipy.spatial import cKDTree, distance as scipy_distance
from shapely import union_all
from shapely.geometry import Polygon, box


CANONICAL_SAMPLES = (
    "Br6660_NAc1_580",
    "Br6660_NAc2_1090",
    "Br6660_NAc3_1580",
    "Br6660_NAc4_2080",
    "Br6660_NAc5_2580",
    "Br6660_NAc6_3080",
    "Br6660_NAc7_3580",
    "Br6660_Nac10_4080",
    "Br6660_NAc8_4580",
    "Br6660_NAc9_5080",
    "Br6660_Nac11_5580",
)
DEFAULT_REGIONS = (
    "lateral",
    "dorsomedial",
    "ventromedial",
    "outside_global_roi",
)
CELLTYPE_COLORS = {"D1_Island_A": "#224B82", "D1_Island_B": "#00A650"}
NN_COLUMNS = (
    "sample", "region", "reference_celltype", "target_celltype",
    "reference_cell_id", "nearest_distance_um",
)
COVERAGE_COLUMNS = (
    "sample", "region", "reference_celltype", "target_celltype",
    "candidate_distance_um", "n_reference_cells", "n_target_cells",
    "coverage_ge_1", "mean_target_neighbors", "median_target_neighbors",
    "p25_target_neighbors", "p75_target_neighbors", "p90_target_neighbors",
    "fraction_with_ge_5", "fraction_with_ge_10",
)
ISLAND_COLUMNS = (
    "sample", "region", "cell_type", "island_id", "component_label",
    "n_cells", "n_occupied_bins", "occupied_area_um2",
    "equivalent_diameter_um", "centroid_x_um", "centroid_y_um",
    "bbox_width_um", "bbox_height_um", "maximum_span_um",
    "maximum_span_method", "eligible_island", "polygon_wkt",
)
OPPOSITE_COLUMNS = (
    "sample", "region", "reference_celltype", "target_celltype",
    "reference_island_id", "nearest_target_island_id",
    "boundary_distance_um", "centroid_distance_um", "touches_or_overlaps",
)
ISLAND_GAP_OUT_OF_RANGE_REASON = (
    "Median A/B island boundary gap exceeds the tested candidate-distance "
    "range; expand the candidate range before making a recommendation."
)
RECOMMENDATION_LOGIC_VERSION = "island_gap_candidate_range_guard_v2"


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Diagnose a global CRAWDAD neighborhood distance for "
            "D1_Island_A/B without permutations."
        )
    )
    inputs = parser.add_argument_group("inputs")
    inputs.add_argument("--input-cell-table", type=Path)
    inputs.add_argument("--aligned-coordinates", type=Path)
    inputs.add_argument("--celltype-annotations", type=Path)
    inputs.add_argument("--roi-assignments", type=Path)
    inputs.add_argument("--coordinate-unit", default="micrometer")
    parser.add_argument("--donor", default="Br6660")
    parser.add_argument("--samples", nargs="+", default=["all"])
    parser.add_argument("--regions", nargs="+", default=list(DEFAULT_REGIONS))
    parser.add_argument("--reference-celltype", default="D1_Island_A")
    parser.add_argument("--target-celltype", default="D1_Island_B")
    parser.add_argument(
        "--candidate-distances-um", nargs="+", type=float,
        default=[25, 50, 75, 100, 150, 200, 300],
    )
    parser.add_argument("--island-bin-size-um", type=float, default=25)
    parser.add_argument("--island-connectivity", type=int, choices=(4, 8), default=8)
    parser.add_argument("--island-closing-bins", type=int, default=1)
    parser.add_argument("--min-island-cells", type=int, default=20)
    parser.add_argument("--min-cells-per-type", type=int, default=20)
    parser.add_argument("--coverage-saturation-threshold", type=float, default=0.90)
    parser.add_argument("--median-count-saturation", type=float, default=10)
    parser.add_argument("--p90-count-saturation", type=float, default=50)
    parser.add_argument("--directionality-threshold", type=float, default=0.30)
    parser.add_argument("--global-support-fraction", type=float, default=0.60)
    parser.add_argument("--global-max-saturation-fraction", type=float, default=0.25)
    parser.add_argument("--elbow-min-score", type=float, default=0.05)
    parser.add_argument("--run-dbscan-sensitivity", action="store_true")
    parser.add_argument("--dbscan-eps-um", type=float, default=50)
    parser.add_argument("--dbscan-min-samples", type=int, default=20)
    parser.add_argument("--n-jobs", type=int, default=11)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--max-plot-cells", type=int, default=150_000)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument(
        "--plot-dir", type=Path,
        help="Root directory for PNG figures; tables and reports remain in --output-dir.",
    )
    parser.add_argument(
        "--regional-job", action="store_true",
        help=(
            "Run exactly one region in a shared output root. Region jobs may "
            "run concurrently; the last compatible job automatically builds "
            "the global 44-task summary."
        ),
    )
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--resume", action="store_true")
    mode.add_argument("--overwrite", action="store_true")
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)

    if args.self_test:
        return args
    separate = (
        args.aligned_coordinates,
        args.celltype_annotations,
        args.roi_assignments,
    )
    if args.input_cell_table and any(x is not None for x in separate):
        parser.error("Use --input-cell-table or the three separate inputs, not both")
    if not args.input_cell_table and not all(x is not None for x in separate):
        parser.error(
            "Provide --input-cell-table or all of --aligned-coordinates, "
            "--celltype-annotations, and --roi-assignments"
        )
    if args.output_dir is None:
        parser.error("--output-dir is required")
    if args.plot_dir is None:
        parser.error("--plot-dir is required")
    if str(args.coordinate_unit).lower() not in {
        "um", "µm", "micrometer", "micrometers", "micron", "microns"
    }:
        parser.error("--coordinate-unit must denote micrometers")
    positive = {
        "candidate distances": args.candidate_distances_um,
        "island bin size": [args.island_bin_size_um],
        "minimum island cells": [args.min_island_cells],
        "minimum cells per type": [args.min_cells_per_type],
        "n-jobs": [args.n_jobs],
        "dpi": [args.dpi],
    }
    for label, values in positive.items():
        if any(not np.isfinite(x) or x <= 0 for x in values):
            parser.error(f"{label} must be positive")
    args.candidate_distances_um = sorted(set(args.candidate_distances_um))
    if len(args.regions) != len(set(args.regions)):
        parser.error("--regions contains duplicates")
    if args.regional_job and len(args.regions) != 1:
        parser.error("--regional-job requires exactly one --regions value")
    if args.reference_celltype == args.target_celltype:
        parser.error("Reference and target cell types must differ")
    if args.island_closing_bins < 0:
        parser.error("--island-closing-bins must be >= 0")
    for name in (
        "coverage_saturation_threshold", "directionality_threshold",
        "global_support_fraction", "global_max_saturation_fraction",
        "elbow_min_score",
    ):
        value = getattr(args, name)
        if not 0 <= value <= 1:
            parser.error(f"--{name.replace('_', '-')} must be in [0, 1]")
    return args


def read_table(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.name.endswith(".csv.gz") or path.suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported table format: {path}")


def canonicalize_columns(frame: pd.DataFrame) -> pd.DataFrame:
    aliases = {
        "cell_id": ("cell_id", "cell", "barcode"),
        "sample": ("sample", "sample_id"),
        "cell_type": ("cell_type", "celltype", "annotation"),
        "x_um": ("x_um", "x_aligned", "x_second_pass", "x"),
        "y_um": ("y_um", "y_aligned", "y_second_pass", "y"),
        "region": ("region", "primary_roi", "analysis_region"),
        "donor": ("donor",),
    }
    rename: dict[str, str] = {}
    for canonical, candidates in aliases.items():
        present = [x for x in candidates if x in frame.columns]
        if len(present) > 1 and canonical not in present:
            raise ValueError(f"Ambiguous columns for {canonical}: {present}")
        if present:
            rename[present[0]] = canonical
    return frame.rename(columns=rename)


def load_cells(args: argparse.Namespace) -> tuple[pd.DataFrame, list[str]]:
    sources: list[str] = []
    if args.input_cell_table:
        cells = canonicalize_columns(read_table(args.input_cell_table))
        sources.append(str(args.input_cell_table.resolve()))
    else:
        tables = []
        for path in (
            args.aligned_coordinates,
            args.celltype_annotations,
            args.roi_assignments,
        ):
            table = canonicalize_columns(read_table(path))
            if "cell_id" not in table:
                raise ValueError(f"{path} has no cell_id column")
            if table["cell_id"].duplicated().any():
                raise ValueError(f"{path} has duplicated cell_id values")
            tables.append(table)
            sources.append(str(path.resolve()))
        cells = tables[0]
        for table in tables[1:]:
            overlap = (set(cells) & set(table)) - {"cell_id"}
            if overlap:
                raise ValueError(
                    "Separate input tables contain overlapping non-key columns: "
                    + ", ".join(sorted(overlap))
                )
            cells = cells.merge(table, on="cell_id", how="inner", validate="one_to_one")

    required = {"cell_id", "sample", "cell_type", "x_um", "y_um", "region"}
    missing = sorted(required - set(cells))
    if missing:
        raise ValueError(f"Cell table is missing required columns: {missing}")
    cells = cells.copy()
    for column in ("cell_id", "sample", "cell_type", "region"):
        if cells[column].isna().any():
            raise ValueError(f"{column} contains missing values")
        cells[column] = cells[column].astype(str)
    for column in ("x_um", "y_um"):
        cells[column] = pd.to_numeric(cells[column], errors="coerce")
        if not np.isfinite(cells[column]).all():
            raise ValueError(f"{column} contains non-finite coordinates")
    if "donor" in cells:
        cells = cells.loc[cells["donor"].astype(str) == args.donor].copy()
    if cells.empty:
        raise ValueError(f"No cells remain for donor {args.donor}")

    duplicated = cells.loc[cells["cell_id"].duplicated(False)]
    if not duplicated.empty:
        conflicts = duplicated.groupby("cell_id", sort=False)[
            ["sample", "cell_type", "region", "x_um", "y_um"]
        ].nunique(dropna=False)
        if (conflicts > 1).any(axis=None):
            raise ValueError("At least one cell_id maps to conflicting attributes")
        cells = cells.drop_duplicates("cell_id")

    available_samples = set(cells["sample"])
    if args.samples == ["all"]:
        missing_canonical = [x for x in CANONICAL_SAMPLES if x not in available_samples]
        if missing_canonical:
            raise ValueError(f"Missing canonical samples: {missing_canonical}")
        args.samples = list(CANONICAL_SAMPLES)
    else:
        unknown = sorted(set(args.samples) - available_samples)
        if unknown:
            raise ValueError(f"Requested samples not found: {unknown}")
    unknown_regions = sorted(set(args.regions) - set(cells["region"]))
    if unknown_regions:
        raise ValueError(f"Requested regions not found: {unknown_regions}")
    cells = cells.loc[
        cells["sample"].isin(args.samples) & cells["region"].isin(args.regions)
    ].copy()
    return cells, sources


def atomic_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_csv(tmp, index=False)
    os.replace(tmp, path)


def atomic_parquet(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    frame.to_parquet(tmp, index=False)
    os.replace(tmp, path)


def atomic_json(value: Any, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    tmp.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(tmp, path)


def empty_frame(columns: Iterable[str]) -> pd.DataFrame:
    return pd.DataFrame(columns=list(columns))


def task_key(sample: str, region: str) -> str:
    return f"{region}__{sample}"


def task_hash(args_dict: dict[str, Any], sample: str, region: str) -> str:
    relevant = {
        key: args_dict[key]
        for key in (
            "reference_celltype", "target_celltype", "candidate_distances_um",
            "island_bin_size_um", "island_connectivity", "island_closing_bins",
            "min_island_cells", "min_cells_per_type",
            "coverage_saturation_threshold", "median_count_saturation",
            "p90_count_saturation", "directionality_threshold", "elbow_min_score",
            "run_dbscan_sensitivity", "dbscan_eps_um", "dbscan_min_samples",
            "seed", "max_plot_cells", "dpi", "input_identity",
        )
    }
    relevant.update(sample=sample, region=region)
    relevant["recommendation_logic_version"] = RECOMMENDATION_LOGIC_VERSION
    return hashlib.sha256(
        json.dumps(relevant, sort_keys=True).encode()
    ).hexdigest()


def nearest_distances(
    reference_xy: np.ndarray, target_xy: np.ndarray, exclude_self: bool
) -> np.ndarray:
    if not len(reference_xy) or not len(target_xy):
        return np.full(len(reference_xy), np.nan)
    tree = cKDTree(target_xy)
    if exclude_self:
        if len(target_xy) < 2:
            return np.full(len(reference_xy), np.nan)
        distances, _ = tree.query(reference_xy, k=2, workers=1)
        return distances[:, 1]
    distances, _ = tree.query(reference_xy, k=1, workers=1)
    return np.asarray(distances)


def radius_counts(
    reference_xy: np.ndarray, target_xy: np.ndarray, radius: float,
    exclude_self: bool,
) -> np.ndarray:
    if not len(reference_xy) or not len(target_xy):
        return np.zeros(len(reference_xy), dtype=int)
    values = cKDTree(target_xy).query_ball_point(
        reference_xy, radius, return_length=True, workers=1
    ).astype(int)
    if exclude_self:
        values = np.maximum(values - 1, 0)
    return values


def connectivity_structure(connectivity: int) -> np.ndarray:
    return ndimage.generate_binary_structure(2, 2 if connectivity == 8 else 1)


def maximum_span(polygon: Polygon) -> tuple[float, str]:
    hull = polygon.convex_hull
    if not hasattr(hull, "exterior"):
        return 0.0, "degenerate"
    coords = np.asarray(hull.exterior.coords[:-1], dtype=float)
    if len(coords) <= 1:
        return 0.0, "degenerate"
    if len(coords) > 2000:
        index = np.linspace(0, len(coords) - 1, 2000, dtype=int)
        coords = coords[index]
        method = "deterministic_2000_vertex_approximation"
    else:
        method = "exact_convex_hull_vertices"
    return float(np.max(scipy_distance.pdist(coords))), method


def build_islands(
    task: pd.DataFrame, cell_type: str, args: dict[str, Any],
    origin: tuple[float, float], shape: tuple[int, int], roi_mask: np.ndarray,
) -> tuple[pd.DataFrame, list[dict[str, Any]], dict[str, int]]:
    selected = task.loc[task["cell_type"] == cell_type]
    bin_size = args["island_bin_size_um"]
    ix = np.floor((selected["x_um"].to_numpy() - origin[0]) / bin_size).astype(int)
    iy = np.floor((selected["y_um"].to_numpy() - origin[1]) / bin_size).astype(int)
    occupied = np.zeros(shape, dtype=bool)
    occupied[iy, ix] = True
    raw_labels, raw_n = ndimage.label(
        occupied, structure=connectivity_structure(args["island_connectivity"])
    )
    if args["island_closing_bins"] > 0:
        closed = ndimage.binary_closing(
            occupied,
            structure=np.ones((3, 3), dtype=bool),
            iterations=args["island_closing_bins"],
            border_value=0,
        )
        closed &= roi_mask
        closed |= occupied
    else:
        closed = occupied.copy()
    labels, n_components = ndimage.label(
        closed, structure=connectivity_structure(args["island_connectivity"])
    )
    cell_labels = labels[iy, ix]
    raw_cell_labels = raw_labels[iy, ix]
    cells_connected = 0
    for closed_label in np.unique(cell_labels):
        members = cell_labels == closed_label
        if len(np.unique(raw_cell_labels[members])) > 1:
            cells_connected += int(members.sum())
    rows: list[dict[str, Any]] = []
    internal: list[dict[str, Any]] = []
    prefix = "A" if cell_type == args["reference_celltype"] else "B"
    for label_value in range(1, n_components + 1):
        ys, xs = np.where(labels == label_value)
        member = selected.iloc[np.flatnonzero(cell_labels == label_value)]
        boxes = [
            box(
                origin[0] + x * bin_size,
                origin[1] + y * bin_size,
                origin[0] + (x + 1) * bin_size,
                origin[1] + (y + 1) * bin_size,
            )
            for y, x in zip(ys, xs, strict=True)
        ]
        polygon = union_all(boxes)
        island_id = f"{prefix}_{label_value:03d}"
        n_cells = len(member)
        area = float(len(boxes) * bin_size**2)
        minx, miny, maxx, maxy = polygon.bounds
        span, span_method = maximum_span(polygon)
        eligible = n_cells >= args["min_island_cells"]
        row = {
            "sample": task["sample"].iloc[0],
            "region": task["region"].iloc[0],
            "cell_type": cell_type,
            "island_id": island_id,
            "component_label": label_value,
            "n_cells": n_cells,
            "n_occupied_bins": len(boxes),
            "occupied_area_um2": area,
            "equivalent_diameter_um": 2 * math.sqrt(area / math.pi),
            "centroid_x_um": polygon.centroid.x,
            "centroid_y_um": polygon.centroid.y,
            "bbox_width_um": maxx - minx,
            "bbox_height_um": maxy - miny,
            "maximum_span_um": span,
            "maximum_span_method": span_method,
            "eligible_island": eligible,
            "polygon_wkt": polygon.wkt,
        }
        rows.append(row)
        internal.append({**row, "polygon": polygon})
    qc = {
        "occupied_bins_before_closing": int(occupied.sum()),
        "occupied_bins_after_closing": int(closed.sum()),
        "bins_added_by_closing": int((closed & ~occupied).sum()),
        "raw_component_count": int(raw_n),
        "closed_component_count": int(n_components),
        "components_connected_by_closing": max(0, int(raw_n - n_components)),
        "cells_connected_by_closing": cells_connected,
    }
    return pd.DataFrame(rows, columns=ISLAND_COLUMNS), internal, qc


def opposite_island_distances(
    source: list[dict[str, Any]], target: list[dict[str, Any]],
    source_type: str, target_type: str, sample: str, region: str,
) -> pd.DataFrame:
    rows: list[dict[str, Any]] = []
    eligible_source = [x for x in source if x["eligible_island"]]
    eligible_target = [x for x in target if x["eligible_island"]]
    for ref in eligible_source:
        if not eligible_target:
            continue
        distances = [ref["polygon"].distance(x["polygon"]) for x in eligible_target]
        chosen = eligible_target[int(np.argmin(distances))]
        rows.append({
            "sample": sample,
            "region": region,
            "reference_celltype": source_type,
            "target_celltype": target_type,
            "reference_island_id": ref["island_id"],
            "nearest_target_island_id": chosen["island_id"],
            "boundary_distance_um": float(ref["polygon"].distance(chosen["polygon"])),
            "centroid_distance_um": float(ref["polygon"].centroid.distance(chosen["polygon"].centroid)),
            "touches_or_overlaps": bool(ref["polygon"].distance(chosen["polygon"]) == 0),
        })
    return pd.DataFrame(rows, columns=OPPOSITE_COLUMNS)


def coverage_summary(
    task: pd.DataFrame, args: dict[str, Any]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    sample = task["sample"].iloc[0]
    region = task["region"].iloc[0]
    types = (args["reference_celltype"], args["target_celltype"])
    subsets = {x: task.loc[task["cell_type"] == x] for x in types}
    coordinates = {x: subsets[x][["x_um", "y_um"]].to_numpy() for x in types}
    nearest_rows = []
    coverage_rows = []
    for ref_type, target_type in (
        (types[0], types[1]), (types[1], types[0]),
        (types[0], types[0]), (types[1], types[1]),
    ):
        exclude_self = ref_type == target_type
        ref = subsets[ref_type]
        ref_xy = coordinates[ref_type]
        target_xy = coordinates[target_type]
        nearest = nearest_distances(ref_xy, target_xy, exclude_self)
        nearest_rows.append(pd.DataFrame({
            "sample": sample,
            "region": region,
            "reference_celltype": ref_type,
            "target_celltype": target_type,
            "reference_cell_id": ref["cell_id"].to_numpy(),
            "nearest_distance_um": nearest,
        }))
        for radius in args["candidate_distances_um"]:
            counts = radius_counts(ref_xy, target_xy, radius, exclude_self)
            coverage_rows.append({
                "sample": sample,
                "region": region,
                "reference_celltype": ref_type,
                "target_celltype": target_type,
                "candidate_distance_um": radius,
                "n_reference_cells": len(ref_xy),
                "n_target_cells": len(target_xy),
                "coverage_ge_1": float(np.mean(counts >= 1)),
                "mean_target_neighbors": float(np.mean(counts)),
                "median_target_neighbors": float(np.median(counts)),
                "p25_target_neighbors": float(np.quantile(counts, 0.25)),
                "p75_target_neighbors": float(np.quantile(counts, 0.75)),
                "p90_target_neighbors": float(np.quantile(counts, 0.90)),
                "fraction_with_ge_5": float(np.mean(counts >= 5)),
                "fraction_with_ge_10": float(np.mean(counts >= 10)),
            })
    return (
        pd.concat(nearest_rows, ignore_index=True),
        pd.DataFrame(coverage_rows, columns=COVERAGE_COLUMNS),
    )


def stable_elbow(curve: pd.DataFrame, min_score: float) -> float | None:
    curve = curve.sort_values("candidate_distance_um")
    x = curve["candidate_distance_um"].to_numpy(float)
    y = curve["coverage_ge_1"].to_numpy(float)
    if len(x) < 3 or x[-1] == x[0] or np.ptp(y) < 0.05:
        return None
    x_norm = (x - x[0]) / (x[-1] - x[0])
    y_norm = (y - y[0]) / np.ptp(y)
    score = y_norm - x_norm
    index = int(np.argmax(score))
    if index in (0, len(x) - 1) or score[index] < min_score:
        return None
    return float(x[index])


def recommend_task(
    coverage: pd.DataFrame, opposite: pd.DataFrame, args: dict[str, Any],
    sample: str, region: str,
) -> dict[str, Any]:
    a, b = args["reference_celltype"], args["target_celltype"]
    cross = coverage.loc[coverage["reference_celltype"] != coverage["target_celltype"]]
    pivot = cross.pivot_table(
        index="candidate_distance_um", columns="reference_celltype",
        values="coverage_ge_1", aggfunc="first",
    ).sort_index()
    mean_curve = pivot.mean(axis=1).rename("coverage_ge_1").reset_index()
    elbow = stable_elbow(mean_curve, args["elbow_min_score"])
    maximum_candidate = float(max(args["candidate_distances_um"]))
    finite_gaps = pd.to_numeric(
        opposite.get("boundary_distance_um", pd.Series(dtype=float)),
        errors="coerce",
    )
    finite_gaps = finite_gaps[np.isfinite(finite_gaps)]
    median_gap = float(finite_gaps.median()) if len(finite_gaps) else None
    gap_candidate: float | None = None
    if median_gap is not None and median_gap <= maximum_candidate:
        candidates = [x for x in args["candidate_distances_um"] if x >= median_gap]
        if candidates:
            gap_candidate = float(min(candidates))
    supports = [x for x in (elbow, gap_candidate) if x is not None]
    support_threshold = max(supports) if supports else None
    status = "eligible"
    reason = "Candidate supported by island-gap and coverage/saturation QC."
    recommendation: float | None = None
    if median_gap is not None and median_gap > maximum_candidate:
        status = "island_gap_exceeds_candidate_range"
        reason = ISLAND_GAP_OUT_OF_RANGE_REASON
    elif median_gap is None:
        status = "insufficient_data"
        reason = (
            "Median A/B island boundary gap is missing or non-finite; "
            "a distance recommendation cannot be made."
        )
    elif support_threshold is None:
        status = "no_stable_support"
        reason = "No stable coverage elbow or in-range island-gap support was available."
    else:
        supported = [x for x in args["candidate_distances_um"] if x >= support_threshold]
        for candidate in supported:
            rows = cross.loc[np.isclose(cross["candidate_distance_um"], candidate)]
            coverage_sat = bool(
                (rows["coverage_ge_1"] >= args["coverage_saturation_threshold"]).any()
            )
            count_sat = bool(
                (rows["median_target_neighbors"] >= args["median_count_saturation"]).any()
                or (rows["p90_target_neighbors"] >= args["p90_count_saturation"]).any()
            )
            if not coverage_sat and not count_sat:
                recommendation = float(candidate)
                status = "recommended"
                reason = "Smallest supported candidate before coverage/count saturation."
                break
        if recommendation is None:
            status = "ambiguous_large_d"
            reason = "Supported candidates occur only at or after coverage/count saturation."
    directionality = cross.pivot_table(
        index="candidate_distance_um", columns="reference_celltype",
        values="coverage_ge_1", aggfunc="first",
    )
    max_difference = (
        float(np.max(np.abs(directionality[a] - directionality[b])))
        if a in directionality and b in directionality else np.nan
    )
    return {
        "sample": sample,
        "region": region,
        "recommendation_status": status,
        "task_qc_status": status,
        "recommendation_reason": reason,
        "coverage_elbow_um": elbow,
        "median_nearest_opposite_island_boundary_gap_um": median_gap,
        "maximum_tested_candidate_distance_um": maximum_candidate,
        "nearest_island_median_gap_candidate_um": gap_candidate,
        "support_threshold_um": support_threshold,
        "recommended_distance_um": recommendation,
        "recommended_task_distance_um": recommendation,
        "maximum_directional_coverage_difference": max_difference,
        "strong_directionality": bool(
            np.isfinite(max_difference)
            and max_difference >= args["directionality_threshold"]
        ),
    }


def deterministic_subset(frame: pd.DataFrame, limit: int, seed: int) -> pd.DataFrame:
    if limit <= 0 or len(frame) <= limit:
        return frame
    return frame.sample(limit, random_state=seed)


def save_figure(fig: plt.Figure, path: Path, dpi: int) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(f".{path.stem}.tmp.{os.getpid()}.png")
    fig.savefig(tmp, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    os.replace(tmp, path)


def plot_task_figures(
    task: pd.DataFrame, nearest: pd.DataFrame, coverage: pd.DataFrame,
    islands: pd.DataFrame, island_internal: list[dict[str, Any]],
    opposite: pd.DataFrame, recommendation: dict[str, Any],
    args: dict[str, Any], task_plot_dir: Path,
) -> None:
    sample = task["sample"].iloc[0]
    region = task["region"].iloc[0]
    a, b = args["reference_celltype"], args["target_celltype"]
    title = f"{sample} | {region}"
    seed = int(hashlib.sha256(f"{args['seed']}|{sample}|{region}".encode()).hexdigest()[:8], 16)

    fig, ax = plt.subplots(figsize=(9, 8))
    display = deterministic_subset(task, args["max_plot_cells"], seed)
    ax.scatter(display["x_um"], display["y_um"], s=0.15, c="#C7C7C7", alpha=0.35, rasterized=True)
    for cell_type in (a, b):
        selected = task.loc[task["cell_type"] == cell_type]
        ax.scatter(selected["x_um"], selected["y_um"], s=2.2, c=CELLTYPE_COLORS.get(cell_type), label=cell_type, rasterized=True)
        eligible = [x for x in island_internal if x["cell_type"] == cell_type and x["eligible_island"]]
        for item in eligible:
            boundaries = [item["polygon"]] if item["polygon"].geom_type == "Polygon" else list(item["polygon"].geoms)
            for polygon in boundaries:
                xy = np.asarray(polygon.exterior.coords)
                ax.plot(xy[:, 0], xy[:, 1], color=CELLTYPE_COLORS.get(cell_type), linewidth=0.7)
        if eligible:
            largest = max(eligible, key=lambda x: x["n_cells"])
            xy = selected[["x_um", "y_um"]].to_numpy()
            center = np.array([largest["centroid_x_um"], largest["centroid_y_um"]])
            representative = xy[np.argmin(np.linalg.norm(xy - center, axis=1))]
            for radius, style in zip((50, 100, 200), ("-", "--", ":"), strict=True):
                ax.add_patch(Circle(representative, radius, fill=False, linestyle=style, linewidth=0.8, color=CELLTYPE_COLORS.get(cell_type)))
    ax.set_aspect("equal")
    ax.set_title(f"Spatial island QC | {title}")
    ax.set_xlabel("Aligned x (µm)")
    ax.set_ylabel("Aligned y (µm)")
    ax.legend(markerscale=4, frameon=False)
    save_figure(fig, task_plot_dir / "spatial_island_qc.png", args["dpi"])

    fig, ax = plt.subplots(figsize=(8, 5.5))
    for (ref, target), group in nearest.groupby(["reference_celltype", "target_celltype"], sort=False):
        values = np.sort(group["nearest_distance_um"].dropna().to_numpy())
        if len(values):
            ax.step(values, np.arange(1, len(values) + 1) / len(values), where="post", label=f"{ref} → {target}")
    for value in args["candidate_distances_um"]:
        ax.axvline(value, color="black", alpha=0.45 if value in (50, 100, 200) else 0.15, linewidth=1.2 if value in (50, 100, 200) else 0.6)
    counts = task["cell_type"].value_counts()
    ax.set(title=f"Directional nearest-neighbor ECDF | {title}\n{a}: {counts.get(a, 0):,}; {b}: {counts.get(b, 0):,}", xlabel="Nearest distance (µm)", ylabel="ECDF", xlim=(0, max(args["candidate_distances_um"]) * 1.15), ylim=(0, 1.01))
    ax.legend(frameon=False, fontsize=8)
    save_figure(fig, task_plot_dir / "nearest_neighbor_ecdf.png", args["dpi"])

    cross = coverage.loc[coverage["reference_celltype"] != coverage["target_celltype"]]
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    for (ref, target), group in cross.groupby(["reference_celltype", "target_celltype"], sort=False):
        label = f"{ref} → {target}"
        axes[0].plot(group["candidate_distance_um"], group["coverage_ge_1"], marker="o", label=label)
        axes[1].plot(group["candidate_distance_um"], group["median_target_neighbors"], marker="o", label=label)
        axes[1].fill_between(group["candidate_distance_um"], group["p25_target_neighbors"], group["p75_target_neighbors"], alpha=0.15)
    axes[0].set(xlabel="Candidate d (µm)", ylabel="Coverage ≥1", ylim=(0, 1.02), title="Directional coverage")
    axes[1].set(xlabel="Candidate d (µm)", ylabel="Target-neighbor count", title="Median and IQR")
    axes[0].legend(frameon=False, fontsize=8)
    fig.suptitle(title)
    save_figure(fig, task_plot_dir / "coverage_and_neighbor_count_curves.png", args["dpi"])

    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8))
    island_values = [islands.loc[(islands["cell_type"] == x) & islands["eligible_island"], "equivalent_diameter_um"].to_numpy() for x in (a, b)]
    axes[0].boxplot(island_values, tick_labels=[a, b], showfliers=True)
    axes[0].set(ylabel="Equivalent diameter (µm)", title="Eligible island sizes")
    if opposite.empty:
        for ax in axes[1:]:
            ax.text(0.5, 0.5, "No eligible opposite-island pairs", ha="center", va="center")
            ax.set_axis_off()
    else:
        axes[1].boxplot([opposite["boundary_distance_um"]], tick_labels=["Boundary"])
        axes[2].boxplot([opposite["centroid_distance_um"]], tick_labels=["Centroid"])
        for ax in axes[1:]:
            for value in (50, 100, 200):
                ax.axhline(value, color="black", alpha=0.2, linewidth=0.8)
            ax.set_ylabel("Distance (µm)")
        axes[1].set_title("Boundary-to-boundary")
        axes[2].set_title("Centroid-to-centroid")
    fig.suptitle(title)
    save_figure(fig, task_plot_dir / "island_size_and_distance_qc.png", args["dpi"])

    metrics = []
    labels = []
    for (ref, target), group in cross.groupby(["reference_celltype", "target_celltype"], sort=False):
        labels.append(f"{ref}→{target} coverage")
        metrics.append(group.sort_values("candidate_distance_um")["coverage_ge_1"].to_numpy())
        labels.append(f"{ref}→{target} median n")
        values = group.sort_values("candidate_distance_um")["median_target_neighbors"].to_numpy()
        metrics.append(np.minimum(values / max(args["median_count_saturation"], 1), 1))
    fig, ax = plt.subplots(figsize=(9, 3.8))
    image = ax.imshow(np.asarray(metrics), aspect="auto", vmin=0, vmax=1, cmap="viridis")
    ax.set_xticks(range(len(args["candidate_distances_um"])), labels=[f"{x:g}" for x in args["candidate_distances_um"]])
    ax.set_yticks(range(len(labels)), labels=labels)
    ax.set_xlabel("Candidate d (µm)")
    ax.set_title(f"Task decision heatmap | {title}\nQC: {recommendation['task_qc_status']}")
    fig.colorbar(image, ax=ax, label="Coverage or count/saturation threshold")
    save_figure(fig, task_plot_dir / "task_decision_heatmap.png", args["dpi"])


def dbscan_sensitivity(task: pd.DataFrame, args: dict[str, Any]) -> pd.DataFrame:
    if not args["run_dbscan_sensitivity"]:
        return pd.DataFrame()
    from sklearn.cluster import DBSCAN

    rows = []
    for cell_type in (args["reference_celltype"], args["target_celltype"]):
        selected = task.loc[task["cell_type"] == cell_type]
        labels = DBSCAN(
            eps=args["dbscan_eps_um"], min_samples=args["dbscan_min_samples"], n_jobs=1
        ).fit_predict(selected[["x_um", "y_um"]])
        for label in sorted(set(labels) - {-1}):
            rows.append({
                "sample": task["sample"].iloc[0], "region": task["region"].iloc[0],
                "cell_type": cell_type, "dbscan_component": int(label),
                "n_cells": int(np.sum(labels == label)),
                "eps_um": args["dbscan_eps_um"],
                "min_samples": args["dbscan_min_samples"],
            })
    return pd.DataFrame(rows)


def analyze_task(payload: tuple[pd.DataFrame, dict[str, Any], str]) -> dict[str, Any]:
    task, args, output_string = payload
    output_dir = Path(output_string)
    sample = task["sample"].iloc[0]
    region = task["region"].iloc[0]
    key = task_key(sample, region)
    task_dir = output_dir / "tasks" / region / sample
    plot_dir = Path(args["plot_dir"]) / "per_task" / region / sample
    task_dir.mkdir(parents=True, exist_ok=True)
    expected_hash = task_hash(args, sample, region)
    complete_path = task_dir / "task_status.json"
    if args["resume"] and complete_path.exists():
        old = json.loads(complete_path.read_text())
        if old.get("parameter_hash") != expected_hash:
            raise RuntimeError(f"{key}: existing task parameters differ; use --overwrite")
        if old.get("status") == "completed":
            return {"task_key": key, "status": "resumed", "task_dir": str(task_dir)}

    a, b = args["reference_celltype"], args["target_celltype"]
    counts = task["cell_type"].value_counts()
    n_a, n_b = int(counts.get(a, 0)), int(counts.get(b, 0))
    eligibility = {
        "sample": sample, "region": region, "n_total_cells": len(task),
        "reference_celltype": a, "n_reference_cells": n_a,
        "target_celltype": b, "n_target_cells": n_b,
        "min_cells_per_type": args["min_cells_per_type"],
        "eligible_cell_counts": n_a >= args["min_cells_per_type"] and n_b >= args["min_cells_per_type"],
        "eligibility_reason": "eligible",
    }
    if not eligibility["eligible_cell_counts"]:
        eligibility["eligibility_reason"] = "insufficient_reference_cells" if n_a < args["min_cells_per_type"] else "insufficient_target_cells"
        atomic_csv(pd.DataFrame([eligibility]), task_dir / "task_eligibility.csv")
        atomic_parquet(empty_frame(NN_COLUMNS), task_dir / "nearest_neighbor_distances.parquet")
        atomic_csv(empty_frame(COVERAGE_COLUMNS), task_dir / "neighborhood_coverage_and_counts.csv")
        atomic_csv(empty_frame(ISLAND_COLUMNS), task_dir / "island_components_all.csv")
        atomic_csv(empty_frame(ISLAND_COLUMNS), task_dir / "island_metrics.csv")
        atomic_csv(empty_frame(OPPOSITE_COLUMNS), task_dir / "nearest_opposite_island_distances.csv")
        atomic_csv(pd.DataFrame([{
            "sample": sample, "region": region,
            "recommendation_status": "insufficient_data",
            "task_qc_status": "insufficient_data",
            "recommendation_reason": (
                "Too few A or B cells for neighborhood-distance diagnosis."
            ),
            "coverage_elbow_um": np.nan,
            "median_nearest_opposite_island_boundary_gap_um": np.nan,
            "maximum_tested_candidate_distance_um": max(args["candidate_distances_um"]),
            "nearest_island_median_gap_candidate_um": np.nan,
            "support_threshold_um": np.nan,
            "recommended_distance_um": np.nan,
            "recommended_task_distance_um": np.nan,
            "maximum_directional_coverage_difference": np.nan,
            "strong_directionality": False,
        }]), task_dir / "task_level_recommendation.csv")
        atomic_json({"status": "completed", "parameter_hash": expected_hash, "task_qc_status": "insufficient_data"}, complete_path)
        return {"task_key": key, "status": "completed", "task_dir": str(task_dir)}

    nearest, coverage = coverage_summary(task, args)
    bin_size = args["island_bin_size_um"]
    origin = (
        math.floor(task["x_um"].min() / bin_size) * bin_size,
        math.floor(task["y_um"].min() / bin_size) * bin_size,
    )
    all_ix = np.floor((task["x_um"].to_numpy() - origin[0]) / bin_size).astype(int)
    all_iy = np.floor((task["y_um"].to_numpy() - origin[1]) / bin_size).astype(int)
    shape = (int(all_iy.max()) + 1, int(all_ix.max()) + 1)
    roi_mask = np.zeros(shape, dtype=bool)
    roi_mask[all_iy, all_ix] = True
    island_frames = []
    internals: dict[str, list[dict[str, Any]]] = {}
    closing_qc = {}
    for cell_type in (a, b):
        frame, internal, qc = build_islands(task, cell_type, args, origin, shape, roi_mask)
        island_frames.append(frame)
        internals[cell_type] = internal
        closing_qc[cell_type] = qc
    islands = pd.concat(island_frames, ignore_index=True)
    opposite = pd.concat([
        opposite_island_distances(internals[a], internals[b], a, b, sample, region),
        opposite_island_distances(internals[b], internals[a], b, a, sample, region),
    ], ignore_index=True)
    n_eligible_a = int(((islands["cell_type"] == a) & islands["eligible_island"]).sum())
    n_eligible_b = int(((islands["cell_type"] == b) & islands["eligible_island"]).sum())
    eligibility.update({
        "n_eligible_reference_islands": n_eligible_a,
        "n_eligible_target_islands": n_eligible_b,
        "eligible_islands": n_eligible_a > 0 and n_eligible_b > 0,
        "reference_bins_added_by_closing": closing_qc[a]["bins_added_by_closing"],
        "target_bins_added_by_closing": closing_qc[b]["bins_added_by_closing"],
        "reference_components_connected_by_closing": closing_qc[a]["components_connected_by_closing"],
        "target_components_connected_by_closing": closing_qc[b]["components_connected_by_closing"],
        "reference_cells_connected_by_closing": closing_qc[a]["cells_connected_by_closing"],
        "target_cells_connected_by_closing": closing_qc[b]["cells_connected_by_closing"],
    })
    if not eligibility["eligible_islands"]:
        eligibility["eligibility_reason"] = "insufficient_eligible_islands"
    recommendation = recommend_task(coverage, opposite, args, sample, region)
    if not eligibility["eligible_islands"]:
        recommendation["recommendation_status"] = "insufficient_data"
        recommendation["task_qc_status"] = "insufficient_data"
        recommendation["recommendation_reason"] = (
            "Too few eligible A or B islands for island-gap diagnosis."
        )
        recommendation["recommended_distance_um"] = None
        recommendation["recommended_task_distance_um"] = None

    atomic_csv(pd.DataFrame([eligibility]), task_dir / "task_eligibility.csv")
    atomic_parquet(nearest, task_dir / "nearest_neighbor_distances.parquet")
    atomic_csv(coverage, task_dir / "neighborhood_coverage_and_counts.csv")
    atomic_csv(islands, task_dir / "island_components_all.csv")
    atomic_csv(islands.loc[islands["eligible_island"]], task_dir / "island_metrics.csv")
    atomic_csv(opposite, task_dir / "nearest_opposite_island_distances.csv")
    atomic_csv(pd.DataFrame([recommendation]), task_dir / "task_level_recommendation.csv")
    dbscan = dbscan_sensitivity(task, args)
    if args["run_dbscan_sensitivity"]:
        atomic_csv(dbscan, task_dir / "dbscan_components_all.csv")
    plot_task_figures(
        task, nearest, coverage, islands,
        internals[a] + internals[b], opposite, recommendation, args, plot_dir,
    )
    atomic_json({
        "status": "completed", "parameter_hash": expected_hash,
        "task_qc_status": recommendation["task_qc_status"],
        "coordinate_unit": "micrometer",
    }, complete_path)
    return {"task_key": key, "status": "completed", "task_dir": str(task_dir)}


def combine_task_outputs(output_dir: Path, tasks: list[tuple[str, str]]) -> dict[str, pd.DataFrame]:
    specs = {
        "task_eligibility": ("task_eligibility.csv", "csv"),
        "nearest_neighbor_distances": ("nearest_neighbor_distances.parquet", "parquet"),
        "neighborhood_coverage_and_counts": ("neighborhood_coverage_and_counts.csv", "csv"),
        "island_components_all": ("island_components_all.csv", "csv"),
        "island_metrics": ("island_metrics.csv", "csv"),
        "nearest_opposite_island_distances": ("nearest_opposite_island_distances.csv", "csv"),
        "task_level_recommendations": ("task_level_recommendation.csv", "csv"),
    }
    combined: dict[str, pd.DataFrame] = {}
    for output_name, (filename, kind) in specs.items():
        frames = []
        for sample, region in tasks:
            path = output_dir / "tasks" / region / sample / filename
            if not path.exists():
                raise FileNotFoundError(f"Missing completed task output: {path}")
            frames.append(pd.read_parquet(path) if kind == "parquet" else pd.read_csv(path))
        nonempty = [frame for frame in frames if not frame.empty]
        if nonempty:
            combined[output_name] = pd.concat(nonempty, ignore_index=True)
        elif frames:
            combined[output_name] = frames[0].copy()
        else:
            combined[output_name] = pd.DataFrame()
        target = output_dir / (f"{output_name}.parquet" if kind == "parquet" else f"{output_name}.csv")
        if kind == "parquet":
            atomic_parquet(combined[output_name], target)
        else:
            atomic_csv(combined[output_name], target)
    return combined


def aggregate_candidate_summary(
    eligibility: pd.DataFrame, coverage: pd.DataFrame,
    recommendations: pd.DataFrame, args: argparse.Namespace,
) -> pd.DataFrame:
    eligible_keys = eligibility.loc[
        eligibility["eligible_cell_counts"].astype(bool) & eligibility.get("eligible_islands", False).fillna(False),
        ["sample", "region"],
    ].drop_duplicates()
    eligible = coverage.merge(eligible_keys, on=["sample", "region"], how="inner")
    cross = eligible.loc[eligible["reference_celltype"] != eligible["target_celltype"]].copy()
    status_column = (
        "recommendation_status"
        if "recommendation_status" in recommendations
        else "task_qc_status"
    )
    distance_column = (
        "recommended_distance_um"
        if "recommended_distance_um" in recommendations
        else "recommended_task_distance_um"
    )
    rows = []
    groups: list[tuple[str, pd.DataFrame]] = [("overall", eligible_keys)]
    groups.extend((region, group[["sample", "region"]]) for region, group in eligible_keys.groupby("region"))
    for scope, keys in groups:
        rec = recommendations.merge(keys, on=["sample", "region"], how="inner")
        excluded = rec.loc[
            rec[status_column] == "island_gap_exceeds_candidate_range",
            ["sample", "region"],
        ].drop_duplicates()
        recommendation_keys = keys.merge(
            excluded.assign(_excluded=True),
            on=["sample", "region"], how="left",
        ).loc[lambda x: x["_excluded"].isna(), ["sample", "region"]]
        rec_for_global = rec.merge(
            recommendation_keys, on=["sample", "region"], how="inner"
        )
        cov = cross.merge(
            recommendation_keys, on=["sample", "region"], how="inner"
        )
        denominator = len(recommendation_keys)
        for candidate in args.candidate_distances_um:
            at_d = cov.loc[np.isclose(cov["candidate_distance_um"], candidate)].copy()
            saturation = at_d.assign(
                saturated=(at_d["coverage_ge_1"] >= args.coverage_saturation_threshold)
                | (at_d["median_target_neighbors"] >= args.median_count_saturation)
                | (at_d["p90_target_neighbors"] >= args.p90_count_saturation)
            ).groupby(["sample", "region"])["saturated"].any()
            n_recommending = int(np.isclose(
                pd.to_numeric(rec_for_global[distance_column], errors="coerce"),
                candidate,
            ).sum())
            rows.append({
                "scope": scope, "candidate_distance_um": candidate,
                "n_island_eligible_tasks": len(keys),
                "n_excluded_island_gap_exceeds_candidate_range": len(excluded),
                "n_eligible_tasks": denominator,
                "n_tasks_recommending_candidate": n_recommending,
                "fraction_tasks_recommending_candidate": n_recommending / denominator if denominator else np.nan,
                "n_saturated_tasks": int(saturation.sum()),
                "fraction_saturated_tasks": float(saturation.mean()) if len(saturation) else np.nan,
            })
    return pd.DataFrame(rows)


def plot_aggregate(
    combined: dict[str, pd.DataFrame], summary: pd.DataFrame,
    args: argparse.Namespace, plot_root: Path,
) -> None:
    figure_dir = plot_root / "aggregate"
    coverage = combined["neighborhood_coverage_and_counts"]
    eligibility = combined["task_eligibility"]
    keys = eligibility.loc[eligibility["eligible_cell_counts"].astype(bool), ["sample", "region"]]
    cross = coverage.merge(keys, on=["sample", "region"]).loc[lambda x: x.reference_celltype != x.target_celltype]
    direction_column = cross["reference_celltype"] + " → " + cross["target_celltype"]
    cross = cross.assign(direction=direction_column)

    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for direction, group in cross.groupby("direction"):
        for _, task in group.groupby(["sample", "region"]):
            axes[0].plot(task["candidate_distance_um"], task["coverage_ge_1"], color="#777777", alpha=0.12, linewidth=0.7)
            axes[1].plot(task["candidate_distance_um"], task["median_target_neighbors"], color="#777777", alpha=0.12, linewidth=0.7)
        stats = group.groupby("candidate_distance_um").agg(
            coverage_median=("coverage_ge_1", "median"),
            coverage_q25=("coverage_ge_1", lambda x: x.quantile(.25)),
            coverage_q75=("coverage_ge_1", lambda x: x.quantile(.75)),
            count_median=("median_target_neighbors", "median"),
            count_q25=("median_target_neighbors", lambda x: x.quantile(.25)),
            count_q75=("median_target_neighbors", lambda x: x.quantile(.75)),
        ).reset_index()
        axes[0].plot(stats["candidate_distance_um"], stats["coverage_median"], marker="o", linewidth=2, label=direction)
        axes[0].fill_between(stats["candidate_distance_um"], stats["coverage_q25"], stats["coverage_q75"], alpha=.18)
        axes[1].plot(stats["candidate_distance_um"], stats["count_median"], marker="o", linewidth=2, label=direction)
        axes[1].fill_between(stats["candidate_distance_um"], stats["count_q25"], stats["count_q75"], alpha=.18)
    axes[0].set(xlabel="Candidate d (µm)", ylabel="Coverage ≥1", ylim=(0, 1.02), title="Coverage: tasks, median, and IQR")
    axes[1].set(xlabel="Candidate d (µm)", ylabel="Median target-neighbor count", title="Neighbor counts: tasks, median, and IQR")
    axes[0].legend(frameon=False)
    save_figure(fig, figure_dir / "aggregate_coverage_and_neighbor_counts.png", args.dpi)

    islands = combined["island_metrics"]
    opposite = combined["nearest_opposite_island_distances"]
    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))
    island_groups = [group["equivalent_diameter_um"].to_numpy() for _, group in islands.groupby(["region", "cell_type"])]
    island_labels = [f"{r}\n{c}" for (r, c), _ in islands.groupby(["region", "cell_type"])]
    if island_groups:
        axes[0].boxplot(island_groups, tick_labels=island_labels, showfliers=False)
        sample_points = islands.groupby(
            ["region", "cell_type", "sample"], as_index=False
        )["equivalent_diameter_um"].median()
        group_keys = list(islands.groupby(["region", "cell_type"]).groups)
        key_to_x = {key: index + 1 for index, key in enumerate(group_keys)}
        point_x = [key_to_x[(row.region, row.cell_type)] for row in sample_points.itertuples()]
        axes[0].scatter(
            point_x, sample_points["equivalent_diameter_um"],
            s=11, c="black", alpha=.55, zorder=3, label="Sample median",
        )
        axes[0].tick_params(axis="x", rotation=35)
        axes[0].legend(frameon=False, fontsize=8)
    axes[0].set(ylabel="Equivalent diameter (µm)", title="Island-size distributions by region/type")
    if not opposite.empty:
        axes[1].hist(opposite["boundary_distance_um"], bins=30, alpha=.7)
        axes[2].hist(opposite["centroid_distance_um"], bins=30, alpha=.7)
        for ax in axes[1:]:
            for value in (50, 100, 200):
                ax.axvline(value, color="black", alpha=.35, linestyle="--")
            ax.set_xlabel("Distance (µm)")
            ax.set_ylabel("Nearest-island observations")
    axes[1].set_title("Boundary-to-boundary distribution")
    axes[2].set_title("Centroid-to-centroid distribution")
    save_figure(fig, figure_dir / "aggregate_island_sizes_and_distances.png", args.dpi)

    task_order = [f"{region} | {sample}" for region in args.regions for sample in args.samples]
    cross = cross.assign(task=cross["region"] + " | " + cross["sample"])
    fig, axes = plt.subplots(3, 2, figsize=(13, max(13, len(task_order) * .28)))
    metrics = [
        (args.reference_celltype, "coverage_ge_1", f"{args.reference_celltype} → {args.target_celltype} coverage"),
        (args.target_celltype, "coverage_ge_1", f"{args.target_celltype} → {args.reference_celltype} coverage"),
        (args.reference_celltype, "median_target_neighbors", f"{args.reference_celltype} → {args.target_celltype} median n"),
        (args.target_celltype, "median_target_neighbors", f"{args.target_celltype} → {args.reference_celltype} median n"),
    ]
    for ax, (reference, value, title) in zip(axes.flat[:4], metrics, strict=True):
        matrix = cross.loc[cross["reference_celltype"] == reference].pivot(index="task", columns="candidate_distance_um", values=value).reindex(task_order)
        image = ax.imshow(matrix.to_numpy(), aspect="auto", cmap="viridis", vmin=0 if value == "coverage_ge_1" else None, vmax=1 if value == "coverage_ge_1" else None)
        ax.set_xticks(range(len(matrix.columns)), labels=[f"{x:g}" for x in matrix.columns])
        ax.set_yticks(range(len(matrix.index)), labels=matrix.index, fontsize=6)
        ax.set_title(title)
        ax.set_xlabel("Candidate d (µm)")
        fig.colorbar(image, ax=ax, fraction=.035)
    status_order = [
        "recommended", "no_stable_support", "ambiguous_large_d",
        "island_gap_exceeds_candidate_range", "insufficient_data",
    ]
    status_map = {status: index for index, status in enumerate(status_order)}
    recommendations = combined["task_level_recommendations"].copy()
    recommendations["task"] = recommendations["region"] + " | " + recommendations["sample"]
    status_column = (
        "recommendation_status"
        if "recommendation_status" in recommendations
        else "task_qc_status"
    )
    status_values = recommendations.set_index("task")[status_column].reindex(task_order)
    status_numeric = status_values.map(status_map).fillna(len(status_order)).to_numpy()[:, None]
    status_ax = axes.flat[4]
    status_image = status_ax.imshow(status_numeric, aspect="auto", cmap="tab10", vmin=0, vmax=9)
    status_ax.set_yticks(range(len(task_order)), labels=task_order, fontsize=6)
    status_ax.set_xticks([0], labels=["QC status"])
    status_ax.set_title("Task-level QC status")
    legend_text = "\n".join(f"{index}: {status}" for status, index in status_map.items())
    status_ax.text(1.15, .5, legend_text, transform=status_ax.transAxes, va="center", fontsize=8)
    axes.flat[5].set_axis_off()
    fig.suptitle("QC decision heatmap across all region × sample tasks")
    save_figure(fig, figure_dir / "qc_decision_heatmap.png", args.dpi)


def write_report(summary: pd.DataFrame, recommendations: pd.DataFrame, args: argparse.Namespace, output_dir: Path) -> None:
    overall = summary.loc[summary["scope"] == "overall"].sort_values("candidate_distance_um")
    supported = overall.loc[
        (overall["fraction_tasks_recommending_candidate"] >= args.global_support_fraction)
        & (overall["fraction_saturated_tasks"] <= args.global_max_saturation_fraction)
    ]
    global_d = float(supported["candidate_distance_um"].min()) if not supported.empty else None
    status_column = (
        "recommendation_status"
        if "recommendation_status" in recommendations
        else "task_qc_status"
    )
    unresolved = recommendations.loc[
        recommendations[status_column] == "island_gap_exceeds_candidate_range"
    ].copy()
    n_excluded = len(unresolved)
    lines = [
        "# CRAWDAD neighborhood-distance diagnostic recommendation",
        "",
        f"Reference: `{args.reference_celltype}`; target: `{args.target_celltype}`.",
        "",
        "This is a permutation-free geometry/QC diagnostic. Coverage is descriptive and can increase because of spatial proximity, target-cell abundance, or both. It is not a significance test and it does not diagnose CRAWDAD null shuffle scale `r`.",
        "",
        "## Automatic advisory conclusion",
        "",
    ]
    if global_d is None:
        lines.append(
            "No single candidate met the pragmatic global rules. Use a prespecified sensitivity analysis rather than selecting a separate distance from each region or slice."
        )
    else:
        lines.append(
            f"The automatic rule supports **d = {global_d:g} µm** as the global primary candidate. Confirm this against the QC figures before changing downstream CRAWDAD settings."
        )
    lines.extend([
        "",
        "The global rule requires at least "
        f"{args.global_support_fraction:.0%} of eligible tasks to recommend the candidate and no more than "
        f"{args.global_max_saturation_fraction:.0%} to be saturated. These are pragmatic QC thresholds, not statistical significance thresholds.",
        "",
        f"Tasks excluded from the global recommendation because the median island gap exceeds the tested range: **{n_excluded}**.",
        "",
        "## Candidate summary",
        "",
        "| d (µm) | island-eligible tasks | excluded: gap out of range | global denominator | recommending | saturated |",
        "|---:|---:|---:|---:|---:|---:|",
    ])
    for _, row in overall.iterrows():
        recommending = (
            f"{row.fraction_tasks_recommending_candidate:.1%}"
            if np.isfinite(row.fraction_tasks_recommending_candidate)
            else "NA"
        )
        saturated = (
            f"{row.fraction_saturated_tasks:.1%}"
            if np.isfinite(row.fraction_saturated_tasks)
            else "NA"
        )
        lines.append(
            f"| {row.candidate_distance_um:g} | {int(row.n_island_eligible_tasks)} | {int(row.n_excluded_island_gap_exceeds_candidate_range)} | {int(row.n_eligible_tasks)} | {recommending} | {saturated} |"
        )
    if n_excluded:
        lines.extend([
            "", "## Unresolved island gaps beyond the tested range", "",
            ISLAND_GAP_OUT_OF_RANGE_REASON,
            "",
            "| sample | region | median boundary gap (µm) | maximum tested d (µm) |",
            "|---|---|---:|---:|",
        ])
        for _, row in unresolved.sort_values(["region", "sample"]).iterrows():
            lines.append(
                f"| {row['sample']} | {row['region']} | "
                f"{row['median_nearest_opposite_island_boundary_gap_um']:.3g} | "
                f"{row['maximum_tested_candidate_distance_um']:.3g} |"
            )
    lines.extend([
        "", "## Interpretation guide", "",
        "- 50 µm corresponds most directly to local cell/boundary contact.",
        "- 100 µm may better capture intermediate A/B island proximity.",
        "- 200 µm is useful only when it reaches neighboring islands without widespread coverage/count saturation.",
        "- If 200 µm saturates coverage or counts, it is too broad for a spatially specific neighborhood.",
        "- If tasks disagree strongly, no single d is adequate across all regions and slices.",
        "", "## Task QC statuses", "",
    ])
    status_counts = recommendations[status_column].value_counts(dropna=False)
    for status, count in status_counts.items():
        lines.append(f"- `{status}`: {count}")
    (output_dir / "neighborhood_distance_recommendation.md").write_text("\n".join(lines) + "\n")


def args_for_worker(args: argparse.Namespace) -> dict[str, Any]:
    value = vars(args).copy()
    for key, item in list(value.items()):
        if isinstance(item, Path):
            value[key] = str(item)
    return value


def run_self_tests() -> None:
    # Touching and 100-µm-separated raster polygons.
    touching_a = box(0, 0, 25, 25)
    touching_b = box(25, 0, 50, 25)
    assert touching_a.distance(touching_b) == 0
    separated = box(125, 0, 150, 25)
    assert abs(touching_a.distance(separated) - 100) <= 25

    # Directional coverage differs under asymmetric geometry/abundance.
    a_xy = np.array([[0.0, 0.0], [1000.0, 0.0]])
    b_xy = np.array([[10.0, 0.0]])
    a_to_b = np.mean(radius_counts(a_xy, b_xy, 50, False) >= 1)
    b_to_a = np.mean(radius_counts(b_xy, a_xy, 50, False) >= 1)
    assert a_to_b != b_to_a

    # Self-exclusion returns positive same-type distances for unique cells.
    same = np.array([[0.0, 0.0], [10.0, 0.0], [30.0, 0.0]])
    assert np.all(nearest_distances(same, same, True) > 0)

    # Island-gap recommendation boundary cases.
    recommendation_args = {
        "reference_celltype": "A", "target_celltype": "B",
        "candidate_distances_um": [50.0, 100.0, 300.0],
        "elbow_min_score": .05, "coverage_saturation_threshold": .9,
        "median_count_saturation": 10.0, "p90_count_saturation": 50.0,
        "directionality_threshold": .3,
    }
    recommendation_coverage_rows = []
    for reference, target in (("A", "B"), ("B", "A")):
        for candidate, value in zip(
            recommendation_args["candidate_distances_um"],
            (.10, .55, .80), strict=True,
        ):
            recommendation_coverage_rows.append({
                "sample": "s1", "region": "lateral",
                "reference_celltype": reference, "target_celltype": target,
                "candidate_distance_um": candidate,
                "coverage_ge_1": value,
                "median_target_neighbors": 1.0,
                "p90_target_neighbors": 2.0,
            })
    recommendation_coverage = pd.DataFrame(recommendation_coverage_rows)

    def recommendation_for_gap(gap: float | None) -> dict[str, Any]:
        if gap is None:
            opposite = pd.DataFrame(columns=["boundary_distance_um"])
        else:
            opposite = pd.DataFrame({"boundary_distance_um": [gap, gap]})
        return recommend_task(
            recommendation_coverage, opposite, recommendation_args,
            "s1", "lateral",
        )

    within = recommendation_for_gap(120.0)
    assert within["recommendation_status"] == "recommended"
    assert within["median_nearest_opposite_island_boundary_gap_um"] == 120.0
    assert within["nearest_island_median_gap_candidate_um"] == 300.0

    exactly_maximum = recommendation_for_gap(300.0)
    assert exactly_maximum["recommendation_status"] == "recommended"
    assert exactly_maximum["recommended_distance_um"] == 300.0

    beyond = recommendation_for_gap(300.01)
    assert beyond["recommendation_status"] == "island_gap_exceeds_candidate_range"
    assert beyond["recommended_distance_um"] is None
    assert beyond["maximum_tested_candidate_distance_um"] == 300.0
    assert beyond["recommendation_reason"] == ISLAND_GAP_OUT_OF_RANGE_REASON

    missing_gap = recommendation_for_gap(None)
    assert missing_gap["recommendation_status"] == "insufficient_data"
    nonfinite_gap = recommendation_for_gap(np.nan)
    assert nonfinite_gap["recommendation_status"] == "insufficient_data"

    # Out-of-range tasks remain in QC counts but leave the global denominator.
    aggregate_eligibility = pd.DataFrame({
        "sample": ["s1", "s2"], "region": ["lateral", "lateral"],
        "eligible_cell_counts": [True, True],
        "eligible_islands": [True, True],
    })
    aggregate_coverage = pd.concat([
        recommendation_coverage.assign(sample="s1"),
        recommendation_coverage.assign(sample="s2"),
    ], ignore_index=True)
    aggregate_recommendations = pd.DataFrame({
        "sample": ["s1", "s2"], "region": ["lateral", "lateral"],
        "recommendation_status": [
            "recommended", "island_gap_exceeds_candidate_range"
        ],
        "recommended_distance_um": [100.0, np.nan],
    })
    aggregate_args = argparse.Namespace(
        candidate_distances_um=[50.0, 100.0, 300.0],
        coverage_saturation_threshold=.9,
        median_count_saturation=10.0,
        p90_count_saturation=50.0,
    )
    aggregate_test = aggregate_candidate_summary(
        aggregate_eligibility, aggregate_coverage,
        aggregate_recommendations, aggregate_args,
    )
    aggregate_row = aggregate_test.loc[
        (aggregate_test["scope"] == "overall")
        & np.isclose(aggregate_test["candidate_distance_um"], 100.0)
    ].iloc[0]
    assert aggregate_row["n_island_eligible_tasks"] == 2
    assert aggregate_row["n_excluded_island_gap_exceeds_candidate_range"] == 1
    assert aggregate_row["n_eligible_tasks"] == 1
    assert aggregate_row["fraction_tasks_recommending_candidate"] == 1.0

    # Repeated calls are deterministic and retain ineligible-task reasoning.
    synthetic = pd.DataFrame({
        "cell_id": ["a1", "a2", "b1"], "sample": "synthetic",
        "region": "lateral", "cell_type": ["A", "A", "B"],
        "x_um": [0.0, 1000.0, 10.0], "y_um": [0.0, 0.0, 0.0],
    })
    test_args = {
        "reference_celltype": "A", "target_celltype": "B",
        "candidate_distances_um": [25.0, 50.0, 100.0],
    }
    first = coverage_summary(synthetic, test_args)[1]
    second = coverage_summary(synthetic, test_args)[1]
    pd.testing.assert_frame_equal(first, second)
    assert len(synthetic.loc[synthetic.cell_type == "B"]) < 20

    ineligible_args = {
        **test_args,
        "island_bin_size_um": 25.0, "island_connectivity": 8,
        "island_closing_bins": 1, "min_island_cells": 20,
        "min_cells_per_type": 20, "coverage_saturation_threshold": .9,
        "median_count_saturation": 10.0, "p90_count_saturation": 50.0,
        "directionality_threshold": .3, "elbow_min_score": .05,
        "run_dbscan_sensitivity": False, "dbscan_eps_um": 50.0,
        "dbscan_min_samples": 20, "seed": 1, "resume": False,
        "max_plot_cells": 1000, "dpi": 100, "input_identity": "synthetic",
    }
    with tempfile.TemporaryDirectory() as temporary:
        ineligible_args["plot_dir"] = temporary
        result = analyze_task((synthetic, ineligible_args, temporary))
        qc = pd.read_csv(
            Path(result["task_dir"]) / "task_eligibility.csv"
        )
        assert qc.loc[0, "eligibility_reason"] in {
            "insufficient_reference_cells", "insufficient_target_cells"
        }
        assert not bool(qc.loc[0, "eligible_cell_counts"])

    # Numerical task operations are worker-count independent because every
    # task uses exact single-threaded spatial queries and sorted aggregation.
    with concurrent.futures.ProcessPoolExecutor(max_workers=2) as executor:
        parallel = list(executor.map(_self_test_numeric_worker, [synthetic, synthetic]))
    assert parallel[0] == parallel[1]
    print("All Module 01 synthetic tests passed (distances in micrometers).")


def _self_test_numeric_worker(frame: pd.DataFrame) -> str:
    args = {
        "reference_celltype": "A", "target_celltype": "B",
        "candidate_distances_um": [25.0, 50.0, 100.0],
    }
    result = coverage_summary(frame, args)[1]
    return hashlib.sha256(result.to_csv(index=False).encode()).hexdigest()


def parameter_manifest(
    args: argparse.Namespace, regions: Sequence[str], samples: Sequence[str]
) -> dict[str, Any]:
    values = {
        key: (str(value.resolve()) if isinstance(value, Path) else value)
        for key, value in vars(args).items()
    }
    values["regions"] = list(regions)
    values["samples"] = list(samples)
    values.update({
        "canonical_samples": list(CANONICAL_SAMPLES),
        "coordinate_unit_resolved": "micrometer",
        "source_files": list(getattr(args, "source_files", [])),
        "input_identity": getattr(args, "input_identity", []),
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "scientific_scope": (
            "neighborhood distance d only; no permutations and no "
            "null-scale diagnosis"
        ),
        "recommendation_logic_version": RECOMMENDATION_LOGIC_VERSION,
    })
    return values


def finalize_outputs(
    output_dir: Path, tasks: list[tuple[str, str]], args: argparse.Namespace
) -> None:
    combined = combine_task_outputs(output_dir, tasks)
    summary = aggregate_candidate_summary(
        combined["task_eligibility"],
        combined["neighborhood_coverage_and_counts"],
        combined["task_level_recommendations"], args,
    )
    atomic_csv(summary, output_dir / "aggregate_candidate_summary.csv")
    plot_aggregate(combined, summary, args, args.plot_dir.resolve())
    write_report(
        summary, combined["task_level_recommendations"], args, output_dir
    )


def try_finalize_global(
    output_dir: Path, args: argparse.Namespace
) -> bool:
    """Serialize global aggregation and run it only after 44 compatible tasks."""
    lock_path = output_dir / ".global_aggregation.lock"
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with lock_path.open("a+") as lock_handle:
        fcntl.flock(lock_handle, fcntl.LOCK_EX)
        aggregate_args = copy.copy(args)
        aggregate_args.regions = list(DEFAULT_REGIONS)
        aggregate_args.samples = list(CANONICAL_SAMPLES)
        expected_args = args_for_worker(aggregate_args)
        tasks = [
            (sample, region)
            for region in DEFAULT_REGIONS
            for sample in CANONICAL_SAMPLES
        ]
        incomplete: list[str] = []
        incompatible: list[str] = []
        for sample, region in tasks:
            status_path = (
                output_dir / "tasks" / region / sample / "task_status.json"
            )
            if not status_path.exists():
                incomplete.append(task_key(sample, region))
                continue
            status = json.loads(status_path.read_text())
            if status.get("status") != "completed":
                incomplete.append(task_key(sample, region))
                continue
            expected = task_hash(expected_args, sample, region)
            if status.get("parameter_hash") != expected:
                incompatible.append(task_key(sample, region))
        if incomplete or incompatible:
            print(
                "Global aggregation deferred: "
                f"{len(incomplete)} incomplete and "
                f"{len(incompatible)} parameter-incompatible tasks.",
                flush=True,
            )
            return False
        atomic_json(
            parameter_manifest(aggregate_args, DEFAULT_REGIONS, CANONICAL_SAMPLES),
            output_dir / "parameters.json",
        )
        finalize_outputs(output_dir, tasks, aggregate_args)
        atomic_json({
            "status": "completed",
            "n_tasks": len(tasks),
            "completed_at_utc": datetime.now(timezone.utc).isoformat(),
        }, output_dir / "global_aggregation_status.json")
        print("Global 44-task aggregation completed.", flush=True)
        return True


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    if args.self_test:
        run_self_tests()
        return 0
    output_dir = args.output_dir.resolve()
    plot_root = args.plot_dir.resolve()
    if args.overwrite and not args.regional_job:
        if output_dir.exists():
            shutil.rmtree(output_dir)
        if plot_root.exists():
            shutil.rmtree(plot_root)
    elif output_dir.exists() and not args.resume:
        raise FileExistsError(f"Output exists; use --resume or --overwrite: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=True)
    (plot_root / "per_task").mkdir(parents=True, exist_ok=True)
    (plot_root / "aggregate").mkdir(parents=True, exist_ok=True)

    cells, sources = load_cells(args)
    input_identity = []
    for source in sources:
        source_path = Path(source)
        stat = source_path.stat()
        input_identity.append({
            "path": str(source_path), "size_bytes": stat.st_size,
            "mtime_ns": stat.st_mtime_ns,
        })
    args.input_identity = input_identity
    args.source_files = sources
    parameters = parameter_manifest(args, args.regions, args.samples)
    if args.regional_job:
        manifest_name = f"{args.regions[0]}.json"
        atomic_json(parameters, output_dir / "run_manifests" / manifest_name)
    else:
        atomic_json(parameters, output_dir / "parameters.json")

    tasks: list[tuple[str, str]] = [
        (sample, region) for region in args.regions for sample in args.samples
    ]
    grouped = {
        (sample, region): group.copy()
        for (sample, region), group in cells.groupby(["sample", "region"], sort=False)
    }
    missing = [x for x in tasks if x not in grouped]
    if missing:
        raise ValueError(f"Requested sample-region tasks have no cells: {missing}")
    if args.regional_job and args.overwrite:
        for sample, region in tasks:
            task_dir = output_dir / "tasks" / region / sample
            plot_dir = plot_root / "per_task" / region / sample
            if task_dir.exists():
                shutil.rmtree(task_dir)
            if plot_dir.exists():
                shutil.rmtree(plot_dir)
    worker_args = args_for_worker(args)
    payloads = [(grouped[key], worker_args, str(output_dir)) for key in tasks]
    print(f"Running {len(tasks)} region × sample tasks with n_jobs={min(args.n_jobs, len(tasks))}", flush=True)
    if args.n_jobs == 1:
        results = [analyze_task(payload) for payload in payloads]
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=min(args.n_jobs, len(tasks))
        ) as executor:
            results = list(executor.map(analyze_task, payloads))
    print(pd.DataFrame(results)["status"].value_counts().to_string(), flush=True)

    global_completed = False
    if args.regional_job:
        global_completed = try_finalize_global(output_dir, args)
    else:
        finalize_outputs(output_dir, tasks, args)
        global_completed = True
    print(f"Output: {output_dir}")
    print(f"Plots: {plot_root}")
    if global_completed:
        print(f"Recommendation: {output_dir / 'neighborhood_distance_recommendation.md'}")
    else:
        print("Recommendation: pending completion of all four regional jobs")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise
