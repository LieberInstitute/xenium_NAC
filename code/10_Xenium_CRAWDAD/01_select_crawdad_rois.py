#!/usr/bin/env python3
"""Define three non-overlapping anatomical ROIs for Xenium CRAWDAD.

The script uses final second-pass Xenium coordinates and six final second-pass
Visium HD coordinate tables that already share the aligned Xenium reference
frame. It constructs six VHD footprints, defines lateral/dorsomedial/
ventromedial primary ROIs, assigns every Xenium cell, and exports per-slice
CRAWDAD input tables and QC deliverables.

The global dorsal-ventral split is the midpoint of the robust 0.05th-to-99.95th
percentile aligned-y range of broad MSN Xenium cells inside the medial ROI.
Larger y is dorsomedial and smaller y is ventromedial.
"""

from __future__ import annotations

import argparse
import importlib.metadata
import json
import logging
import os
import re
import subprocess
import sys
from collections.abc import Iterable, Sequence
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_crawdad_rois")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Polygon as MplPolygon
from shapely import (
    MultiPoint,
    concave_hull,
    distance,
    intersects_xy,
    make_valid,
    points,
    unary_union,
)
from shapely.geometry import (
    GeometryCollection,
    MultiPolygon,
    Polygon,
    box,
    mapping,
    shape,
)


LOGGER = logging.getLogger("select_crawdad_rois")

PRIMARY_ROIS = ("lateral", "dorsomedial", "ventromedial")
OUTSIDE_ROI = "outside_global_roi"
ROI_COLORS = {
    "lateral": "#E68613",
    "dorsomedial": "#4C78A8",
    "ventromedial": "#59A14F",
    OUTSIDE_ROI: "#666666",
}
D1_OUTSIDE_ROI_COLORS = {
    "A": "#E31A1C",
    "B": "#984EA3",
}
MSN_LABELS = (
    "DRD1_MSN",
    "DRD2_MSN",
    "D1_Island_A",
    "D1_Island_B",
    "MSN_Oligo",
)
SPLIT_QUANTILES = (0.0005, 0.9995)
DEFAULT_LATERAL_ROI = "H1_M3TCP9V_D1"
DEFAULT_NEIGHBORHOOD_DISTANCES = (50.0, 100.0)
DEFAULT_TILE_SIZES = (100.0, 200.0, 400.0)
REQUIRED_XENIUM_COLUMNS = {
    "donor",
    "sample",
    "cell_id",
    "cell_type",
    "x_second_pass",
    "y_second_pass",
    "z_height",
}
REQUIRED_VHD_COLUMNS = {
    "donor",
    "vhd_sample",
    "vhd_obs_id",
    "x_second_pass",
    "y_second_pass",
}


def git_root() -> Path:
    """Return the repository root."""
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def safe_name(value: str) -> str:
    """Convert a value to a filename- and column-safe component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    if not cleaned:
        raise ValueError(f"Cannot construct a safe name from {value!r}")
    return cleaned


def canonical_roi_name(value: str) -> str:
    """Convert a VHD sample label to the ROI naming convention."""
    value = str(value).strip().replace("-", "_")
    return value[4:] if value.startswith("VHD_") else value


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Construct three mutually exclusive CRAWDAD ROIs from aligned "
            "Xenium and Visium HD coordinates."
        )
    )
    parser.add_argument("--donor", default="Br6660")
    parser.add_argument(
        "--xenium-coordinate-csv",
        type=Path,
        help="Final donor-level Xenium second-pass coordinate CSV.",
    )
    parser.add_argument(
        "--vhd-coordinate-dir",
        type=Path,
        help=(
            "Module 02 VHD second-pass base directory. All pair-specific "
            "coordinate CSVs for --donor are discovered recursively."
        ),
    )
    parser.add_argument(
        "--vhd-coordinate-files",
        type=Path,
        nargs="+",
        help="Explicit final VHD second-pass coordinate CSVs.",
    )
    parser.add_argument(
        "--vhd-roi-file",
        type=Path,
        help=(
            "Optional precomputed six-ROI GeoJSON or vertex CSV. When given, "
            "VHD coordinate files are not used to construct hulls."
        ),
    )
    parser.add_argument("--roi-name-column", default="roi_name")
    parser.add_argument(
        "--lateral-roi-name", default=DEFAULT_LATERAL_ROI
    )
    parser.add_argument(
        "--hull-method",
        choices=("concave_rectangle", "concave", "convex", "rectangle"),
        default="concave_rectangle",
        help=(
            "VHD footprint construction method. 'concave_rectangle' builds "
            "a concave hull and then its minimum rotated rectangle."
        ),
    )
    parser.add_argument(
        "--hull-ratio",
        type=float,
        default=0.05,
        help="Shapely concave-hull ratio in [0,1].",
    )
    parser.add_argument(
        "--tissue-polygon",
        type=Path,
        help="Optional Xenium tissue polygon in GeoJSON or vertex CSV format.",
    )
    parser.add_argument(
        "--coordinate-unit",
        default="micrometer",
        help="Physical unit of the aligned coordinates; must be micrometers.",
    )
    parser.add_argument(
        "--geometry-tolerance-um2",
        type=float,
        default=1e-6,
        help="Maximum accepted positive-area ROI overlap or reconstruction error.",
    )
    parser.add_argument(
        "--d1a-label", default="D1_Island_A"
    )
    parser.add_argument(
        "--d1b-label", default="D1_Island_B"
    )
    parser.add_argument(
        "--neighborhood-distances-um",
        nargs="+",
        type=float,
        default=list(DEFAULT_NEIGHBORHOOD_DISTANCES),
    )
    parser.add_argument(
        "--tile-sizes-um",
        nargs="+",
        type=float,
        default=list(DEFAULT_TILE_SIZES),
    )
    parser.add_argument("--low-reference-cell-threshold", type=int, default=20)
    parser.add_argument("--small-roi-area-um2", type=float, default=100_000.0)
    parser.add_argument(
        "--high-boundary-fraction-threshold", type=float, default=0.5
    )
    parser.add_argument(
        "--large-outside-fraction-warning", type=float, default=0.5
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help="Processed-data output directory.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help="Plot output directory; all figures are written as PNG.",
    )
    parser.add_argument(
        "--max-plot-cells",
        type=int,
        default=500_000,
        help="Deterministic display cap per figure; 0 plots all cells.",
    )
    parser.add_argument(
        "--preview-only",
        action="store_true",
        help=(
            "Write resolved config, geometry QC, and preview figures, but no "
            "final cell assignments or CRAWDAD inputs."
        ),
    )
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--log-level", default="INFO")
    args = parser.parse_args(argv)

    if args.vhd_coordinate_files and args.vhd_roi_file:
        parser.error(
            "--vhd-coordinate-files and --vhd-roi-file are mutually exclusive"
        )
    if not 0 <= args.hull_ratio <= 1:
        parser.error("--hull-ratio must be between 0 and 1")
    for name in ("geometry_tolerance_um2", "small_roi_area_um2"):
        value = getattr(args, name)
        if value is not None and value < 0:
            parser.error(f"--{name.replace('_', '-')} must be non-negative")
    normalized_unit = str(args.coordinate_unit).strip().lower()
    if normalized_unit not in {
        "um",
        "µm",
        "micron",
        "microns",
        "micrometer",
        "micrometers",
    }:
        parser.error("--coordinate-unit must describe micrometers")
    args.coordinate_unit = "micrometer"
    if args.low_reference_cell_threshold < 0:
        parser.error("--low-reference-cell-threshold must be non-negative")
    if not 0 <= args.high_boundary_fraction_threshold <= 1:
        parser.error("--high-boundary-fraction-threshold must be in [0,1]")
    if not 0 <= args.large_outside_fraction_warning <= 1:
        parser.error("--large-outside-fraction-warning must be in [0,1]")
    if any(value <= 0 for value in args.neighborhood_distances_um):
        parser.error("--neighborhood-distances-um values must be positive")
    if any(value <= 0 for value in args.tile_sizes_um):
        parser.error("--tile-sizes-um values must be positive")
    if args.max_plot_cells < 0:
        parser.error("--max-plot-cells must be non-negative")
    return args


def resolve_paths(args: argparse.Namespace, root: Path) -> None:
    """Fill pipeline-aware input and output defaults."""
    args.xenium_coordinate_csv = args.xenium_coordinate_csv or (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / f"module02_Spateo_second_pass_{args.donor}"
        / f"{args.donor}_second_pass_coordinates.csv"
    )
    args.vhd_coordinate_dir = args.vhd_coordinate_dir or (
        root
        / "processed-data"
        / "11_Xenium_VisiumHD_alignment"
        / "module_02_Spateo_second_pass"
    )
    module_name = "module01_select_crawdad_rois"
    args.output_dir = args.output_dir or (
        root / "processed-data" / "10_Xenium_CRAWDAD" / module_name
    )
    args.plot_dir = args.plot_dir or (
        root / "plots" / "10_Xenium_CRAWDAD" / module_name
    )


def prepare_output_dirs(args: argparse.Namespace) -> dict[str, Path]:
    """Create module output directories and enforce overwrite protection."""
    output_dirs = {
        "root": args.output_dir,
        "polygons": args.output_dir / "roi_polygons",
        "assignments": args.output_dir / "cell_assignments",
        "crawdad": args.output_dir / "crawdad_inputs",
        "qc": args.output_dir / "qc",
        "plots": args.plot_dir,
    }
    protected = [
        output_dirs["root"] / "config_resolved.yaml",
        output_dirs["root"] / "metadata.json",
        output_dirs["polygons"] / "original_six_vhd_rois.geojson",
        output_dirs["polygons"] / "primary_three_rois.geojson",
        output_dirs["assignments"] / "all_cells_roi_assignment.csv.gz",
        output_dirs["assignments"] / "all_cells_roi_assignment.parquet",
        output_dirs["qc"] / "validation_report.txt",
    ]
    existing = [path for path in protected if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Output files already exist; rerun with --overwrite: "
            + ", ".join(map(str, existing))
        )
    for path in output_dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    if args.overwrite:
        for pdf_path in output_dirs["plots"].glob("*.pdf"):
            pdf_path.unlink()
    return output_dirs


def load_xenium_coordinates(
    path: Path, donor: str
) -> pd.DataFrame:
    """Load and validate final donor-level aligned Xenium coordinates."""
    if not path.is_file():
        raise FileNotFoundError(path)
    data = pd.read_csv(
        path,
        dtype={
            "donor": "string",
            "sample": "string",
            "cell_id": "string",
            "cell_type": "string",
        },
    )
    missing = REQUIRED_XENIUM_COLUMNS - set(data.columns)
    if missing:
        raise ValueError(
            f"Xenium coordinate CSV is missing columns: {sorted(missing)}"
        )
    data = data.loc[data["donor"].astype(str).eq(donor)].copy()
    if data.empty:
        raise ValueError(f"No {donor} cells found in {path}")
    if data["cell_id"].isna().any() or data["cell_id"].duplicated().any():
        raise ValueError("Xenium cell_id values must be present and unique")
    required_nonmissing = [
        "sample",
        "cell_type",
        "x_second_pass",
        "y_second_pass",
        "z_height",
    ]
    if data[required_nonmissing].isna().any().any():
        raise ValueError("Xenium coordinates or annotations contain missing values")
    numeric = ["x_second_pass", "y_second_pass", "z_height"]
    data[numeric] = data[numeric].apply(pd.to_numeric, errors="raise")
    if not np.isfinite(data[numeric].to_numpy()).all():
        raise ValueError("Xenium coordinates contain non-finite values")
    sample_depths = data.groupby("sample", observed=True)["z_height"].nunique()
    if (sample_depths != 1).any():
        raise ValueError("Every Xenium sample must have exactly one z_height")
    return data


def discover_vhd_coordinate_files(
    args: argparse.Namespace,
) -> list[Path]:
    """Discover the donor's final VHD coordinate tables."""
    if args.vhd_coordinate_files:
        files = list(args.vhd_coordinate_files)
    else:
        if not args.vhd_coordinate_dir.is_dir():
            raise FileNotFoundError(args.vhd_coordinate_dir)
        files = sorted(
            args.vhd_coordinate_dir.glob(
                f"{args.donor}*/*_second_pass_coordinates.csv"
            )
        )
    missing = [path for path in files if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Missing VHD coordinate files: " + ", ".join(map(str, missing))
        )
    if len(files) != 6:
        raise ValueError(
            f"Expected exactly six aligned VHD coordinate files for "
            f"{args.donor}; found {len(files)}"
        )
    return files


def repair_polygonal(geometry: Any) -> tuple[Any, bool]:
    """Repair a geometry and retain all polygonal components."""
    repaired = False
    if geometry.is_empty:
        raise ValueError("Encountered an empty ROI geometry")
    if not geometry.is_valid:
        geometry = make_valid(geometry)
        repaired = True
    if isinstance(geometry, Polygon | MultiPolygon):
        polygonal = geometry
    elif isinstance(geometry, GeometryCollection):
        parts = [
            item
            for item in geometry.geoms
            if isinstance(item, Polygon | MultiPolygon) and not item.is_empty
        ]
        if not parts:
            raise ValueError("Geometry repair produced no polygonal components")
        polygonal = unary_union(parts)
    else:
        raise ValueError(
            f"Expected Polygon/MultiPolygon; observed {geometry.geom_type}"
        )
    if not polygonal.is_valid:
        polygonal = make_valid(polygonal)
        repaired = True
    if polygonal.is_empty:
        raise ValueError("Polygon repair produced an empty geometry")
    return polygonal, repaired


def polygon_components(geometry: Any) -> list[Polygon]:
    """Return all polygon components without dropping small pieces."""
    if isinstance(geometry, Polygon):
        return [geometry]
    if isinstance(geometry, MultiPolygon):
        return list(geometry.geoms)
    return [
        item
        for item in getattr(geometry, "geoms", [])
        if isinstance(item, Polygon)
    ]


def load_polygon_file(
    path: Path, roi_name_column: str
) -> dict[str, Any]:
    """Load named polygons from GeoJSON or a vertex CSV."""
    if not path.is_file():
        raise FileNotFoundError(path)
    if path.suffix.lower() in {".json", ".geojson"}:
        payload = json.loads(path.read_text())
        features = (
            payload.get("features", [])
            if payload.get("type") == "FeatureCollection"
            else [payload]
        )
        result: dict[str, Any] = {}
        for feature in features:
            properties = feature.get("properties", {})
            name = properties.get(roi_name_column)
            if name is None:
                raise ValueError(
                    f"GeoJSON feature is missing property {roi_name_column!r}"
                )
            geometry, _ = repair_polygonal(shape(feature["geometry"]))
            result[str(name)] = geometry
        return result
    if path.suffix.lower() == ".csv":
        vertices = pd.read_csv(path)
        required = {roi_name_column, "x", "y"}
        missing = required - set(vertices.columns)
        if missing:
            raise ValueError(
                f"Polygon CSV is missing columns: {sorted(missing)}"
            )
        order_column = (
            "vertex_order" if "vertex_order" in vertices.columns else None
        )
        result = {}
        for name, frame in vertices.groupby(roi_name_column, sort=False):
            if order_column:
                frame = frame.sort_values(order_column, kind="stable")
            geometry, _ = repair_polygonal(
                Polygon(frame[["x", "y"]].to_numpy(dtype=float))
            )
            result[str(name)] = geometry
        return result
    raise ValueError("Polygon input must be GeoJSON/JSON or CSV")


def construct_vhd_footprints(
    files: Sequence[Path],
    donor: str,
    hull_method: str,
    hull_ratio: float,
) -> tuple[dict[str, Any], dict[str, Any], list[dict[str, Any]]]:
    """Construct named VHD footprints in the aligned Xenium frame."""
    footprints: dict[str, Any] = {}
    concave_hulls: dict[str, Any] = {}
    records: list[dict[str, Any]] = []
    for path in files:
        data = pd.read_csv(
            path,
            usecols=list(REQUIRED_VHD_COLUMNS),
            dtype={
                "donor": "string",
                "vhd_sample": "string",
                "vhd_obs_id": "string",
            },
        )
        observed_donors = set(data["donor"].dropna().astype(str))
        if observed_donors != {donor}:
            raise ValueError(
                f"{path}: expected donor {donor}; observed "
                f"{sorted(observed_donors)}"
            )
        sample_values = data["vhd_sample"].dropna().astype(str).unique()
        if len(sample_values) != 1:
            raise ValueError(f"{path}: expected exactly one vhd_sample")
        if data["vhd_obs_id"].isna().any() or data["vhd_obs_id"].duplicated().any():
            raise ValueError(f"{path}: VHD observation IDs must be unique")
        coordinates = data[
            ["x_second_pass", "y_second_pass"]
        ].to_numpy(dtype=float)
        if not np.isfinite(coordinates).all():
            raise ValueError(f"{path}: non-finite aligned VHD coordinates")
        roi_name = canonical_roi_name(sample_values[0])
        if roi_name in footprints:
            raise ValueError(f"Duplicate VHD ROI name: {roi_name}")
        concave_repaired = False
        if hull_method == "rectangle":
            min_x, min_y = coordinates.min(axis=0)
            max_x, max_y = coordinates.max(axis=0)
            geometry = box(min_x, min_y, max_x, max_y)
        elif hull_method in {"concave", "concave_rectangle"}:
            cloud = MultiPoint(coordinates)
            concave_geometry = concave_hull(
                cloud, ratio=hull_ratio, allow_holes=False
            )
            concave_geometry, concave_repaired = repair_polygonal(
                concave_geometry
            )
            concave_hulls[roi_name] = concave_geometry
            geometry = concave_geometry
            if hull_method == "concave_rectangle":
                geometry = geometry.minimum_rotated_rectangle
        else:
            cloud = MultiPoint(coordinates)
            geometry = cloud.convex_hull
        geometry, repaired = repair_polygonal(geometry)
        footprints[roi_name] = geometry
        records.append(
            {
                "roi_name": roi_name,
                "source_file": str(path.resolve()),
                "n_vhd_observations": len(data),
                "hull_method": hull_method,
                "hull_ratio": (
                    hull_ratio
                    if hull_method in {"concave", "concave_rectangle"}
                    else np.nan
                ),
                "concave_hull_geometry_repaired": concave_repaired,
                "geometry_repaired": repaired,
            }
        )
    return footprints, concave_hulls, records


def construct_primary_rois(
    footprints: dict[str, Any],
    lateral_name: str,
    tissue_polygon: Any | None,
    xenium: pd.DataFrame,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Build lateral, dorsomedial, and ventromedial polygons."""
    if lateral_name not in footprints:
        raise ValueError(
            f"Lateral ROI {lateral_name!r} is absent; available="
            f"{sorted(footprints)}"
        )
    if len(footprints) != 6:
        raise ValueError(
            f"Exactly six original VHD ROIs are required; found {len(footprints)}"
        )
    lateral_raw = footprints[lateral_name]
    medial_names = [name for name in footprints if name != lateral_name]
    medial_union_raw = unary_union([footprints[name] for name in medial_names])
    overlap_assigned_lateral = medial_union_raw.intersection(lateral_raw)
    medial_raw = medial_union_raw.difference(lateral_raw)

    lateral = lateral_raw
    medial = medial_raw
    if tissue_polygon is not None:
        lateral = lateral.intersection(tissue_polygon)
        medial = medial.intersection(tissue_polygon)
    lateral, lateral_repaired = repair_polygonal(lateral)
    medial, medial_repaired = repair_polygonal(medial)

    min_x, medial_y_min, max_x, medial_y_max = medial.bounds
    x = xenium["x_second_pass"].to_numpy(dtype=float)
    y = xenium["y_second_pass"].to_numpy(dtype=float)
    finite_xy = np.isfinite(x) & np.isfinite(y)
    in_medial = intersects_xy(medial, x, y)
    is_broad_msn = (
        xenium["cell_type"].astype(str).isin(MSN_LABELS).to_numpy()
    )
    split_reference = finite_xy & in_medial & is_broad_msn
    split_reference_indices = np.flatnonzero(split_reference)
    split_reference_cell_count = int(len(split_reference_indices))
    if split_reference_cell_count == 0:
        raise ValueError(
            "No broad MSN Xenium cells were found inside medial_roi; "
            "cannot calculate the dorsal-ventral dividing line."
        )
    y_reference = y[split_reference]
    y_low, y_high = np.quantile(y_reference, SPLIT_QUANTILES)
    y_split = (y_low + y_high) / 2.0
    if (
        not np.isfinite(y_split)
        or not y_low < y_high
        or not y_low < y_split < y_high
    ):
        raise ValueError(
            "Invalid Xenium-derived dorsal-ventral range: "
            f"y_low={y_low}, y_high={y_high}, y_split={y_split}"
        )

    margin = max(
        max_x - min_x,
        medial_y_max - medial_y_min,
    ) + 1.0
    dorsal_half = box(
        min_x - margin,
        y_split,
        max_x + margin,
        medial_y_max + margin,
    )
    ventral_half = box(
        min_x - margin,
        medial_y_min - margin,
        max_x + margin,
        y_split,
    )
    dorsomedial_raw = medial.intersection(dorsal_half)
    ventromedial_raw = medial.intersection(ventral_half)
    if dorsomedial_raw.is_empty or ventromedial_raw.is_empty:
        raise ValueError(
            "The Xenium-derived dividing line produced an empty medial "
            "subdivision."
        )
    dorsomedial, dm_repaired = repair_polygonal(dorsomedial_raw)
    ventromedial, vm_repaired = repair_polygonal(ventromedial_raw)

    primary = {
        "lateral": lateral,
        "dorsomedial": dorsomedial,
        "ventromedial": ventromedial,
    }
    assigned_medial = unary_union([dorsomedial, ventromedial])
    missing_medial = medial.difference(assigned_medial)
    extra_medial = assigned_medial.difference(medial)
    context = {
        "medial_roi_names": medial_names,
        "medial_union_raw": medial_union_raw,
        "lateral_raw": lateral_raw,
        "overlap_assigned_lateral": overlap_assigned_lateral,
        "medial_after_lateral_removal": medial,
        "medial_y_min": float(medial_y_min),
        "medial_y_max": float(medial_y_max),
        "split_quantile_low": SPLIT_QUANTILES[0],
        "split_quantile_high": SPLIT_QUANTILES[1],
        "split_reference_celltypes": list(MSN_LABELS),
        "split_reference_cell_count": split_reference_cell_count,
        "split_reference_indices": split_reference_indices,
        "y_low": float(y_low),
        "y_high": float(y_high),
        "y_split": float(y_split),
        "dorsal_direction": "greater_y",
        "missing_medial_after_split": missing_medial,
        "extra_medial_after_split": extra_medial,
        "repair_flags": {
            "lateral": lateral_repaired,
            "medial": medial_repaired,
            "dorsomedial": dm_repaired,
            "ventromedial": vm_repaired,
        },
    }
    return primary, context


def build_boolean_operation_qc(
    primary: dict[str, Any], context: dict[str, Any]
) -> pd.DataFrame:
    """Report areas used to construct and validate the primary ROIs."""
    return pd.DataFrame(
        [
            {
                "measurement": "five_medial_vhd_union",
                "area_um2": context["medial_union_raw"].area,
            },
            {
                "measurement": "medial_lateral_overlap_assigned_to_lateral",
                "area_um2": context["overlap_assigned_lateral"].area,
            },
            {
                "measurement": "medial_after_lateral_removal_and_tissue_clip",
                "area_um2": context["medial_after_lateral_removal"].area,
            },
            {
                "measurement": "dorsomedial_assigned",
                "area_um2": primary["dorsomedial"].area,
            },
            {
                "measurement": "ventromedial_assigned",
                "area_um2": primary["ventromedial"].area,
            },
            {
                "measurement": "medial_unassigned_after_split",
                "area_um2": context["missing_medial_after_split"].area,
            },
            {
                "measurement": "extra_medial_after_split",
                "area_um2": context["extra_medial_after_split"].area,
            },
        ]
    )


def write_geojson(
    geometries: dict[str, Any],
    path: Path,
    extra_properties: dict[str, dict[str, Any]] | None = None,
) -> None:
    """Write named polygon geometries as a GeoJSON FeatureCollection."""
    features = []
    for name, geometry in geometries.items():
        properties = {"roi_name": name}
        if extra_properties and name in extra_properties:
            properties.update(extra_properties[name])
        features.append(
            {
                "type": "Feature",
                "properties": properties,
                "geometry": mapping(geometry),
            }
        )
    path.write_text(
        json.dumps(
            {"type": "FeatureCollection", "features": features}, indent=2
        )
        + "\n"
    )


def geometry_record(
    name: str,
    stage: str,
    geometry: Any,
    *,
    area_before_clip: float,
    area_after_clip: float,
    repaired: bool,
    fragment_removed_area: float,
) -> dict[str, Any]:
    """Create one geometry-QC record."""
    min_x, min_y, max_x, max_y = geometry.bounds
    return {
        "stage": stage,
        "roi_name": name,
        "geometry_valid": bool(geometry.is_valid),
        "geometry_type": geometry.geom_type,
        "total_area_um2": float(geometry.area),
        "area_before_tissue_clip_um2": float(area_before_clip),
        "area_after_tissue_clip_um2": float(area_after_clip),
        "n_components": len(polygon_components(geometry)),
        "min_x": min_x,
        "min_y": min_y,
        "max_x": max_x,
        "max_y": max_y,
        "width_um": max_x - min_x,
        "height_um": max_y - min_y,
        "geometry_repaired": repaired,
        "small_fragment_area_removed_um2": fragment_removed_area,
    }


def build_geometry_qc(
    footprints: dict[str, Any],
    primary: dict[str, Any],
    context: dict[str, Any],
    tissue_polygon: Any | None,
) -> pd.DataFrame:
    """Summarize original and derived geometry properties."""
    records = []
    for name, geometry in footprints.items():
        clipped = (
            geometry.intersection(tissue_polygon)
            if tissue_polygon is not None
            else geometry
        )
        records.append(
            geometry_record(
                name,
                "original_vhd",
                geometry,
                area_before_clip=geometry.area,
                area_after_clip=clipped.area,
                repaired=False,
                fragment_removed_area=0.0,
            )
        )
    for name, geometry in primary.items():
        raw = (
            context["lateral_raw"]
            if name == "lateral"
            else context["medial_after_lateral_removal"]
        )
        records.append(
            geometry_record(
                name,
                "primary",
                geometry,
                area_before_clip=raw.area,
                area_after_clip=geometry.area,
                repaired=context["repair_flags"][name],
                fragment_removed_area=0.0,
            )
        )
    result = pd.DataFrame.from_records(records)
    result["medial_y_min"] = context["medial_y_min"]
    result["medial_y_max"] = context["medial_y_max"]
    result["split_reference_cell_count"] = context[
        "split_reference_cell_count"
    ]
    result["y_low"] = context["y_low"]
    result["y_high"] = context["y_high"]
    result["y_split"] = context["y_split"]
    result["dorsal_direction"] = context["dorsal_direction"]
    return result


def validate_primary_geometry(
    primary: dict[str, Any],
    context: dict[str, Any],
    tolerance_um2: float,
) -> pd.DataFrame:
    """Validate pairwise overlap and reconstruction of the medial ROI."""
    rows = []
    for index, first in enumerate(PRIMARY_ROIS):
        for second in PRIMARY_ROIS[index + 1 :]:
            area = float(primary[first].intersection(primary[second]).area)
            rows.append(
                {
                    "roi_1": first,
                    "roi_2": second,
                    "intersection_area_um2": area,
                    "tolerance_area_um2": tolerance_um2,
                    "passes": area <= tolerance_um2,
                }
            )
    result = pd.DataFrame(rows)
    if not result["passes"].all():
        raise ValueError(
            "Final primary ROI overlap exceeds tolerance:\n"
            + result.to_string(index=False)
        )
    missing_area = float(context["missing_medial_after_split"].area)
    extra_area = float(context["extra_medial_after_split"].area)
    if missing_area > tolerance_um2 or extra_area > tolerance_um2:
        raise ValueError(
            "Dorsomedial and ventromedial polygons do not reconstruct the "
            f"medial ROI within tolerance: missing={missing_area:.6f} µm², "
            f"extra={extra_area:.6f} µm²"
        )
    return result


def plot_indices(n_rows: int, maximum: int) -> np.ndarray:
    """Return deterministic, evenly spaced plotting indices."""
    if maximum == 0 or n_rows <= maximum:
        return np.arange(n_rows)
    return np.linspace(0, n_rows - 1, maximum, dtype=np.int64)


def draw_geometry(
    axis: plt.Axes,
    geometry: Any,
    *,
    edgecolor: str,
    facecolor: str = "none",
    linewidth: float = 1.5,
    alpha: float = 1.0,
    zorder: int = 3,
) -> None:
    """Draw every Polygon component, including holes."""
    for polygon in polygon_components(geometry):
        exterior = np.asarray(polygon.exterior.coords)
        axis.add_patch(
            MplPolygon(
                exterior,
                closed=True,
                facecolor=facecolor,
                edgecolor=edgecolor,
                linewidth=linewidth,
                alpha=alpha,
                zorder=zorder,
            )
        )
        for interior in polygon.interiors:
            hole = np.asarray(interior.coords)
            axis.add_patch(
                MplPolygon(
                    hole,
                    closed=True,
                    facecolor="white",
                    edgecolor=edgecolor,
                    linewidth=max(0.5, linewidth * 0.7),
                    alpha=1.0,
                    zorder=zorder + 0.1,
                )
            )


def save_figure(
    figure: plt.Figure, plot_dir: Path, stem: str, dpi: int = 300
) -> None:
    """Save a server-safe PNG figure."""
    figure.savefig(
        plot_dir / f"{stem}.png",
        dpi=dpi,
        bbox_inches="tight",
    )
    plt.close(figure)


def set_geometry_limits(axis: plt.Axes, geometries: Iterable[Any]) -> None:
    """Limit an axis to the supplied geometries with a small visual margin."""
    min_x, min_y, max_x, max_y = unary_union(list(geometries)).bounds
    x_margin = 0.05 * max(max_x - min_x, 1.0)
    y_margin = 0.05 * max(max_y - min_y, 1.0)
    axis.set_xlim(min_x - x_margin, max_x + x_margin)
    axis.set_ylim(min_y - y_margin, max_y + y_margin)


def plot_original_footprints(
    cells: pd.DataFrame,
    footprints: dict[str, Any],
    lateral_name: str,
    plot_dir: Path,
    maximum_cells: int,
    *,
    title: str = "Original aligned VHD footprints over Xenium",
    stem: str = "figure01_original_vhd_footprints",
) -> None:
    """Figure 1: aligned Xenium cells and all original VHD footprints."""
    indices = plot_indices(len(cells), maximum_cells)
    figure, axis = plt.subplots(figsize=(9, 9))
    axis.scatter(
        cells["x_aligned"].to_numpy()[indices],
        cells["y_aligned"].to_numpy()[indices],
        s=0.08,
        c="#D0D0D0",
        alpha=0.35,
        linewidths=0,
        rasterized=True,
    )
    colors = plt.get_cmap("tab10")(
        np.linspace(0, 1, len(footprints))
    )
    handles = []
    for color, (name, geometry) in zip(colors, sorted(footprints.items())):
        linewidth = 3.0 if name == lateral_name else 1.7
        draw_geometry(
            axis,
            geometry,
            edgecolor=color,
            linewidth=linewidth,
        )
        label = f"{name} (lateral)" if name == lateral_name else name
        handles.append(Line2D([0], [0], color=color, lw=linewidth, label=label))
    axis.legend(handles=handles, loc="center left", bbox_to_anchor=(1.02, 0.5))
    axis.set_title(title)
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_aspect("equal")
    save_figure(figure, plot_dir, stem)


def plot_split_preview(
    cells: pd.DataFrame,
    footprints: dict[str, Any],
    context: dict[str, Any],
    lateral_name: str,
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Preview: proposed dorsal-ventral split before assignments."""
    indices = plot_indices(len(cells), maximum_cells)
    figure, axis = plt.subplots(figsize=(9, 9))
    axis.scatter(
        cells["x_aligned"].to_numpy()[indices],
        cells["y_aligned"].to_numpy()[indices],
        s=0.08,
        c="#D0D0D0",
        alpha=0.3,
        linewidths=0,
        rasterized=True,
    )
    for name, geometry in sorted(footprints.items()):
        draw_geometry(
            axis,
            geometry,
            edgecolor="#E68613" if name == lateral_name else "#4C78A8",
            linewidth=2.0 if name == lateral_name else 1.0,
        )
    axis.axhline(
        context["y_split"],
        color="black",
        linestyle="--",
        linewidth=2,
        label=f"Broad-MSN Xenium y_split = {context['y_split']:.2f} µm",
    )
    axis.legend(loc="upper right")
    axis.set_title("Dorsal–ventral split preview")
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_aspect("equal")
    set_geometry_limits(axis, footprints.values())
    save_figure(figure, plot_dir, "roi_split_preview")


def plot_roi_construction(
    cells: pd.DataFrame,
    footprints: dict[str, Any],
    lateral_name: str,
    context: dict[str, Any],
    primary: dict[str, Any],
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Figure 2: auditable Boolean ROI construction."""
    figure, axis = plt.subplots(figsize=(9, 9))
    draw_geometry(
        axis,
        context["medial_union_raw"],
        edgecolor="#4C78A8",
        facecolor="#4C78A8",
        alpha=0.18,
    )
    draw_geometry(
        axis,
        context["lateral_raw"],
        edgecolor=ROI_COLORS["lateral"],
        facecolor=ROI_COLORS["lateral"],
        alpha=0.22,
        linewidth=2.2,
    )
    reference = cells.iloc[context["split_reference_indices"]]
    reference = reference.iloc[
        plot_indices(len(reference), maximum_cells)
    ]
    axis.scatter(
        reference["x_aligned"],
        reference["y_aligned"],
        s=0.2,
        c="#7A0177",
        alpha=0.45,
        linewidths=0,
        rasterized=True,
        zorder=3.5,
    )
    if not context["overlap_assigned_lateral"].is_empty:
        draw_geometry(
            axis,
            context["overlap_assigned_lateral"],
            edgecolor="#D62728",
            facecolor="#D62728",
            alpha=0.45,
            linewidth=1.5,
            zorder=4,
        )
    for name, geometry in sorted(footprints.items()):
        draw_geometry(
            axis,
            geometry,
            edgecolor=(
                ROI_COLORS["lateral"] if name == lateral_name else "#6F6F6F"
            ),
            linewidth=1.2 if name == lateral_name else 0.8,
            zorder=4,
        )
    draw_geometry(
        axis,
        primary["lateral"],
        edgecolor=ROI_COLORS["lateral"],
        linewidth=2.5,
        zorder=5,
    )
    draw_geometry(
        axis,
        primary["dorsomedial"],
        edgecolor=ROI_COLORS["dorsomedial"],
        linewidth=2.5,
        zorder=5,
    )
    draw_geometry(
        axis,
        primary["ventromedial"],
        edgecolor=ROI_COLORS["ventromedial"],
        linewidth=2.5,
        zorder=5,
    )
    bounds = unary_union(list(primary.values())).bounds
    y_split = context["y_split"]
    axis.axhline(
        context["y_low"],
        color="#7A0177",
        linestyle=":",
        linewidth=1.8,
    )
    axis.axhline(
        y_split,
        color="black",
        linestyle="--",
        linewidth=2,
    )
    axis.axhline(
        context["y_high"],
        color="#7A0177",
        linestyle=":",
        linewidth=1.8,
    )
    medial_bounds = context["medial_after_lateral_removal"].bounds
    label_x = medial_bounds[0] + 0.02 * max(
        medial_bounds[2] - medial_bounds[0],
        1.0,
    )
    label_offset = 0.04 * max(bounds[3] - bounds[1], 1.0)
    axis.text(
        label_x,
        y_split + label_offset,
        "dorsomedial (larger y)",
        ha="left",
        va="bottom",
        fontsize=10,
        weight="bold",
    )
    axis.text(
        label_x,
        y_split - label_offset,
        "ventromedial (smaller y)",
        ha="left",
        va="top",
        fontsize=10,
        weight="bold",
    )
    handles = [
        Line2D(
            [0],
            [0],
            color="#6F6F6F",
            lw=1,
            label="Six original VHD footprints",
        ),
        Patch(facecolor="#4C78A8", alpha=0.25, label="Five-footprint medial union"),
        Patch(facecolor=ROI_COLORS["lateral"], alpha=0.3, label="Lateral footprint"),
        Patch(
            facecolor="#D62728",
            alpha=0.5,
            label="Medial–lateral overlap assigned to lateral",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color="#7A0177",
            label="Broad MSN Xenium split reference",
        ),
        Line2D([0], [0], color=ROI_COLORS["dorsomedial"], lw=2.5, label="Dorsomedial"),
        Line2D([0], [0], color=ROI_COLORS["ventromedial"], lw=2.5, label="Ventromedial"),
        Line2D(
            [0],
            [0],
            color="#7A0177",
            ls=":",
            lw=1.8,
            label="Broad MSN y 0.05th/99.95th percentiles",
        ),
        Line2D([0], [0], color="black", ls="--", lw=2, label="Dorsal–ventral split"),
    ]
    axis.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        fontsize=8,
    )
    axis.set_title("ROI construction and Xenium-derived broad-MSN split")
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_aspect("equal")
    set_geometry_limits(axis, primary.values())
    save_figure(figure, plot_dir, "figure02_roi_construction")


def assign_cells(
    xenium: pd.DataFrame,
    footprints: dict[str, Any],
    primary: dict[str, Any],
    context: dict[str, Any],
) -> pd.DataFrame:
    """Assign cells to original footprints and one primary ROI."""
    result = xenium.rename(
        columns={
            "sample": "sample",
            "cell_type": "cell_type",
            "x_second_pass": "x_aligned",
            "y_second_pass": "y_aligned",
        }
    ).copy()
    x = result["x_aligned"].to_numpy(dtype=float)
    y = result["y_aligned"].to_numpy(dtype=float)

    for name, geometry in footprints.items():
        result[f"in_{safe_name(name)}"] = intersects_xy(geometry, x, y)

    in_lateral = intersects_xy(primary["lateral"], x, y)
    medial = context["medial_after_lateral_removal"]
    in_medial = intersects_xy(medial, x, y) & ~in_lateral
    y_split = float(context["y_split"])
    on_boundary = y == y_split
    in_dm = in_medial & (y >= y_split)
    in_vm = in_medial & (y < y_split)

    membership_count = (
        in_lateral.astype(np.int8)
        + in_dm.astype(np.int8)
        + in_vm.astype(np.int8)
    )
    if membership_count.max(initial=0) > 1:
        raise AssertionError("A Xenium cell was assigned to multiple primary ROIs")
    labels = np.full(len(result), OUTSIDE_ROI, dtype=object)
    labels[in_lateral] = "lateral"
    labels[in_dm] = "dorsomedial"
    labels[in_vm] = "ventromedial"
    result["primary_roi"] = labels
    result["on_dorsal_ventral_boundary"] = on_boundary & in_medial

    boundary_distance = np.full(len(result), np.nan, dtype=float)
    for roi_name in PRIMARY_ROIS:
        mask = labels == roi_name
        if mask.any():
            boundary_distance[mask] = distance(
                points(x[mask], y[mask]), primary[roi_name].boundary
            )
    result["distance_to_roi_boundary_um"] = boundary_distance

    keep_first = [
        "cell_id",
        "sample",
        "cell_type",
        "x_aligned",
        "y_aligned",
        "z_height",
        "primary_roi",
        "distance_to_roi_boundary_um",
        "on_dorsal_ventral_boundary",
    ]
    membership_columns = sorted(
        column for column in result if column.startswith("in_")
    )
    remaining = [
        column
        for column in result
        if column not in set(keep_first + membership_columns)
    ]
    return result[keep_first + membership_columns + remaining]


def plot_final_rois(
    assignments: pd.DataFrame,
    primary: dict[str, Any],
    context: dict[str, Any],
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Figure 3: all Xenium cells colored by primary ROI."""
    indices = plot_indices(len(assignments), maximum_cells)
    subset = assignments.iloc[indices]
    figure, axis = plt.subplots(figsize=(9, 9))
    for roi_name in (OUTSIDE_ROI, *PRIMARY_ROIS):
        frame = subset.loc[subset["primary_roi"].eq(roi_name)]
        axis.scatter(
            frame["x_aligned"],
            frame["y_aligned"],
            s=0.12 if roi_name == OUTSIDE_ROI else 0.18,
            c=ROI_COLORS[roi_name],
            alpha=0.35 if roi_name == OUTSIDE_ROI else 0.55,
            linewidths=0,
            label=roi_name,
            rasterized=True,
        )
    for roi_name, geometry in primary.items():
        draw_geometry(
            axis,
            geometry,
            edgecolor=ROI_COLORS[roi_name],
            linewidth=2.2,
        )
    axis.axhline(
        context["y_split"],
        color="black",
        linestyle="--",
        linewidth=1.5,
    )
    axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=5)
    axis.set_title("Final mutually exclusive Xenium ROIs")
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_aspect("equal")
    save_figure(figure, plot_dir, "figure03_final_primary_rois")


def plot_rois_by_slice(
    assignments: pd.DataFrame,
    primary: dict[str, Any],
    context: dict[str, Any],
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Figure 4: final ROI assignments in one panel per Xenium slice."""
    sample_order = (
        assignments.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    ncols = 4
    nrows = int(np.ceil(len(sample_order) / ncols))
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    per_panel_max = (
        0
        if maximum_cells == 0
        else max(1, maximum_cells // max(1, len(sample_order)))
    )
    for axis, sample in zip(axes.ravel(), sample_order):
        frame = assignments.loc[assignments["sample"].astype(str).eq(sample)]
        indices = plot_indices(len(frame), per_panel_max)
        frame = frame.iloc[indices]
        for roi_name in (OUTSIDE_ROI, *PRIMARY_ROIS):
            selected = frame.loc[frame["primary_roi"].eq(roi_name)]
            axis.scatter(
                selected["x_aligned"],
                selected["y_aligned"],
                s=0.08,
                c=ROI_COLORS[roi_name],
                alpha=0.4 if roi_name == OUTSIDE_ROI else 0.6,
                linewidths=0,
                rasterized=True,
            )
        for roi_name, geometry in primary.items():
            draw_geometry(
                axis,
                geometry,
                edgecolor=ROI_COLORS[roi_name],
                linewidth=0.8,
            )
        axis.axhline(
            context["y_split"],
            color="black",
            linestyle="--",
            linewidth=0.7,
        )
        depth = assignments.loc[
            assignments["sample"].astype(str).eq(sample), "z_height"
        ].iloc[0]
        axis.set_title(f"{sample}\n{depth:g} µm", fontsize=9)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    for axis in axes.ravel()[len(sample_order) :]:
        axis.axis("off")
    handles = [
        Patch(color=ROI_COLORS[name], label=name)
        for name in (*PRIMARY_ROIS, OUTSIDE_ROI)
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, -0.01),
    )
    figure.suptitle("Final CRAWDAD ROIs by Xenium slice", y=1.01)
    figure.tight_layout()
    save_figure(figure, plot_dir, "figure04_final_rois_by_slice")


def plot_d1_overlay(
    assignments: pd.DataFrame,
    primary: dict[str, Any],
    d1a_label: str,
    d1b_label: str,
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Figure 5: fixed-ROI D1 Island A/B descriptive QC."""
    inside = assignments["primary_roi"].isin(PRIMARY_ROIS)
    d1 = assignments.loc[
        inside & assignments["cell_type"].isin([d1a_label, d1b_label])
    ]
    indices = plot_indices(len(d1), maximum_cells)
    d1 = d1.iloc[indices]
    figure, axis = plt.subplots(figsize=(9, 9))
    colors = {d1a_label: "#224B82", d1b_label: "#00C853"}
    for label in (d1a_label, d1b_label):
        frame = d1.loc[d1["cell_type"].eq(label)]
        axis.scatter(
            frame["x_aligned"],
            frame["y_aligned"],
            s=0.3,
            color=colors[label],
            alpha=0.7,
            linewidths=0,
            label=label,
            rasterized=True,
        )
    for roi_name, geometry in primary.items():
        draw_geometry(
            axis,
            geometry,
            edgecolor=ROI_COLORS[roi_name],
            linewidth=2.0,
        )
    axis.legend(loc="center left", bbox_to_anchor=(1.02, 0.5), markerscale=4)
    axis.set_title("D1 Island A/B distribution within fixed ROIs (QC only)")
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_aspect("equal")
    save_figure(figure, plot_dir, "figure05_d1_island_roi_qc")


def plot_d1_overlay_by_slice(
    assignments: pd.DataFrame,
    primary: dict[str, Any],
    context: dict[str, Any],
    d1a_label: str,
    d1b_label: str,
    plot_dir: Path,
    maximum_cells: int,
) -> None:
    """Figure 5b: D1 Island cells and outside-ROI D1 cells by slice."""
    sample_order = (
        assignments.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    d1 = assignments.loc[
        assignments["cell_type"].isin([d1a_label, d1b_label])
    ]
    ncols = 4
    nrows = int(np.ceil(len(sample_order) / ncols))
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(4.2 * ncols, 4.0 * nrows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    per_panel_max = (
        0
        if maximum_cells == 0
        else max(1, maximum_cells // max(1, len(sample_order)))
    )
    colors = {
        d1a_label: "#224B82",
        d1b_label: "#00A65A",
    }
    for axis, sample in zip(axes.ravel(), sample_order):
        frame = d1.loc[d1["sample"].astype(str).eq(sample)]
        indices = plot_indices(len(frame), per_panel_max)
        frame = frame.iloc[indices]
        outside = frame.loc[frame["primary_roi"].eq(OUTSIDE_ROI)]
        for label, color_key in ((d1a_label, "A"), (d1b_label, "B")):
            selected = outside.loc[outside["cell_type"].eq(label)]
            axis.scatter(
                selected["x_aligned"],
                selected["y_aligned"],
                s=0.35,
                c=D1_OUTSIDE_ROI_COLORS[color_key],
                alpha=0.8,
                linewidths=0,
                rasterized=True,
            )
        inside = frame.loc[frame["primary_roi"].isin(PRIMARY_ROIS)]
        for label in (d1a_label, d1b_label):
            selected = inside.loc[inside["cell_type"].eq(label)]
            axis.scatter(
                selected["x_aligned"],
                selected["y_aligned"],
                s=0.25,
                c=colors[label],
                alpha=0.7,
                linewidths=0,
                rasterized=True,
            )
        for roi_name, geometry in primary.items():
            draw_geometry(
                axis,
                geometry,
                edgecolor=ROI_COLORS[roi_name],
                linewidth=0.8,
            )
        axis.axhline(
            context["y_split"],
            color="black",
            linestyle="--",
            linewidth=0.7,
        )
        depth = assignments.loc[
            assignments["sample"].astype(str).eq(sample), "z_height"
        ].iloc[0]
        axis.set_title(f"{sample}\n{depth:g} µm", fontsize=9)
        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
    for axis in axes.ravel()[len(sample_order) :]:
        axis.axis("off")
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=colors[d1a_label],
            label=f"{d1a_label} (inside global ROI)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=colors[d1b_label],
            label=f"{d1b_label} (inside global ROI)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=D1_OUTSIDE_ROI_COLORS["A"],
            label=f"{d1a_label} (outside global ROI)",
        ),
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=D1_OUTSIDE_ROI_COLORS["B"],
            label=f"{d1b_label} (outside global ROI)",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        ncol=4,
        frameon=False,
        bbox_to_anchor=(0.5, -0.005),
    )
    figure.suptitle("D1 Island A/B and outside-global-ROI cells by slice", y=1.01)
    figure.tight_layout()
    save_figure(figure, plot_dir, "figure05b_d1_island_roi_qc_by_slice")


def complete_sample_roi_index(
    assignments: pd.DataFrame,
) -> pd.MultiIndex:
    """Return every sample × primary-ROI combination."""
    samples = (
        assignments.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
        .index.astype(str)
    )
    return pd.MultiIndex.from_product(
        [samples, PRIMARY_ROIS], names=["sample", "primary_roi"]
    )


def build_qc_tables(
    assignments: pd.DataFrame,
    primary: dict[str, Any],
    d1a_label: str,
    d1b_label: str,
    neighborhood_distances: Sequence[float],
    tile_sizes: Sequence[float],
    low_reference_threshold: int,
    small_roi_area: float,
    high_boundary_fraction_threshold: float,
) -> dict[str, pd.DataFrame]:
    """Build cell-count, boundary, area, tile, and eligibility QC tables."""
    inside = assignments.loc[
        assignments["primary_roi"].isin(PRIMARY_ROIS)
    ].copy()
    full_index = complete_sample_roi_index(assignments)
    sample_order = (
        full_index.get_level_values("sample").drop_duplicates().tolist()
    )

    total_counts = (
        inside.groupby(["sample", "primary_roi"], observed=True)
        .size()
        .reindex(full_index, fill_value=0)
        .rename("n_cells")
        .reset_index()
    )
    sample_totals = assignments.groupby("sample", observed=True).size()
    total_counts["fraction_of_sample_cells"] = [
        count / sample_totals.loc[sample]
        for sample, count in zip(
            total_counts["sample"], total_counts["n_cells"]
        )
    ]

    celltype_counts = (
        inside.groupby(
            ["sample", "primary_roi", "cell_type"], observed=True
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    roi_totals = total_counts.set_index(
        ["sample", "primary_roi"]
    )["n_cells"]
    celltype_counts["fraction_within_roi"] = [
        count / roi_totals.loc[(sample, roi)] if roi_totals.loc[(sample, roi)] else np.nan
        for sample, roi, count in zip(
            celltype_counts["sample"],
            celltype_counts["primary_roi"],
            celltype_counts["n_cells"],
        )
    ]

    d1 = inside.loc[inside["cell_type"].isin([d1a_label, d1b_label])]
    d1_counts = (
        d1.groupby(
            ["sample", "primary_roi", "cell_type"], observed=True
        )
        .size()
        .unstack("cell_type", fill_value=0)
        .reindex(full_index, fill_value=0)
        .reset_index()
    )
    for label in (d1a_label, d1b_label):
        if label not in d1_counts:
            d1_counts[label] = 0
    d1_counts = d1_counts[
        ["sample", "primary_roi", d1a_label, d1b_label]
    ].rename(
        columns={d1a_label: "n_d1a", d1b_label: "n_d1b"}
    )
    d1_counts["n_d1_total"] = d1_counts["n_d1a"] + d1_counts["n_d1b"]

    area_rows = []
    for sample in sample_order:
        for roi_name in PRIMARY_ROIS:
            geometry = primary[roi_name]
            min_x, min_y, max_x, max_y = geometry.bounds
            area_rows.append(
                {
                    "sample": str(sample),
                    "primary_roi": roi_name,
                    "area_um2": geometry.area,
                    "area_kind": "footprint_area",
                    "width_um": max_x - min_x,
                    "height_um": max_y - min_y,
                }
            )
    area = pd.DataFrame(area_rows)

    boundary_rows = []
    reference_labels = {d1a_label, d1b_label}
    for (sample, roi_name), frame in inside.groupby(
        ["sample", "primary_roi"], observed=True
    ):
        reference = frame.loc[frame["cell_type"].isin(reference_labels)]
        for threshold in neighborhood_distances:
            boundary_rows.append(
                {
                    "sample": str(sample),
                    "primary_roi": str(roi_name),
                    "neighborhood_distance_um": threshold,
                    "n_cells": len(frame),
                    "n_cells_near_boundary": int(
                        frame["distance_to_roi_boundary_um"].lt(threshold).sum()
                    ),
                    "fraction_cells_near_boundary": float(
                        frame["distance_to_roi_boundary_um"].lt(threshold).mean()
                    ),
                    "n_d1_reference_cells": len(reference),
                    "n_d1_reference_near_boundary": int(
                        reference["distance_to_roi_boundary_um"].lt(threshold).sum()
                    ),
                    "fraction_d1_reference_near_boundary": (
                        float(
                            reference["distance_to_roi_boundary_um"]
                            .lt(threshold)
                            .mean()
                        )
                        if len(reference)
                        else np.nan
                    ),
                }
            )
    boundary = pd.DataFrame(boundary_rows)
    boundary_grid = pd.MultiIndex.from_product(
        [
            sample_order,
            PRIMARY_ROIS,
            neighborhood_distances,
        ],
        names=[
            "sample",
            "primary_roi",
            "neighborhood_distance_um",
        ],
    )
    if boundary.empty:
        boundary = pd.DataFrame(index=boundary_grid).reset_index()
    else:
        boundary = (
            boundary.set_index(
                [
                    "sample",
                    "primary_roi",
                    "neighborhood_distance_um",
                ]
            )
            .reindex(boundary_grid)
            .reset_index()
        )
        for column in (
            "n_cells",
            "n_cells_near_boundary",
            "n_d1_reference_cells",
            "n_d1_reference_near_boundary",
        ):
            boundary[column] = boundary[column].fillna(0).astype(int)

    tile_rows = []
    for _, row in area.iterrows():
        for tile_size in tile_sizes:
            n_x = int(np.floor(row["width_um"] / tile_size))
            n_y = int(np.floor(row["height_um"] / tile_size))
            tile_rows.append(
                {
                    **row.to_dict(),
                    "tile_size_um": tile_size,
                    "approx_tiles_x": n_x,
                    "approx_tiles_y": n_y,
                    "approx_tiles_total": n_x * n_y,
                    "tile_size_comparable_to_roi": (
                        tile_size >= min(row["width_um"], row["height_um"])
                    ),
                }
            )
    tile = pd.DataFrame(tile_rows)

    eligibility = (
        total_counts.merge(d1_counts, on=["sample", "primary_roi"])
        .merge(area, on=["sample", "primary_roi"])
    )
    eligibility["flag_zero_cells"] = eligibility["n_cells"].eq(0)
    eligibility["flag_zero_d1a"] = eligibility["n_d1a"].eq(0)
    eligibility["flag_zero_d1b"] = eligibility["n_d1b"].eq(0)
    eligibility["flag_low_d1a"] = eligibility["n_d1a"].lt(
        low_reference_threshold
    )
    eligibility["flag_low_d1b"] = eligibility["n_d1b"].lt(
        low_reference_threshold
    )
    eligibility["flag_small_roi_area"] = eligibility["area_um2"].lt(
        small_roi_area
    )
    max_distance = max(neighborhood_distances)
    max_boundary = boundary.loc[
        boundary["neighborhood_distance_um"].eq(max_distance),
        [
            "sample",
            "primary_roi",
            "fraction_d1_reference_near_boundary",
        ],
    ]
    eligibility = eligibility.merge(
        max_boundary, on=["sample", "primary_roi"], how="left"
    )
    eligibility["boundary_distance_for_flag_um"] = max_distance
    eligibility["flag_high_d1_boundary_fraction"] = (
        eligibility["fraction_d1_reference_near_boundary"]
        .fillna(1.0)
        .gt(high_boundary_fraction_threshold)
    )
    tile_flag = (
        tile.groupby(["sample", "primary_roi"], observed=True)[
            "tile_size_comparable_to_roi"
        ]
        .any()
        .rename("flag_any_tile_size_comparable_to_roi")
        .reset_index()
    )
    eligibility = eligibility.merge(
        tile_flag, on=["sample", "primary_roi"], how="left"
    )
    return {
        "roi_cell_counts_by_sample": total_counts,
        "roi_celltype_counts_by_sample": celltype_counts,
        "roi_d1_island_counts_by_sample": d1_counts,
        "roi_area_by_sample": area,
        "boundary_distance_summary": boundary,
        "tile_size_qc": tile,
        "crawdad_eligibility_qc": eligibility,
    }


def plot_counts_and_boundary(
    qc: dict[str, pd.DataFrame],
    d1a_label: str,
    d1b_label: str,
    plot_dir: Path,
) -> None:
    """Figure 6: counts and D1 boundary diagnostics."""
    totals = qc["roi_cell_counts_by_sample"]
    d1 = qc["roi_d1_island_counts_by_sample"]
    boundary = qc["boundary_distance_summary"]
    samples = totals["sample"].drop_duplicates().astype(str).tolist()
    maximum_boundary_distance = max(
        boundary["neighborhood_distance_um"].dropna().max(), 1.0
    )
    x = np.arange(len(samples))
    roi_offsets = np.linspace(-0.25, 0.25, len(PRIMARY_ROIS))
    width = 0.22
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), squeeze=False)
    for offset, roi_name in zip(roi_offsets, PRIMARY_ROIS):
        frame = totals.loc[totals["primary_roi"].eq(roi_name)].set_index("sample")
        axes[0, 0].bar(
            x + offset,
            frame.reindex(samples)["n_cells"],
            width=width,
            color=ROI_COLORS[roi_name],
            label=roi_name,
        )
        d1_frame = d1.loc[d1["primary_roi"].eq(roi_name)].set_index("sample")
        axes[0, 1].bar(
            x + offset,
            d1_frame.reindex(samples)["n_d1a"],
            width=width,
            color=ROI_COLORS[roi_name],
            label=roi_name,
        )
        axes[1, 0].bar(
            x + offset,
            d1_frame.reindex(samples)["n_d1b"],
            width=width,
            color=ROI_COLORS[roi_name],
            label=roi_name,
        )
    axes[0, 0].set_title("Total cells")
    axes[0, 1].set_title(f"{d1a_label} cells")
    axes[1, 0].set_title(f"{d1b_label} cells")
    for axis in (axes[0, 0], axes[0, 1], axes[1, 0]):
        axis.set_xticks(x)
        axis.set_xticklabels(samples, rotation=90, fontsize=7)
        axis.set_ylabel("Cell count")
    for roi_name in PRIMARY_ROIS:
        frame = boundary.loc[boundary["primary_roi"].eq(roi_name)]
        thresholds = sorted(
            frame["neighborhood_distance_um"].dropna().unique()
        )
        for threshold_index, threshold in enumerate(thresholds):
            values = (
                frame.loc[
                    frame["neighborhood_distance_um"].eq(threshold)
                ]
                .set_index("sample")
                .reindex(samples)["fraction_d1_reference_near_boundary"]
            )
            axes[1, 1].plot(
                x,
                values,
                marker="o" if threshold_index == 0 else "s",
                linestyle="--" if threshold_index == 0 else "-",
                linewidth=1.2,
                color=ROI_COLORS[roi_name],
                alpha=min(
                    1.0,
                    0.5 + threshold / (2 * maximum_boundary_distance),
                ),
                label=f"{roi_name}, {threshold:g} µm",
            )
    axes[1, 1].set_title("D1 reference cells near ROI boundary")
    axes[1, 1].set_xticks(x)
    axes[1, 1].set_xticklabels(samples, rotation=90, fontsize=7)
    axes[1, 1].set_ylabel("Fraction")
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].legend(
        fontsize=7,
        ncol=2,
        title="ROI, distance cutoff",
        title_fontsize=7,
    )
    handles = [
        Patch(color=ROI_COLORS[name], label=name) for name in PRIMARY_ROIS
    ]
    figure.legend(
        handles=handles,
        loc="upper center",
        ncol=3,
        frameon=False,
        bbox_to_anchor=(0.5, 0.965),
    )
    figure.suptitle("CRAWDAD ROI counts and boundary diagnostics", y=0.995)
    figure.tight_layout(rect=(0, 0, 1, 0.92))
    save_figure(figure, plot_dir, "figure06_counts_and_boundary_diagnostics")


def write_assignments_and_crawdad_inputs(
    assignments: pd.DataFrame,
    output_dirs: dict[str, Path],
) -> None:
    """Write all-cell assignments and sample × ROI CRAWDAD inputs."""
    assignments.to_csv(
        output_dirs["assignments"] / "all_cells_roi_assignment.csv.gz",
        index=False,
        compression="gzip",
    )
    assignments.to_parquet(
        output_dirs["assignments"] / "all_cells_roi_assignment.parquet",
        index=False,
    )
    sample_order = (
        assignments.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    for sample in sample_order:
        sample_dir = output_dirs["crawdad"] / safe_name(str(sample))
        sample_dir.mkdir(parents=True, exist_ok=True)
        for roi_name in PRIMARY_ROIS:
            frame = assignments.loc[
                assignments["sample"].astype(str).eq(sample)
                & assignments["primary_roi"].eq(roi_name)
            ]
            output = frame[
                [
                    "cell_id",
                    "x_aligned",
                    "y_aligned",
                    "cell_type",
                    "sample",
                    "primary_roi",
                    "distance_to_roi_boundary_um",
                ]
            ].rename(
                columns={
                    "x_aligned": "x",
                    "y_aligned": "y",
                    "cell_type": "celltype",
                    "primary_roi": "roi",
                }
            )
            output.to_csv(
                sample_dir / f"{roi_name}.csv.gz",
                index=False,
                compression="gzip",
            )


def write_validation_report(
    path: Path,
    assignments: pd.DataFrame,
    footprints: dict[str, Any],
    lateral_name: str,
    primary: dict[str, Any],
    context: dict[str, Any],
    overlap_qc: pd.DataFrame,
    qc: dict[str, pd.DataFrame] | None,
    outside_warning_threshold: float,
    preview_only: bool,
) -> None:
    """Write a concise text validation report."""
    lines = [
        "CRAWDAD ROI validation report",
        "==============================",
        f"Mode: {'preview-only' if preview_only else 'final'}",
        f"Original VHD ROIs: {len(footprints)}",
        f"Lateral ROI: {lateral_name}",
        f"Medial VHD ROIs: {', '.join(context['medial_roi_names'])}",
        f"Medial union area (um^2): {context['medial_union_raw'].area:.6f}",
        f"Medial–lateral overlap assigned to lateral (um^2): {context['overlap_assigned_lateral'].area:.6f}",
        f"Medial y min (um): {context['medial_y_min']:.6f}",
        f"Medial y max (um): {context['medial_y_max']:.6f}",
        "Split method: aligned_xenium_broad_msn_robust_range_midpoint",
        f"Split reference cell types: {', '.join(context['split_reference_celltypes'])}",
        f"Split reference cell count: {context['split_reference_cell_count']}",
        f"Broad MSN y 0.05th percentile (um): {context['y_low']:.6f}",
        f"Broad MSN y 99.95th percentile (um): {context['y_high']:.6f}",
        f"Xenium-derived global y split (um): {context['y_split']:.6f}",
        "Dorsal direction: greater_y",
        f"Dorsomedial area (um^2): {primary['dorsomedial'].area:.6f}",
        f"Ventromedial area (um^2): {primary['ventromedial'].area:.6f}",
        f"Missing medial area after split (um^2): {context['missing_medial_after_split'].area:.6f}",
        f"Extra medial area after split (um^2): {context['extra_medial_after_split'].area:.6f}",
        "",
        "Pairwise primary-ROI overlap:",
        overlap_qc.to_string(index=False),
    ]
    if not preview_only:
        counts = assignments["primary_roi"].value_counts(dropna=False)
        outside_fraction = float(
            assignments["primary_roi"].eq(OUTSIDE_ROI).mean()
        )
        lines.extend(
            [
                "",
                "Cell assignments:",
                counts.to_string(),
                f"Outside-global-ROI fraction: {outside_fraction:.6f}",
            ]
        )
        if outside_fraction > outside_warning_threshold:
            lines.append(
                "WARNING: outside-global-ROI fraction exceeds configured threshold."
            )
        if qc is not None:
            eligibility = qc["crawdad_eligibility_qc"]
            flag_columns = [
                column for column in eligibility if column.startswith("flag_")
            ]
            lines.extend(
                [
                    "",
                    "Eligibility flag totals:",
                    eligibility[flag_columns].sum().to_string(),
                ]
            )
    path.write_text("\n".join(lines) + "\n")


def software_versions() -> dict[str, str]:
    """Return versions needed to reproduce geometry and output behavior."""
    versions = {"python": sys.version.split()[0]}
    for package in (
        "numpy",
        "pandas",
        "matplotlib",
        "shapely",
        "PyYAML",
        "pyarrow",
    ):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not installed"
    return versions


def main(argv: Sequence[str] | None = None) -> int:
    args = parse_args(argv)
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper(), logging.INFO),
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    root = git_root()
    resolve_paths(args, root)
    output_dirs = prepare_output_dirs(args)

    LOGGER.info("Loading aligned Xenium coordinates: %s", args.xenium_coordinate_csv)
    xenium = load_xenium_coordinates(
        args.xenium_coordinate_csv, args.donor
    )
    plotting_cells = xenium.rename(
        columns={
            "x_second_pass": "x_aligned",
            "y_second_pass": "y_aligned",
        }
    )

    footprint_sources: list[dict[str, Any]]
    if args.vhd_roi_file:
        LOGGER.info("Loading precomputed VHD ROI polygons: %s", args.vhd_roi_file)
        footprints = load_polygon_file(
            args.vhd_roi_file, args.roi_name_column
        )
        concave_hulls: dict[str, Any] = {}
        footprint_sources = [
            {
                "roi_name": name,
                "source_file": str(args.vhd_roi_file.resolve()),
                "hull_method": "precomputed_polygon",
                "hull_ratio": np.nan,
                "geometry_repaired": False,
            }
            for name in footprints
        ]
        vhd_files: list[Path] = []
    else:
        vhd_files = discover_vhd_coordinate_files(args)
        LOGGER.info(
            "Constructing six VHD footprints from %d aligned coordinate tables",
            len(vhd_files),
        )
        footprints, concave_hulls, footprint_sources = construct_vhd_footprints(
            vhd_files,
            args.donor,
            args.hull_method,
            args.hull_ratio,
        )

    tissue_polygon = None
    if args.tissue_polygon:
        tissue_geometries = load_polygon_file(
            args.tissue_polygon, args.roi_name_column
        )
        tissue_polygon = unary_union(list(tissue_geometries.values()))
        tissue_polygon, _ = repair_polygonal(tissue_polygon)

    primary, context = construct_primary_rois(
        footprints,
        args.lateral_roi_name,
        tissue_polygon,
        xenium,
    )
    LOGGER.info(
        "Global broad-MSN Xenium split: n=%d, y_low=%.6f, "
        "y_high=%.6f, y_split=%.6f µm; dorsal=greater_y",
        context["split_reference_cell_count"],
        context["y_low"],
        context["y_high"],
        context["y_split"],
    )
    overlap_qc = validate_primary_geometry(
        primary,
        context,
        args.geometry_tolerance_um2,
    )
    boolean_operation_qc = build_boolean_operation_qc(primary, context)
    geometry_qc = build_geometry_qc(
        footprints, primary, context, tissue_polygon
    )
    footprint_source_frame = pd.DataFrame(footprint_sources)
    geometry_qc = geometry_qc.merge(
        footprint_source_frame,
        on="roi_name",
        how="left",
        suffixes=("", "_source"),
    )

    write_geojson(
        footprints,
        output_dirs["polygons"] / "original_six_vhd_rois.geojson",
    )
    write_geojson(
        primary,
        output_dirs["polygons"] / "primary_three_rois.geojson",
    )
    geometry_qc.to_csv(
        output_dirs["polygons"] / "roi_geometry_qc.csv", index=False
    )
    overlap_qc.to_csv(
        output_dirs["qc"] / "primary_roi_overlap_qc.csv", index=False
    )
    boolean_operation_qc.to_csv(
        output_dirs["qc"] / "roi_boolean_operation_qc.csv", index=False
    )

    resolved_config = {
        "donor": args.donor,
        "coordinate_unit": args.coordinate_unit,
        "xenium_coordinate_csv": str(args.xenium_coordinate_csv.resolve()),
        "vhd_coordinate_files": [str(path.resolve()) for path in vhd_files],
        "vhd_roi_file": (
            str(args.vhd_roi_file.resolve()) if args.vhd_roi_file else None
        ),
        "lateral_roi_name": args.lateral_roi_name,
        "hull_method": args.hull_method,
        "hull_ratio": (
            args.hull_ratio
            if args.hull_method in {"concave", "concave_rectangle"}
            else None
        ),
        "tissue_polygon": (
            str(args.tissue_polygon.resolve()) if args.tissue_polygon else None
        ),
        "split_method": "aligned_xenium_broad_msn_robust_range_midpoint",
        "split_quantile_low": context["split_quantile_low"],
        "split_quantile_high": context["split_quantile_high"],
        "split_reference_celltypes": context[
            "split_reference_celltypes"
        ],
        "split_reference_cell_count": context[
            "split_reference_cell_count"
        ],
        "medial_y_min": context["medial_y_min"],
        "medial_y_max": context["medial_y_max"],
        "y_low": context["y_low"],
        "y_high": context["y_high"],
        "y_split": context["y_split"],
        "dorsal_direction": context["dorsal_direction"],
        "geometry_tolerance_um2": args.geometry_tolerance_um2,
        "lateral_area_um2": float(primary["lateral"].area),
        "dorsomedial_area_um2": float(primary["dorsomedial"].area),
        "ventromedial_area_um2": float(primary["ventromedial"].area),
        "pairwise_overlap_area_um2": {
            f"{row.roi_1}__{row.roi_2}": float(row.intersection_area_um2)
            for row in overlap_qc.itertuples(index=False)
        },
        "missing_medial_area_um2": float(
            context["missing_medial_after_split"].area
        ),
        "extra_medial_area_um2": float(
            context["extra_medial_after_split"].area
        ),
        "d1a_label": args.d1a_label,
        "d1b_label": args.d1b_label,
        "neighborhood_distances_um": args.neighborhood_distances_um,
        "tile_sizes_um": args.tile_sizes_um,
        "low_reference_cell_threshold": args.low_reference_cell_threshold,
        "small_roi_area_um2": args.small_roi_area_um2,
        "high_boundary_fraction_threshold": (
            args.high_boundary_fraction_threshold
        ),
        "large_outside_fraction_warning": (
            args.large_outside_fraction_warning
        ),
        "preview_only": args.preview_only,
        "output_dir": str(args.output_dir.resolve()),
        "plot_dir": str(args.plot_dir.resolve()),
        "software_versions": software_versions(),
    }
    (output_dirs["root"] / "config_resolved.yaml").write_text(
        yaml.safe_dump(resolved_config, sort_keys=False)
    )
    (output_dirs["root"] / "metadata.json").write_text(
        json.dumps(resolved_config, indent=2) + "\n"
    )

    if concave_hulls:
        plot_original_footprints(
            plotting_cells,
            concave_hulls,
            args.lateral_roi_name,
            output_dirs["plots"],
            args.max_plot_cells,
            title="Concave VHD hulls before minimum rotated rectangles",
            stem="figure01_concave_hulls_before_rectangle",
        )
    plot_original_footprints(
        plotting_cells,
        footprints,
        args.lateral_roi_name,
        output_dirs["plots"],
        args.max_plot_cells,
        title=(
            "Minimum rotated VHD rectangles over Xenium"
            if args.hull_method == "concave_rectangle"
            else "Original aligned VHD footprints over Xenium"
        ),
    )
    plot_split_preview(
        plotting_cells,
        footprints,
        context,
        args.lateral_roi_name,
        output_dirs["plots"],
        args.max_plot_cells,
    )
    plot_roi_construction(
        plotting_cells,
        footprints,
        args.lateral_roi_name,
        context,
        primary,
        output_dirs["plots"],
        args.max_plot_cells,
    )

    if args.preview_only:
        write_validation_report(
            output_dirs["qc"] / "validation_report.txt",
            plotting_cells,
            footprints,
            args.lateral_roi_name,
            primary,
            context,
            overlap_qc,
            qc=None,
            outside_warning_threshold=args.large_outside_fraction_warning,
            preview_only=True,
        )
        LOGGER.info("Preview complete. Plots: %s", args.plot_dir)
        return 0

    LOGGER.info("Assigning %d Xenium cells to primary ROIs", len(xenium))
    assignments = assign_cells(xenium, footprints, primary, context)
    write_assignments_and_crawdad_inputs(assignments, output_dirs)
    qc = build_qc_tables(
        assignments,
        primary,
        args.d1a_label,
        args.d1b_label,
        args.neighborhood_distances_um,
        args.tile_sizes_um,
        args.low_reference_cell_threshold,
        args.small_roi_area_um2,
        args.high_boundary_fraction_threshold,
    )
    for name, frame in qc.items():
        frame.to_csv(output_dirs["qc"] / f"{name}.csv", index=False)

    plot_final_rois(
        assignments,
        primary,
        context,
        output_dirs["plots"],
        args.max_plot_cells,
    )
    plot_rois_by_slice(
        assignments,
        primary,
        context,
        output_dirs["plots"],
        args.max_plot_cells,
    )
    plot_d1_overlay(
        assignments,
        primary,
        args.d1a_label,
        args.d1b_label,
        output_dirs["plots"],
        args.max_plot_cells,
    )
    plot_d1_overlay_by_slice(
        assignments,
        primary,
        context,
        args.d1a_label,
        args.d1b_label,
        output_dirs["plots"],
        args.max_plot_cells,
    )
    plot_counts_and_boundary(
        qc,
        args.d1a_label,
        args.d1b_label,
        output_dirs["plots"],
    )
    write_validation_report(
        output_dirs["qc"] / "validation_report.txt",
        assignments,
        footprints,
        args.lateral_roi_name,
        primary,
        context,
        overlap_qc,
        qc,
        args.large_outside_fraction_warning,
        preview_only=False,
    )
    LOGGER.info("Processed-data outputs: %s", args.output_dir)
    LOGGER.info("Figures: %s", args.plot_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/10_Xenium_CRAWDAD
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# Preview the automatic global y split without final cell assignments:
# python 01_select_crawdad_rois.py \
#   --donor Br6660 \
#   --hull-method concave_rectangle \
#   --preview-only
#
# Run the final Br6660 assignment after reviewing the preview:
# python 01_select_crawdad_rois.py \
#   --donor Br6660 \
#   --hull-method concave_rectangle
#
# Equivalent launcher preview:
# bash run_01_select_crawdad_rois_Br6660.sh --preview-only
#
# Equivalent launcher final run:
# bash run_01_select_crawdad_rois_Br6660.sh
