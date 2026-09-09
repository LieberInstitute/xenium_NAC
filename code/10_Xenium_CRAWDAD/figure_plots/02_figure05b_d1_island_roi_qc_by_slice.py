#!/usr/bin/env python3
"""Draw depth-ordered cell-type QC panels from saved Module 01 results."""

from __future__ import annotations

import argparse
import importlib.util
import itertools
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib_crawdad_figure_plots")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import shape
from shapely.ops import unary_union

from manuscript_style import (
    ALL_CELLTYPE_ROI_COLOR,
    ALL_CELLTYPE_ROI_HALO_LINEWIDTH,
    ALL_CELLTYPE_ROI_LINEWIDTH,
    ROI_COLORS,
    ROI_HALO_COLOR,
    ROI_HALO_LINEWIDTH,
    ROI_LINEWIDTH,
)


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = CODE_ROOT.parent
MODULE_DIR = CODE_ROOT / "10_Xenium_CRAWDAD"
DEFAULT_INPUT_DIR = (
    PROJECT_ROOT
    / "processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois"
)
DEFAULT_PLOT_DIR = PROJECT_ROOT / "plots/10_Xenium_CRAWDAD/figure_plots"

PLOT_VERSIONS = {
    "figure05b_d1_island_ab": ("D1_Island_A", "D1_Island_B"),
    "figure05b_d1_island_ab_inh_pvalb": (
        "D1_Island_A",
        "D1_Island_B",
        "Inh_PVALB",
    ),
}
ALL_CELLTYPES_VERSION = "figure05b_all_celltypes"
PLOT_VERSION_NAMES = (*PLOT_VERSIONS, ALL_CELLTYPES_VERSION)
PANEL_B_DIRNAME = "panel_B_all_celltypes"
PANEL_F_DIRNAME = "panel_F_celltype_across_slices"
PANEL_B_DM_LINEWIDTH_SCALE = 1.6
SHARED_BOUNDARY_TOLERANCE_UM = 1.0
MIN_SHARED_BOUNDARY_LENGTH_UM = 10.0


def version_output_dir(plot_dir: Path, version: str) -> Path:
    """Route each plot version to its renamed manuscript-panel directory."""
    if version == ALL_CELLTYPES_VERSION:
        return plot_dir / PANEL_B_DIRNAME
    return plot_dir / PANEL_F_DIRNAME / version


def iter_line_geometries(geometry):
    """Yield LineStrings from a potentially multipart boundary geometry."""
    if geometry.is_empty:
        return
    if geometry.geom_type in {"LineString", "LinearRing"}:
        yield geometry
        return
    if hasattr(geometry, "geoms"):
        for part in geometry.geoms:
            yield from iter_line_geometries(part)


def draw_nonshared_roi_boundaries(
    axis, rois: dict, roi_colors: dict[str, str], roi_linewidth: float
) -> None:
    """Draw colored ROI outline portions not shared with another ROI."""
    for roi_name, roi in rois.items():
        nonshared = roi.boundary
        for other_name, other_roi in rois.items():
            if other_name == roi_name:
                continue
            nonshared = nonshared.difference(
                other_roi.boundary.buffer(SHARED_BOUNDARY_TOLERANCE_UM)
            )
        for line in iter_line_geometries(nonshared):
            x, y = line.xy
            axis.plot(
                x,
                y,
                color=roi_colors[roi_name],
                linewidth=roi_linewidth,
                linestyle="-",
                solid_capstyle="butt",
                solid_joinstyle="miter",
                zorder=3,
            )


def draw_shared_roi_boundaries(
    axis,
    rois: dict,
    module01,
    roi_colors: dict[str, str],
    roi_linewidth: float,
) -> None:
    """Split shared strokes across the interiors of adjacent ROIs."""
    for first_name, second_name in itertools.combinations(rois, 2):
        shared = rois[first_name].boundary.intersection(
            rois[second_name].boundary.buffer(
                SHARED_BOUNDARY_TOLERANCE_UM
            )
        )
        for line in iter_line_geometries(shared):
            if line.length < MIN_SHARED_BOUNDARY_LENGTH_UM:
                continue
            x, y = line.xy
            common = dict(
                linewidth=roi_linewidth,
                linestyle="-",
                solid_capstyle="butt",
            )
            axis.plot(
                x,
                y,
                color=roi_colors[second_name],
                zorder=4,
                **common,
            )
            for polygon in module01.polygon_components(rois[first_name]):
                clip_patch = MplPolygon(
                    polygon.exterior.coords,
                    closed=True,
                    transform=axis.transData,
                )
                artist, = axis.plot(
                    x,
                    y,
                    color=roi_colors[first_name],
                    zorder=4.1,
                    **common,
                )
                artist.set_clip_path(clip_patch)


def load_script_module(name: str, path: Path):
    """Load constants and plotting helpers from an existing module."""
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Draw 4-by-2, column-major D1 Island A/B, D1 Island A/B "
            "with Inh_PVALB, and all-cell-type panels from saved Module 01 "
            "assignments."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--point-size", type=float, default=0.45)
    parser.add_argument("--point-alpha", type=float, default=0.75)
    parser.add_argument(
        "--depth-range-um",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(1580.0, 5080.0),
    )
    parser.add_argument(
        "--versions",
        nargs="+",
        choices=PLOT_VERSION_NAMES,
        default=list(PLOT_VERSION_NAMES),
        help="Plot versions to generate; all three are generated by default.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    if args.point_size <= 0:
        parser.error("--point-size must be positive")
    if not 0 < args.point_alpha <= 1:
        parser.error("--point-alpha must be in (0, 1]")
    if args.depth_range_um[0] > args.depth_range_um[1]:
        parser.error("--depth-range-um MIN must be <= MAX")
    return args


def load_saved_results(input_dir: Path, selected_types: list[str]):
    assignment_path = (
        input_dir / "cell_assignments/all_cells_roi_assignment.parquet"
    )
    polygon_path = input_dir / "roi_polygons/primary_three_rois.geojson"
    metadata_path = input_dir / "metadata.json"
    for path in (assignment_path, polygon_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    cells = pd.read_parquet(
        assignment_path,
        columns=[
            "sample",
            "cell_type",
            "x_aligned",
            "y_aligned",
            "z_height",
        ],
        filters=[("cell_type", "in", selected_types)],
    )
    if cells[["x_aligned", "y_aligned", "z_height"]].isna().any().any():
        raise ValueError("Selected cells contain missing coordinates or depths")
    depth_counts = cells.groupby("sample", observed=True)["z_height"].nunique()
    if not depth_counts.eq(1).all():
        raise ValueError("Every sample must have exactly one z_height")
    sample_depths = cells.groupby("sample", observed=True)["z_height"].first()
    sample_depths = sample_depths.sort_values()
    if len(sample_depths) != 11:
        raise ValueError(f"Expected 11 Xenium slices, found {len(sample_depths)}")
    observed = set(cells["cell_type"].astype(str))
    missing_types = sorted(set(selected_types) - observed)
    if missing_types:
        raise ValueError(f"Cell types absent from assignments: {missing_types}")

    feature_collection = json.loads(polygon_path.read_text())
    rois = {
        feature["properties"]["roi_name"]: shape(feature["geometry"])
        for feature in feature_collection["features"]
    }
    expected_rois = {"lateral", "dorsomedial", "ventromedial"}
    if set(rois) != expected_rois:
        raise ValueError(
            f"Expected final ROIs {sorted(expected_rois)}, found {sorted(rois)}"
        )

    slice_table = sample_depths.rename("z_height").reset_index()
    slice_table["sample"] = slice_table["sample"].astype(str)
    slice_table["slice_number"] = range(1, len(slice_table) + 1)
    metadata = json.loads(metadata_path.read_text())
    return cells, slice_table, rois, float(metadata["y_split"])


def shared_limits(cells: pd.DataFrame, rois: dict) -> tuple[float, ...]:
    roi_bounds = unary_union(list(rois.values())).bounds
    min_x = min(float(cells["x_aligned"].min()), roi_bounds[0])
    min_y = min(float(cells["y_aligned"].min()), roi_bounds[1])
    max_x = max(float(cells["x_aligned"].max()), roi_bounds[2])
    max_y = max(float(cells["y_aligned"].max()), roi_bounds[3])
    x_margin = 0.03 * max(max_x - min_x, 1.0)
    y_margin = 0.03 * max(max_y - min_y, 1.0)
    return (
        min_x - x_margin,
        max_x + x_margin,
        min_y - y_margin,
        max_y + y_margin,
    )


def draw_panel(
    axis,
    frame: pd.DataFrame,
    celltypes: tuple[str, ...],
    rois: dict,
    y_split: float,
    limits: tuple[float, ...],
    args: argparse.Namespace,
    module01,
    palette: dict[str, str],
    roi_colors: dict[str, str],
    roi_halo_linewidth: float,
    roi_linewidth: float,
    show_split: bool,
    slice_label: str | None,
    dorsomedial_linewidth_scale: float = 1.0,
    split_shared_borders: bool = False,
) -> None:
    for celltype in celltypes:
        selected = frame.loc[frame["cell_type"].eq(celltype)]
        axis.scatter(
            selected["x_aligned"],
            selected["y_aligned"],
            s=args.point_size,
            color=palette[celltype],
            alpha=args.point_alpha,
            linewidths=0,
            rasterized=True,
            zorder=1,
        )
    # Draw the common white halo independently of the colored border logic.
    for roi_name in ("lateral", "ventromedial", "dorsomedial"):
        module01.draw_geometry(
            axis,
            rois[roi_name],
            edgecolor=ROI_HALO_COLOR,
            linewidth=roi_halo_linewidth,
            zorder=2.5,
        )
    if split_shared_borders:
        draw_nonshared_roi_boundaries(
            axis, rois, roi_colors, roi_linewidth
        )
        draw_shared_roi_boundaries(
            axis, rois, module01, roi_colors, roi_linewidth
        )
    else:
        # Panel B keeps its established black outlines with a thicker green
        # dorsomedial border drawn last at shared edges.
        for roi_name in ("lateral", "ventromedial", "dorsomedial"):
            module01.draw_geometry(
                axis,
                rois[roi_name],
                edgecolor=roi_colors[roi_name],
                linewidth=(
                    roi_linewidth * dorsomedial_linewidth_scale
                    if roi_name == "dorsomedial"
                    else roi_linewidth
                ),
                zorder=3,
            )
    if show_split:
        axis.axhline(
            y_split,
            color="black",
            linestyle="--",
            linewidth=1,
            zorder=2,
        )
    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_aspect("equal")
    axis.set_axis_off()
    if slice_label is not None:
        axis.text(
            0.02,
            0.97,
            slice_label,
            transform=axis.transAxes,
            ha="left",
            va="top",
            fontsize=20,
            weight="bold",
        )


def draw_multislice(
    cells: pd.DataFrame,
    slice_table: pd.DataFrame,
    celltypes: tuple[str, ...],
    rois: dict,
    y_split: float,
    limits: tuple[float, ...],
    output_dir: Path,
    stem: str,
    args: argparse.Namespace,
    module01,
    palette: dict[str, str],
    roi_colors: dict[str, str],
    roi_halo_linewidth: float,
    roi_linewidth: float,
    show_split: bool,
    show_slice_labels: bool,
    dorsomedial_linewidth_scale: float = 1.0,
    split_shared_borders: bool = False,
) -> Path:
    nrows = 4
    ncols = (len(slice_table) + nrows - 1) // nrows
    figure, axes = plt.subplots(
        nrows,
        ncols,
        figsize=(6.5 * ncols, 5.2 * nrows),
        squeeze=False,
    )
    column_major_axes = [
        axes[row_index, column_index]
        for column_index in range(ncols)
        for row_index in range(nrows)
    ]
    for axis, row in zip(
        column_major_axes, slice_table.itertuples(index=False)
    ):
        frame = cells.loc[cells["sample"].astype(str).eq(row.sample)]
        draw_panel(
            axis,
            frame,
            celltypes,
            rois,
            y_split,
            limits,
            args,
            module01,
            palette,
            roi_colors,
            roi_halo_linewidth,
            roi_linewidth,
            show_split,
            slice_label=(
                f"Slice {row.slice_number}" if show_slice_labels else None
            ),
            dorsomedial_linewidth_scale=dorsomedial_linewidth_scale,
            split_shared_borders=split_shared_borders,
        )
    for axis in column_major_axes[len(slice_table) :]:
        axis.set_axis_off()
    figure.subplots_adjust(hspace=0.02, wspace=0.02)
    module01.save_figure(figure, output_dir, stem, dpi=args.dpi)
    return output_dir / f"{stem}.png"


def draw_single_slice(
    frame: pd.DataFrame,
    celltypes: tuple[str, ...],
    rois: dict,
    y_split: float,
    limits: tuple[float, ...],
    output_dir: Path,
    stem: str,
    args: argparse.Namespace,
    module01,
    palette: dict[str, str],
    roi_colors: dict[str, str],
    roi_halo_linewidth: float,
    roi_linewidth: float,
    show_split: bool,
) -> Path:
    figure, axis = plt.subplots(figsize=(9, 9))
    draw_panel(
        axis,
        frame,
        celltypes,
        rois,
        y_split,
        limits,
        args,
        module01,
        palette,
        roi_colors,
        roi_halo_linewidth,
        roi_linewidth,
        show_split,
        slice_label=None,
    )
    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    module01.save_figure(figure, output_dir, stem, dpi=args.dpi)
    return output_dir / f"{stem}.png"


def draw_legend(
    celltypes: tuple[str, ...],
    output_dir: Path,
    dpi: int,
    module01,
    palette: dict[str, str],
    roi_colors: dict[str, str],
    roi_linewidth: float,
    show_split: bool,
) -> Path:
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=palette[celltype],
            markersize=6,
            label=celltype,
        )
        for celltype in celltypes
    ]
    handles.extend(
        Line2D(
            [0],
            [0],
            color=roi_colors[roi_name],
            lw=roi_linewidth,
            label=roi_name.capitalize(),
        )
        for roi_name in ("lateral", "dorsomedial", "ventromedial")
    )
    if show_split:
        handles.append(
            Line2D(
                [0],
                [0],
                color="black",
                linestyle="--",
                lw=1,
                label="Dorsal–ventral split",
            )
        )
    many_entries = len(handles) > 8
    figure = plt.figure(figsize=(9, 5) if many_entries else (7, 2.2))
    figure.legend(
        handles=handles,
        loc="center",
        frameon=False,
        fontsize=9,
        ncol=4 if many_entries else 3,
    )
    module01.save_figure(figure, output_dir, "legend", dpi=dpi)
    return output_dir / "legend.png"


def draw_celltype_legend_3cols(
    celltypes: tuple[str, ...],
    output_dir: Path,
    dpi: int,
    module01,
    palette: dict[str, str],
) -> Path:
    alphabetical = sorted(celltypes, key=str.casefold)
    ncols = 3
    nrows = (len(alphabetical) + ncols - 1) // ncols
    column_major = [
        alphabetical[row * ncols + column]
        for column in range(ncols)
        for row in range(nrows)
        if row * ncols + column < len(alphabetical)
    ]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=palette[celltype],
            markersize=7,
            label=celltype,
        )
        for celltype in column_major
    ]
    figure = plt.figure(figsize=(8, 6.5))
    figure.legend(
        handles=handles,
        loc="center",
        frameon=False,
        fontsize=10,
        ncol=ncols,
    )
    module01.save_figure(
        figure, output_dir, "legend_celltypes_3cols", dpi=dpi
    )
    return output_dir / "legend_celltypes_3cols.png"


def main() -> int:
    args = parse_args()
    module01 = load_script_module(
        "module01_select_crawdad_rois",
        MODULE_DIR / "01_select_crawdad_rois.py",
    )
    module06 = load_script_module(
        "module06_celltype_overlays",
        MODULE_DIR / "06_celltype_overlays.py",
    )
    plot_versions = dict(PLOT_VERSIONS)
    plot_versions[ALL_CELLTYPES_VERSION] = tuple(module06.CELLTYPE_COLORS)
    selected_versions = {
        name: plot_versions[name] for name in args.versions
    }
    selected_types = sorted(
        {
            celltype
            for values in selected_versions.values()
            for celltype in values
        }
    )
    cells, slice_table, rois, y_split = load_saved_results(
        args.input_dir, selected_types
    )
    minimum_depth, maximum_depth = args.depth_range_um
    selected_slices = slice_table.loc[
        slice_table["z_height"].between(minimum_depth, maximum_depth)
    ].copy()
    if selected_slices.empty:
        raise ValueError("No slices fall within --depth-range-um")
    selected_samples = set(selected_slices["sample"])
    selected_cells = cells.loc[
        cells["sample"].astype(str).isin(selected_samples)
    ].copy()
    limits = shared_limits(selected_cells, rois)
    depth_token = (
        f"{minimum_depth:g}_{maximum_depth:g}".replace(".", "p")
    )
    plot_stem = f"figure05b_depth{depth_token}"

    targets = [
        version_output_dir(args.plot_dir, version) / f"{plot_stem}.png"
        for version in selected_versions
    ]
    targets.extend(
        version_output_dir(args.plot_dir, version) / "legend.png"
        for version in selected_versions
    )
    if ALL_CELLTYPES_VERSION in selected_versions:
        targets.extend(
            [
                version_output_dir(args.plot_dir, ALL_CELLTYPES_VERSION)
                / "figure05b_slice8_all_celltypes.png",
                version_output_dir(args.plot_dir, ALL_CELLTYPES_VERSION)
                / "legend_celltypes_3cols.png",
            ]
        )
    existing = [path for path in targets if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            f"{len(existing)} output plots already exist; use --overwrite"
        )

    written = []
    for version, celltypes in selected_versions.items():
        output_dir = version_output_dir(args.plot_dir, version)
        output_dir.mkdir(parents=True, exist_ok=True)
        show_split = version != ALL_CELLTYPES_VERSION
        roi_colors = (
            {roi_name: ALL_CELLTYPE_ROI_COLOR for roi_name in ROI_COLORS}
            if version == ALL_CELLTYPES_VERSION
            else ROI_COLORS
        )
        roi_halo_linewidth = (
            ALL_CELLTYPE_ROI_HALO_LINEWIDTH
            if version == ALL_CELLTYPES_VERSION
            else ROI_HALO_LINEWIDTH
        )
        roi_linewidth = (
            ALL_CELLTYPE_ROI_LINEWIDTH
            if version == ALL_CELLTYPES_VERSION
            else ROI_LINEWIDTH
        )
        multislice_roi_colors = dict(roi_colors)
        if version == ALL_CELLTYPES_VERSION:
            # Panel B uses black lateral/ventromedial outlines and the same
            # deep green dorsomedial outline as panel_A_roi_construction.png.
            multislice_roi_colors["dorsomedial"] = ROI_COLORS["dorsomedial"]
        written.append(
            draw_multislice(
                selected_cells,
                selected_slices,
                celltypes,
                rois,
                y_split,
                limits,
                output_dir,
                plot_stem,
                args,
                module01,
                module06.CELLTYPE_COLORS,
                multislice_roi_colors,
                roi_halo_linewidth,
                roi_linewidth,
                show_split,
                show_slice_labels=version != ALL_CELLTYPES_VERSION,
                dorsomedial_linewidth_scale=(
                    PANEL_B_DM_LINEWIDTH_SCALE
                    if version == ALL_CELLTYPES_VERSION
                    else 1.0
                ),
                split_shared_borders=(
                    version != ALL_CELLTYPES_VERSION
                ),
            )
        )
        written.append(
            draw_legend(
                celltypes,
                output_dir,
                args.dpi,
                module01,
                module06.CELLTYPE_COLORS,
                roi_colors,
                roi_linewidth,
                show_split,
            )
        )
        if version == ALL_CELLTYPES_VERSION:
            slice8 = slice_table.loc[slice_table["slice_number"].eq(8)]
            if len(slice8) != 1:
                raise ValueError(
                    f"Expected one numeric Slice 8, found {len(slice8)}"
                )
            slice8_sample = str(slice8.iloc[0]["sample"])
            slice8_cells = cells.loc[
                cells["sample"].astype(str).eq(slice8_sample)
            ]
            written.append(
                draw_single_slice(
                    slice8_cells,
                    celltypes,
                    rois,
                    y_split,
                    limits,
                    output_dir,
                    "figure05b_slice8_all_celltypes",
                    args,
                    module01,
                    module06.CELLTYPE_COLORS,
                    roi_colors,
                    roi_halo_linewidth,
                    roi_linewidth,
                    show_split=False,
                )
            )
            written.append(
                draw_celltype_legend_3cols(
                    celltypes,
                    output_dir,
                    args.dpi,
                    module01,
                    module06.CELLTYPE_COLORS,
                )
            )

    slice_labels = ", ".join(
        f"Slice {number} ({depth:g} µm)"
        for number, depth in zip(
            selected_slices["slice_number"], selected_slices["z_height"]
        )
    )
    print(f"Selected {len(selected_cells):,} cells from: {slice_labels}")
    print(f"Saved {len(written)} plots under {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# python 10_Xenium_CRAWDAD/figure_plots/02_figure05b_d1_island_roi_qc_by_slice.py
# python 10_Xenium_CRAWDAD/figure_plots/02_figure05b_d1_island_roi_qc_by_slice.py --depth-range-um 1580 5080
# python 10_Xenium_CRAWDAD/figure_plots/02_figure05b_d1_island_roi_qc_by_slice.py --versions figure05b_all_celltypes
# python 10_Xenium_CRAWDAD/figure_plots/02_figure05b_d1_island_roi_qc_by_slice.py --overwrite
