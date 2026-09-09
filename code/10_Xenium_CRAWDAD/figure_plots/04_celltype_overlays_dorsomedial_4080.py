#!/usr/bin/env python3
"""Draw the three saved-ROI cell-type overlays for the 4080-um slice."""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-figure-plots")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Polygon as MplPolygon
from shapely.geometry import shape


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = CODE_ROOT.parent
MODULE06_PATH = CODE_ROOT / "10_Xenium_CRAWDAD/06_celltype_overlays.py"
DEFAULT_INPUT_ROOT = (
    PROJECT_ROOT
    / "processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois"
)
DEFAULT_PLOT_DIR = (
    PROJECT_ROOT
    / "plots/10_Xenium_CRAWDAD/figure_plots"
    / "panel_D_celltype_overlays_Br6660_Nac10_4080_dorsomedial"
)
SAMPLE = "Br6660_Nac10_4080"
REGION = "dorsomedial"
DM_BORDER_COLOR = "#006D5B"
DM_BORDER_LINEWIDTH = 1.5
LEGEND_LAYOUTS = {
    "fibroblast_microglia": (2, 2),
    "msn_island_astro": (3, 2),
    "astro_wm_island_oligo": (2, 2),
}


def load_module06():
    spec = importlib.util.spec_from_file_location(
        "module06_celltype_overlays", MODULE06_PATH
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load {MODULE06_PATH}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--background-color", default="#AFAFAF")
    parser.add_argument("--background-size", type=float, default=0.55)
    parser.add_argument("--highlight-size", type=float, default=2.2)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument(
        "--legend-only",
        action="store_true",
        help="Generate only the three standalone cell-type legends.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.background_size <= 0 or args.highlight_size <= 0:
        parser.error("Point sizes must be positive")
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    return args


def load_dm_roi(input_root: Path):
    polygon_path = input_root / "roi_polygons/primary_three_rois.geojson"
    if not polygon_path.is_file():
        raise FileNotFoundError(polygon_path)
    feature_collection = json.loads(polygon_path.read_text())
    rois = {
        feature["properties"]["roi_name"]: shape(feature["geometry"])
        for feature in feature_collection["features"]
    }
    if REGION not in rois:
        raise ValueError(f"Missing {REGION!r} ROI in {polygon_path}")
    return rois[REGION]


def polygon_components(geometry):
    if geometry.geom_type == "Polygon":
        yield geometry
    elif geometry.geom_type == "MultiPolygon":
        yield from geometry.geoms
    else:
        raise TypeError(f"Expected Polygon or MultiPolygon, got {geometry.geom_type}")


def shared_limits(cells, dm_roi) -> tuple[float, float, float, float]:
    roi_min_x, roi_min_y, roi_max_x, roi_max_y = dm_roi.bounds
    min_x = min(float(cells["x"].min()), roi_min_x)
    max_x = max(float(cells["x"].max()), roi_max_x)
    min_y = min(float(cells["y"].min()), roi_min_y)
    max_y = max(float(cells["y"].max()), roi_max_y)
    x_span = max(max_x - min_x, 1.0)
    y_span = max(max_y - min_y, 1.0)
    return (
        min_x - 0.02 * x_span,
        max_x + 0.02 * x_span,
        min_y - 0.02 * y_span,
        max_y + 0.02 * y_span,
    )


def draw_overlay(
    cells,
    selected_types: tuple[str, ...],
    target: Path,
    limits: tuple[float, float, float, float],
    args: argparse.Namespace,
    module06,
    dm_roi,
) -> None:
    absent = sorted(set(selected_types) - set(cells["celltype"].unique()))
    if absent:
        raise ValueError(f"Cell types absent from {SAMPLE}/{REGION}: {absent}")

    highlighted = cells["celltype"].isin(selected_types)
    background = cells.loc[~highlighted]
    figure, axis = plt.subplots(figsize=(9, 9))
    axis.scatter(
        background["x"],
        background["y"],
        s=args.background_size,
        c=args.background_color,
        alpha=0.34,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )
    for celltype in selected_types:
        selected = cells.loc[cells["celltype"].eq(celltype)]
        axis.scatter(
            selected["x"],
            selected["y"],
            s=args.highlight_size,
            c=module06.CELLTYPE_COLORS[celltype],
            alpha=0.9,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )

    for polygon in polygon_components(dm_roi):
        axis.add_patch(
            MplPolygon(
                polygon.exterior.coords,
                closed=True,
                facecolor="none",
                edgecolor=DM_BORDER_COLOR,
                linewidth=DM_BORDER_LINEWIDTH,
                zorder=3,
            )
        )

    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_aspect("equal")
    axis.set_axis_off()
    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    module06.atomic_save(figure, target, args.dpi)


def legend_handle_order(
    selected_types: tuple[str, ...],
    nrows: int,
    ncols: int,
) -> tuple[str, ...]:
    """Order handles so the rendered legend reads across rows.

    Matplotlib fills multi-column legends down each column. Constructing the
    requested row-by-column grid first and then flattening it column-wise
    preserves the cell-type order used by the corresponding overlay.
    """
    if len(selected_types) > nrows * ncols:
        raise ValueError(
            f"A {nrows}x{ncols} legend cannot hold "
            f"{len(selected_types)} cell types"
        )
    grid = np.full((nrows, ncols), None, dtype=object)
    grid.flat[: len(selected_types)] = selected_types
    return tuple(value for value in grid.T.flat if value is not None)


def draw_legend_only(
    selected_types: tuple[str, ...],
    target: Path,
    layout: tuple[int, int],
    args: argparse.Namespace,
    module06,
) -> None:
    """Save an alphabetized standalone legend for one overlay group."""
    nrows, ncols = layout
    alphabetical_types = tuple(sorted(selected_types, key=str.casefold))
    ordered_types = legend_handle_order(alphabetical_types, nrows, ncols)
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="",
            markerfacecolor=module06.CELLTYPE_COLORS[celltype],
            markeredgecolor="none",
            markersize=10,
            label=celltype,
        )
        for celltype in ordered_types
    ]
    figure = plt.figure(figsize=(3.0 * ncols, 0.85 * nrows))
    figure.legend(
        handles=handles,
        loc="center",
        ncol=ncols,
        frameon=False,
        fontsize=14,
        handletextpad=0.45,
        columnspacing=1.2,
        labelspacing=0.35,
    )
    module06.atomic_save(figure, target, args.dpi)


def main() -> int:
    args = parse_args()
    module06 = load_module06()
    groups = module06.DEFAULT_GROUPS
    targets = [
        args.plot_dir / f"overlay{index:02d}_{module06.safe_token(name)}.png"
        for index, (name, _) in enumerate(groups, start=1)
    ]
    legend_targets = [
        args.plot_dir
        / f"overlay{index:02d}_{module06.safe_token(name)}_legend.png"
        for index, (name, _) in enumerate(groups, start=1)
    ]
    requested_targets = (
        legend_targets if args.legend_only else targets + legend_targets
    )
    existing = [target for target in requested_targets if target.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            f"{len(existing)} output plots already exist; use --overwrite"
        )

    args.plot_dir.mkdir(parents=True, exist_ok=True)
    for legend_target, (name, selected_types) in zip(legend_targets, groups):
        try:
            layout = LEGEND_LAYOUTS[name]
        except KeyError as error:
            raise ValueError(
                f"No standalone legend layout defined for {name}"
            ) from error
        draw_legend_only(
            selected_types,
            legend_target,
            layout,
            args,
            module06,
        )

    if args.legend_only:
        print(f"Saved {len(legend_targets)} legends under {args.plot_dir}")
        return 0

    cells, input_path = module06.read_cells(
        args.input_root, [SAMPLE], [REGION]
    )
    cells = cells.loc[
        cells["sample"].eq(SAMPLE) & cells["roi"].eq(REGION)
    ].copy()
    dm_roi = load_dm_roi(args.input_root)
    limits = shared_limits(cells, dm_roi)
    for target, (_, selected_types) in zip(targets, groups):
        draw_overlay(
            cells,
            selected_types,
            target,
            limits,
            args,
            module06,
            dm_roi,
        )

    print(f"Read {len(cells):,} ROI cells from {input_path}")
    print(
        f"Saved {len(targets)} overlays and {len(legend_targets)} legends "
        f"under {args.plot_dir}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# python 10_Xenium_CRAWDAD/figure_plots/04_celltype_overlays_dorsomedial_4080.py
# Add --overwrite when replacing the existing plots.
#
# Generate only the three standalone legends:
# python 10_Xenium_CRAWDAD/figure_plots/04_celltype_overlays_dorsomedial_4080.py \
#     --legend-only \
#     --overwrite
