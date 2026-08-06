#!/usr/bin/env python3
"""Draw the manuscript version of the Module 01 ROI construction figure."""

from __future__ import annotations

import argparse
import importlib.util
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
from shapely.geometry import shape

from manuscript_style import (
    ROI_COLORS,
    ROI_HALO_COLOR,
    ROI_HALO_LINEWIDTH,
    ROI_LINEWIDTH,
)


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = CODE_ROOT.parent
DEFAULT_INPUT_DIR = (
    PROJECT_ROOT
    / "processed-data/10_Xenium_CRAWDAD/module01_select_crawdad_rois"
)
DEFAULT_PLOT_DIR = PROJECT_ROOT / "plots/10_Xenium_CRAWDAD/figure_plots"
CELL_COLOR = "#A65AA3"
VHD_COLOR = "#4D4D4D"


def load_module01_plotting():
    """Load the established Module 01 plotting helpers and color palette."""
    module_path = CODE_ROOT / "10_Xenium_CRAWDAD/01_select_crawdad_rois.py"
    spec = importlib.util.spec_from_file_location(
        "module01_select_crawdad_rois", module_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load plotting helpers from {module_path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Redraw Module 01 Figure 2 from saved cell assignments and final "
            "ROI polygons without rerunning ROI construction."
        )
    )
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--point-size", type=float, default=0.08)
    parser.add_argument("--point-alpha", type=float, default=0.35)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    if args.point_size <= 0:
        parser.error("--point-size must be positive")
    if not 0 < args.point_alpha <= 1:
        parser.error("--point-alpha must be in (0, 1]")
    return args


def load_saved_results(input_dir: Path):
    assignment_path = (
        input_dir / "cell_assignments/all_cells_roi_assignment.parquet"
    )
    polygon_path = input_dir / "roi_polygons/primary_three_rois.geojson"
    footprint_path = input_dir / "roi_polygons/original_six_vhd_rois.geojson"
    metadata_path = input_dir / "metadata.json"
    for path in (
        assignment_path,
        polygon_path,
        footprint_path,
        metadata_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    cells = pd.read_parquet(
        assignment_path, columns=["sample", "x_aligned", "y_aligned"]
    )
    if cells[["x_aligned", "y_aligned"]].isna().any().any():
        raise ValueError("Aligned cell coordinates contain missing values")
    if cells["sample"].nunique() != 11:
        raise ValueError(
            f"Expected 11 Xenium slices, found {cells['sample'].nunique()}"
        )

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

    footprint_collection = json.loads(footprint_path.read_text())
    footprints = {
        feature["properties"]["roi_name"]: shape(feature["geometry"])
        for feature in footprint_collection["features"]
    }
    if len(footprints) != 6:
        raise ValueError(f"Expected six VHD footprints, found {len(footprints)}")

    metadata = json.loads(metadata_path.read_text())
    return cells, footprints, rois, float(metadata["y_split"])


def draw_figure(
    cells: pd.DataFrame,
    footprints: dict,
    rois: dict,
    y_split: float,
    output_dir: Path,
    dpi: int,
    point_size: float,
    point_alpha: float,
    overwrite: bool,
) -> Path:
    module01 = load_module01_plotting()
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / "figure02_roi_construction.png"
    if output_path.exists() and not overwrite:
        raise FileExistsError(
            f"Output exists: {output_path}; use --overwrite to replace it"
        )

    figure, axis = plt.subplots(figsize=(9, 9))
    axis.scatter(
        cells["x_aligned"],
        cells["y_aligned"],
        s=point_size,
        c=CELL_COLOR,
        alpha=point_alpha,
        linewidths=0,
        rasterized=True,
        zorder=1,
    )

    for geometry in footprints.values():
        module01.draw_geometry(
            axis,
            geometry,
            edgecolor=VHD_COLOR,
            linewidth=0.8,
            zorder=2,
        )

    for roi_name in ("lateral", "dorsomedial", "ventromedial"):
        module01.draw_geometry(
            axis,
            rois[roi_name],
            edgecolor=ROI_HALO_COLOR,
            linewidth=ROI_HALO_LINEWIDTH,
            zorder=3.5,
        )
        module01.draw_geometry(
            axis,
            rois[roi_name],
            edgecolor=ROI_COLORS[roi_name],
            linewidth=ROI_LINEWIDTH,
            zorder=4,
        )

    axis.axhline(
        y_split,
        color="black",
        linestyle="--",
        linewidth=2,
        zorder=3,
    )

    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            color=CELL_COLOR,
            markersize=4,
            label="All Xenium cells (11 slices)",
        ),
        Line2D(
            [0],
            [0],
            color=VHD_COLOR,
            lw=0.8,
            label="Six original VHD footprints",
        ),
        Line2D(
            [0],
            [0],
            color=ROI_COLORS["lateral"],
            lw=ROI_LINEWIDTH,
            label="Lateral",
        ),
        Line2D(
            [0],
            [0],
            color=ROI_COLORS["dorsomedial"],
            lw=ROI_LINEWIDTH,
            label="Dorsomedial",
        ),
        Line2D(
            [0],
            [0],
            color=ROI_COLORS["ventromedial"],
            lw=ROI_LINEWIDTH,
            label="Ventromedial",
        ),
        Line2D(
            [0],
            [0],
            color="black",
            ls="--",
            lw=2,
            label="Dorsal–ventral split",
        ),
    ]
    axis.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
        fontsize=9,
    )
    axis.set_aspect("equal")
    axis.set_axis_off()

    module01.save_figure(figure, output_dir, output_path.stem, dpi=dpi)
    return output_path


def main() -> int:
    args = parse_args()
    cells, footprints, rois, y_split = load_saved_results(args.input_dir)
    output_path = draw_figure(
        cells,
        footprints,
        rois,
        y_split,
        args.plot_dir,
        args.dpi,
        args.point_size,
        args.point_alpha,
        args.overwrite,
    )
    print(f"Read {len(cells):,} cells from {cells['sample'].nunique()} slices")
    print(f"Saved {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# python 10_Xenium_CRAWDAD/figure_plots/01_figure02_roi_construction.py
# python 10_Xenium_CRAWDAD/figure_plots/01_figure02_roi_construction.py --overwrite
