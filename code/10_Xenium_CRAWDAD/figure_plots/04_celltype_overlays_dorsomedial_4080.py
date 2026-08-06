#!/usr/bin/env python3
"""Draw the three saved-ROI cell-type overlays for the 4080-um slice."""

from __future__ import annotations

import argparse
import importlib.util
import os
import sys
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-figure-plots")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt


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
    / "celltype_overlays_Br6660_Nac10_4080_dorsomedial"
)
SAMPLE = "Br6660_Nac10_4080"
REGION = "dorsomedial"


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
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.background_size <= 0 or args.highlight_size <= 0:
        parser.error("Point sizes must be positive")
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    return args


def shared_limits(cells) -> tuple[float, float, float, float]:
    x_span = max(float(cells["x"].max() - cells["x"].min()), 1.0)
    y_span = max(float(cells["y"].max() - cells["y"].min()), 1.0)
    return (
        float(cells["x"].min() - 0.02 * x_span),
        float(cells["x"].max() + 0.02 * x_span),
        float(cells["y"].min() - 0.02 * y_span),
        float(cells["y"].max() + 0.02 * y_span),
    )


def draw_overlay(
    cells,
    selected_types: tuple[str, ...],
    target: Path,
    limits: tuple[float, float, float, float],
    args: argparse.Namespace,
    module06,
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

    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_aspect("equal")
    axis.set_axis_off()
    figure.subplots_adjust(left=0, right=1, bottom=0, top=1)
    module06.atomic_save(figure, target, args.dpi)


def main() -> int:
    args = parse_args()
    module06 = load_module06()
    groups = module06.DEFAULT_GROUPS
    targets = [
        args.plot_dir / f"overlay{index:02d}_{module06.safe_token(name)}.png"
        for index, (name, _) in enumerate(groups, start=1)
    ]
    existing = [target for target in targets if target.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            f"{len(existing)} output plots already exist; use --overwrite"
        )

    cells, input_path = module06.read_cells(
        args.input_root, [SAMPLE], [REGION]
    )
    cells = cells.loc[
        cells["sample"].eq(SAMPLE) & cells["roi"].eq(REGION)
    ].copy()
    limits = shared_limits(cells)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    for target, (_, selected_types) in zip(targets, groups):
        draw_overlay(cells, selected_types, target, limits, args, module06)

    print(f"Read {len(cells):,} ROI cells from {input_path}")
    print(f"Saved {len(targets)} plots under {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# python 10_Xenium_CRAWDAD/figure_plots/04_celltype_overlays_dorsomedial_4080.py
# Add --overwrite when replacing the existing plots.
