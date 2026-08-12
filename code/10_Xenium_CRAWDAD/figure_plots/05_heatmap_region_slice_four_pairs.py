#!/usr/bin/env python3
"""Redraw selected Module 05 relationship heatmaps by region and slice number."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-figure-plots")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = CODE_ROOT.parent
DEFAULT_INPUT = (
    PROJECT_ROOT
    / "processed-data/10_Xenium_CRAWDAD/module05_heatmap_region_depth"
    / "neighdist_50/relationship_points.parquet"
)
DEFAULT_PLOT_DIR = (
    PROJECT_ROOT
    / "plots/10_Xenium_CRAWDAD/figure_plots"
    / "heatmap_region_slice_neighdist_50_four_pairs"
)
REGIONS = (
    "lateral",
    "dorsomedial",
    "ventromedial",
    "outside_global_roi",
)
REGION_LABELS = {
    "lateral": "Lateral",
    "dorsomedial": "Dorsomedial",
    "ventromedial": "Ventromedial",
    "outside_global_roi": "Outside global ROI",
}
PAIRS = (
    ("D1_Island_A", "D1_Island_A"),
    ("D1_Island_A", "D1_Island_B"),
    ("D1_Island_B", "D1_Island_B"),
    ("D1_Island_B", "D1_Island_A"),
    ("Inh_PVALB", "D1_Island_B"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--z-color-limit", type=float, default=8.0)
    parser.add_argument("--minimum-point-area", type=float, default=180.0)
    parser.add_argument("--maximum-point-area", type=float, default=620.0)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.z_color_limit <= 0 or not np.isfinite(args.z_color_limit):
        parser.error("--z-color-limit must be positive and finite")
    if (
        args.minimum_point_area <= 0
        or args.maximum_point_area < args.minimum_point_area
    ):
        parser.error("Point areas must satisfy 0 < minimum <= maximum")
    if args.dpi < 72:
        parser.error("--dpi must be at least 72")
    return args


def point_areas(
    scales: np.ndarray,
    scale_levels: np.ndarray,
    minimum: float,
    maximum: float,
) -> np.ndarray:
    indices = np.array(
        [int(np.argmin(np.abs(scale_levels - value))) for value in scales]
    )
    if len(scale_levels) == 1:
        return np.full(len(scales), maximum)
    fraction = indices / (len(scale_levels) - 1)
    areas = maximum - fraction * (maximum - minimum)
    areas[indices == 0] *= 1.35
    return areas


def load_points(path: Path) -> tuple[pd.DataFrame, np.ndarray, float]:
    if not path.is_file():
        raise FileNotFoundError(path)
    columns = [
        "analysis_region",
        "sample",
        "reference",
        "neighbor",
        "depth",
        "neighborhood_distance_um",
        "reference_eligible",
        "first_significant_scale_um",
        "z_at_first_significant_scale",
        "ever_significant",
        "threshold_used",
    ]
    data = pd.read_parquet(path, columns=columns)
    distances = data["neighborhood_distance_um"].dropna().unique()
    if len(distances) != 1 or not np.isclose(float(distances[0]), 50):
        raise ValueError(f"Expected only neighDist 50 data, found {distances}")
    if set(data["analysis_region"].astype(str).unique()) != set(REGIONS):
        raise ValueError("Input does not contain the expected four regions")
    thresholds = data["threshold_used"].dropna().unique()
    if len(thresholds) != 1:
        raise ValueError("Expected one shared significance threshold")

    depth_table = (
        data[["sample", "depth"]]
        .drop_duplicates()
        .sort_values("depth")
        .reset_index(drop=True)
    )
    if len(depth_table) != 11 or depth_table["sample"].duplicated().any():
        raise ValueError("Expected exactly 11 unique depth-ordered slices")
    depth_table["slice_number"] = np.arange(1, 12)
    data = data.merge(
        depth_table, on=["sample", "depth"], how="left", validate="many_to_one"
    )
    scale_levels = np.arange(200.0, 1000.0 + 1, 100.0)
    observed_scales = data["first_significant_scale_um"].dropna().to_numpy(float)
    invalid = [
        value
        for value in observed_scales
        if not np.isclose(scale_levels, value, atol=1e-8, rtol=0).any()
    ]
    if invalid:
        raise ValueError(
            "First-significant scales fall outside 200–1000 µm: "
            f"{sorted(set(invalid))}"
        )
    return data, scale_levels, float(thresholds[0])


def save_figure(figure: plt.Figure, target: Path, dpi: int) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.stem}.tmp.{os.getpid()}.png")
    figure.savefig(
        temporary,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": Path(__file__).name},
    )
    plt.close(figure)
    temporary.replace(target)


def render_pair(
    data: pd.DataFrame,
    reference: str,
    neighbor: str,
    scale_levels: np.ndarray,
    target: Path,
    args: argparse.Namespace,
) -> int:
    pair = data.loc[
        data["reference"].eq(reference) & data["neighbor"].eq(neighbor)
    ].copy()
    if len(pair) != len(REGIONS) * 11:
        raise ValueError(
            f"Expected 44 rows for {reference} -> {neighbor}, found {len(pair)}"
        )
    significant = pair.loc[pair["ever_significant"].eq(True)].copy()
    x_positions = {region: index for index, region in enumerate(REGIONS)}

    figure, axis = plt.subplots(figsize=(6.1, 9.0))
    norm = Normalize(vmin=-args.z_color_limit, vmax=args.z_color_limit)
    if not significant.empty:
        areas = point_areas(
            significant["first_significant_scale_um"].to_numpy(float),
            scale_levels,
            args.minimum_point_area,
            args.maximum_point_area,
        )
        axis.scatter(
            significant["analysis_region"].map(x_positions),
            significant["slice_number"],
            c=np.clip(
                significant["z_at_first_significant_scale"].to_numpy(float),
                -args.z_color_limit,
                args.z_color_limit,
            ),
            s=areas,
            cmap="RdBu_r",
            norm=norm,
            edgecolors="#222222",
            linewidths=0.55,
            zorder=3,
        )

    axis.set_xticks(range(len(REGIONS)))
    axis.set_xticklabels(
        [REGION_LABELS[region] for region in REGIONS],
        rotation=35,
        ha="right",
        rotation_mode="anchor",
    )
    axis.set_yticks(range(1, 12))
    axis.set_yticklabels([str(value) for value in range(1, 12)])
    axis.set_xlim(-0.52, len(REGIONS) - 0.48)
    axis.set_ylim(11.52, 0.48)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("Region", fontsize=12, fontweight="bold")
    axis.set_ylabel("Slice number", fontsize=12, fontweight="bold")
    axis.set_title(
        f"{reference} → {neighbor}",
        fontsize=12,
        fontweight="bold",
        pad=12,
    )
    axis.grid(color="#d3d3d3", linewidth=0.8, zorder=0)
    axis.set_axisbelow(True)
    axis.tick_params(labelsize=10, pad=3)

    scalar = plt.cm.ScalarMappable(norm=norm, cmap="RdBu_r")
    scalar.set_array([])
    colorbar_axis = figure.add_axes([0.76, 0.53, 0.045, 0.34])
    colorbar = figure.colorbar(scalar, cax=colorbar_axis)
    colorbar.set_label("Mean Z at first significant scale", fontsize=9)
    legend_scales = np.array([200.0, 600.0, 1000.0])
    legend_areas = point_areas(
        legend_scales,
        scale_levels,
        args.minimum_point_area,
        args.maximum_point_area,
    )
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markerfacecolor="#9e9e9e",
            markeredgecolor="#222222",
            markeredgewidth=0.55,
            markersize=np.sqrt(area),
            label=f"{scale:g} µm",
        )
        for scale, area in zip(legend_scales, legend_areas)
    ]
    figure.legend(
        handles=handles,
        title="First significant scale",
        loc="upper center",
        bbox_to_anchor=(0.782, 0.47),
        frameon=False,
        fontsize=9,
        title_fontsize=9,
        labelspacing=2.1,
    )
    figure.subplots_adjust(left=0.16, right=0.70, bottom=0.14, top=0.98)
    save_figure(figure, target, args.dpi)
    return len(significant)


def main() -> int:
    args = parse_args()
    targets = [
        args.plot_dir / f"{reference}__to__{neighbor}.png"
        for reference, neighbor in PAIRS
    ]
    existing = [target for target in targets if target.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            f"{len(existing)} output plots already exist; use --overwrite"
        )

    data, scale_levels, _ = load_points(args.input)
    counts = []
    for (reference, neighbor), target in zip(PAIRS, targets):
        n_points = render_pair(
            data,
            reference,
            neighbor,
            scale_levels,
            target,
            args,
        )
        counts.append((reference, neighbor, n_points))

    print(f"Read saved Module 05 points from {args.input}")
    for reference, neighbor, count in counts:
        print(f"{reference} -> {neighbor}: {count}/44 significant positions")
    print(f"Saved {len(targets)} plots under {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# python 10_Xenium_CRAWDAD/figure_plots/05_heatmap_region_slice_four_pairs.py
# Add --overwrite when replacing the existing plots.
