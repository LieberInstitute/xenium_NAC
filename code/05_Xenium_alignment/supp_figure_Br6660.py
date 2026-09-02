#!/usr/bin/env python3
"""Generate Br6660 supplementary Xenium-alignment figures.

This plotting-only script creates:

1. Unaligned and final Spateo consecutive-slice overlays in 2-by-5 layouts.
2. A side-by-side PCC/MI comparison of unaligned and final Spateo scores.
3. Three 11-by-1 unaligned selected-cell-type slice figures and three compact
   versions containing only slices 3, 7, 8, 9, and 11.

The unaligned coordinates follow Module 00 exactly: raw spatial coordinates in
micrometers with Y reflected independently within each slice. Final coordinates
and similarity scores are read from existing Module 02 and Module 04 outputs;
no alignment or feature-similarity calculation is performed here.
"""

from __future__ import annotations

import argparse
import os
import subprocess
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import anndata as ad
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd


DONOR = "Br6660"
SAMPLE_ORDER = (
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

REFERENCE_COLOR = "#4C78A8"
MOVING_COLOR = "#F58518"
FINAL_COLOR = "#54A24B"
POINT_SIZE = 0.2
POINT_ALPHA = 0.45
CELLTYPE_POINT_SIZE = 0.6
CELLTYPE_POINT_ALPHA = 0.8

CELLTYPE_PALETTE = {
    "Excitatory": "#DB1C5F",
    "D1_Island_B": "#00FF0D",
    "DRD2_MSN": "#7AAA16",
    "DRD1_MSN": "#9400FF",
    "D1_Island_A": "#224B82",
    "WM": "orange",
}

CELLTYPE_FIGURES = (
    (
        ("D1_Island_A", "D1_Island_B", "WM", "Excitatory"),
        "Br6660_unaligned_D1_Island_A_B_WM_Excitatory.png",
    ),
    (
        ("D1_Island_A", "D1_Island_B", "DRD1_MSN", "DRD2_MSN"),
        "Br6660_unaligned_D1_Island_A_B_DRD1_MSN_DRD2_MSN.png",
    ),
    (
        (
            "D1_Island_A",
            "D1_Island_B",
            "WM",
            "Excitatory",
            "DRD1_MSN",
            "DRD2_MSN",
        ),
        "Br6660_unaligned_D1_Island_A_B_WM_Excitatory_DRD1_MSN_DRD2_MSN.png",
    ),
)

SELECTED_SLICE_NUMBERS = (3, 7, 8, 9, 11)

METRICS = (
    ("pcc_all", "PCC", "PCC (all multiscale median)"),
    ("mi_all", "MI", "MI (all multiscale median)"),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of existing supplementary PNG files.",
    )
    return parser.parse_args()


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def require_output(path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {path}")


def load_unaligned_coordinates(
    root: Path,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Load raw Br6660 XY/cell types and reproduce Module 00's Y flip."""
    input_path = (
        root
        / "processed-data"
        / "02_build_spe"
        / "h5ad"
        / "spe_NormCounts_nucleus_normcounts.h5ad"
    )
    adata = ad.read_h5ad(input_path, backed="r")
    try:
        donor_mask = adata.obs["Donor"].astype(str).eq(DONOR).to_numpy()
        donor_samples = adata.obs.loc[donor_mask, "Sample"].astype(str).to_numpy()
        donor_cell_ids = (
            adata.obs.loc[donor_mask, "cell_id"].astype(str).to_numpy()
        )
        donor_xy = np.asarray(adata.obsm["spatial"])[donor_mask, :2].astype(
            float, copy=True
        )
    finally:
        if getattr(adata, "file", None) is not None:
            adata.file.close()

    observed = set(donor_samples)
    missing = [sample for sample in SAMPLE_ORDER if sample not in observed]
    if missing:
        raise ValueError(f"Missing Br6660 samples in source h5ad: {missing}")

    annotation_path = (
        root
        / "processed-data"
        / "05_Clustering"
        / "Banksy_cell_types.csv"
    )
    annotations = pd.read_csv(annotation_path, usecols=["cell_id", "CellType"])
    annotations["cell_id"] = annotations["cell_id"].astype(str)
    if annotations["cell_id"].duplicated().any():
        raise ValueError("Cell-type annotation table contains duplicate cell IDs")
    celltype_map = annotations.set_index("cell_id")["CellType"]
    donor_celltypes = pd.Series(donor_cell_ids).map(celltype_map)
    if donor_celltypes.isna().any():
        examples = donor_cell_ids[donor_celltypes.isna().to_numpy()][:10]
        raise ValueError(
            "Missing cell-type annotations for Br6660 cells; examples="
            f"{examples.tolist()}"
        )
    donor_celltypes = donor_celltypes.astype(str).to_numpy()

    coordinates: dict[str, np.ndarray] = {}
    celltypes: dict[str, np.ndarray] = {}
    for sample in SAMPLE_ORDER:
        sample_mask = donor_samples == sample
        xy = donor_xy[sample_mask].copy()
        if xy.size == 0 or not np.isfinite(xy).all():
            raise ValueError(f"Invalid unaligned coordinates for {sample}")
        y = xy[:, 1].copy()
        xy[:, 1] = y.max() + y.min() - y
        coordinates[sample] = xy
        celltypes[sample] = donor_celltypes[sample_mask]
    return coordinates, celltypes


def load_final_coordinates(root: Path) -> dict[str, np.ndarray]:
    """Load final second-pass coordinates in anatomical depth order."""
    input_path = (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / "module02_Spateo_second_pass_Br6660"
        / "Br6660_second_pass_coordinates.csv"
    )
    frame = pd.read_csv(
        input_path,
        usecols=["donor", "sample", "cell_id", "x_second_pass", "y_second_pass"],
    )
    if frame["cell_id"].astype(str).duplicated().any():
        raise ValueError("Final coordinate CSV contains duplicate cell IDs")
    if set(frame["donor"].dropna().astype(str)) != {DONOR}:
        raise ValueError("Final coordinate CSV contains an unexpected donor")

    observed = set(frame["sample"].dropna().astype(str))
    missing = [sample for sample in SAMPLE_ORDER if sample not in observed]
    if missing:
        raise ValueError(f"Missing Br6660 samples in final coordinate CSV: {missing}")

    coordinates = {}
    for sample in SAMPLE_ORDER:
        selected = frame.loc[
            frame["sample"].astype(str).eq(sample),
            ["x_second_pass", "y_second_pass"],
        ]
        xy = selected.to_numpy(dtype=float)
        if xy.size == 0 or not np.isfinite(xy).all():
            raise ValueError(f"Invalid final coordinates for {sample}")
        coordinates[sample] = xy
    return coordinates


def plot_consecutive_overlays(
    coordinates: dict[str, np.ndarray],
    output_path: Path,
    figure_title: str,
    overwrite: bool,
) -> None:
    """Plot ten consecutive pairs in a 2-by-5 supplementary layout."""
    require_output(output_path, overwrite)
    figure, axes = plt.subplots(2, 5, figsize=(23, 9), squeeze=False)

    for pair_index, axis in enumerate(axes.ravel()):
        reference_number = pair_index + 1
        moving_number = pair_index + 2
        reference = coordinates[SAMPLE_ORDER[pair_index]]
        moving = coordinates[SAMPLE_ORDER[pair_index + 1]]
        reference_label = f"Slice {reference_number}"
        moving_label = f"Slice {moving_number}"

        axis.scatter(
            reference[:, 0],
            reference[:, 1],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            color=REFERENCE_COLOR,
            linewidths=0,
            label=reference_label,
            rasterized=True,
        )
        axis.scatter(
            moving[:, 0],
            moving[:, 1],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            color=MOVING_COLOR,
            linewidths=0,
            label=moving_label,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (µm)")
        axis.set_ylabel("y (µm)")
        axis.set_title(
            f"Slice {reference_number} → Slice {moving_number}", fontsize=10
        )
        axis.legend(
            loc="upper right",
            fontsize=7,
            markerscale=8,
            frameon=False,
        )

    figure.suptitle(figure_title, y=1.002, fontsize=14)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_unaligned_selected_celltypes(
    coordinates: dict[str, np.ndarray],
    celltypes_by_sample: dict[str, np.ndarray],
    selected_celltypes: tuple[str, ...],
    output_path: Path,
    overwrite: bool,
    slice_numbers: tuple[int, ...] = tuple(range(1, 12)),
    show_box: bool = True,
    show_slice_number: bool = True,
    vertical_spacing: float = 0.10,
    height_per_slice: float = 3.8,
    extra_gap_after_first: float = 0.0,
) -> None:
    """Plot selected cell types for requested slices from top to bottom."""
    require_output(output_path, overwrite)
    unknown = set(selected_celltypes) - set(CELLTYPE_PALETTE)
    if unknown:
        raise ValueError(f"Missing palette entries for: {sorted(unknown)}")
    invalid_slices = [
        number for number in slice_numbers
        if number < 1 or number > len(SAMPLE_ORDER)
    ]
    if invalid_slices:
        raise ValueError(f"Invalid slice numbers: {invalid_slices}")

    figure, axes = plt.subplots(
        len(slice_numbers),
        1,
        figsize=(7.0, height_per_slice * len(slice_numbers)),
        squeeze=False,
    )
    axes_flat = axes.ravel()

    for axis, slice_number in zip(axes_flat, slice_numbers):
        sample = SAMPLE_ORDER[slice_number - 1]
        xy = coordinates[sample]
        labels = celltypes_by_sample[sample]
        if len(xy) != len(labels):
            raise ValueError(f"Coordinate/label length mismatch for {sample}")

        for celltype in selected_celltypes:
            mask = labels == celltype
            axis.scatter(
                xy[mask, 0],
                xy[mask, 1],
                s=CELLTYPE_POINT_SIZE,
                alpha=CELLTYPE_POINT_ALPHA,
                color=CELLTYPE_PALETTE[celltype],
                linewidths=0,
                label=celltype,
                rasterized=True,
            )

        axis.set_aspect("equal")
        axis.set_xticks([])
        axis.set_yticks([])
        if not show_box:
            axis.set_frame_on(False)
            for spine in axis.spines.values():
                spine.set_visible(False)
        if show_slice_number:
            axis.text(
                0.02,
                0.96,
                f"Slice {slice_number}",
                transform=axis.transAxes,
                ha="left",
                va="top",
                fontsize=12,
                fontweight="bold",
            )

    legend_handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            linestyle="none",
            markersize=7,
            markerfacecolor=CELLTYPE_PALETTE[celltype],
            markeredgewidth=0,
            label=celltype,
        )
        for celltype in selected_celltypes
    ]
    figure.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.82, 0.995),
        frameon=False,
        fontsize=10,
    )
    figure.subplots_adjust(
        left=0.04,
        right=0.80,
        top=0.995 - extra_gap_after_first,
        bottom=0.01,
        hspace=vertical_spacing,
    )
    if extra_gap_after_first > 0:
        first_position = axes_flat[0].get_position()
        axes_flat[0].set_position(
            [
                first_position.x0,
                first_position.y0 + extra_gap_after_first,
                first_position.width,
                first_position.height,
            ]
        )
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def load_final_metric_comparison(root: Path) -> pd.DataFrame:
    input_path = (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / "module04_feature_similarity_comparison_Br6660"
        / "Br6660_unaligned_vs_final_spateo_multiscale_median_by_pair.csv"
    )
    frame = pd.read_csv(input_path)
    required = {
        "pair_number",
        "metric",
        "summary_statistic",
        "unaligned_score",
        "final_spateo_score",
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Metric comparison CSV is missing columns: {sorted(missing)}")
    selected = frame.loc[frame["metric"].isin(metric for metric, _, _ in METRICS)]
    if set(selected["metric"]) != {metric for metric, _, _ in METRICS}:
        raise ValueError("Metric comparison CSV does not contain both PCC and MI")
    if set(selected["summary_statistic"].astype(str)) != {"multiscale_median"}:
        raise ValueError("Expected multiscale_median comparison values")
    return selected.copy()


def plot_metric_comparison(
    comparison: pd.DataFrame,
    output_path: Path,
    overwrite: bool,
) -> None:
    """Plot only PCC and MI with one shared external legend."""
    require_output(output_path, overwrite)
    figure, axes = plt.subplots(1, 2, figsize=(15, 5.8), squeeze=False)
    legend_handles = None
    legend_labels = None

    for axis, (metric, panel_title, ylabel) in zip(axes.ravel(), METRICS):
        group = comparison.loc[comparison["metric"].eq(metric)].sort_values(
            "pair_number"
        )
        expected_pairs = np.arange(1, len(SAMPLE_ORDER))
        observed_pairs = group["pair_number"].to_numpy(dtype=int)
        if not np.array_equal(observed_pairs, expected_pairs):
            raise ValueError(
                f"Unexpected pair order for {metric}: {observed_pairs.tolist()}"
            )

        positions = np.arange(len(group))
        width = 0.38
        axis.bar(
            positions - width / 2,
            group["unaligned_score"],
            width,
            label="Unaligned",
            color=REFERENCE_COLOR,
        )
        axis.bar(
            positions + width / 2,
            group["final_spateo_score"],
            width,
            label="Final Spateo alignment",
            color=FINAL_COLOR,
        )
        axis.set_xticks(
            positions,
            [f"{number}–{number + 1}" for number in expected_pairs],
            fontsize=9,
        )
        axis.set_xlabel("Consecutive slice pair (ordered by depth)")
        axis.set_ylabel(ylabel)
        axis.set_title(panel_title)
        axis.axhline(0, color="black", linewidth=0.6, alpha=0.4)
        if legend_handles is None:
            legend_handles, legend_labels = axis.get_legend_handles_labels()

    figure.legend(
        legend_handles,
        legend_labels,
        loc="center left",
        bbox_to_anchor=(0.87, 0.5),
        frameon=False,
        fontsize=10,
    )
    figure.suptitle(
        "Unaligned vs final Spateo alignment: multiscale median",
        y=0.99,
        fontsize=14,
    )
    figure.subplots_adjust(left=0.08, right=0.84, bottom=0.16, top=0.86, wspace=0.28)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    root = git_root()
    output_dir = root / "plots" / "05_Xenium_alignment" / "supp_fig"

    unaligned_output = output_dir / "Br6660_overlay_slices_2d_unaligned.png"
    final_output = output_dir / "Br6660_overlay_slices_2d_second_spateo.png"
    metric_output = (
        output_dir
        / "Br6660_unaligned_vs_final_spateo_multiscale_median.png"
    )

    unaligned, unaligned_celltypes = load_unaligned_coordinates(root)
    final = load_final_coordinates(root)
    comparison = load_final_metric_comparison(root)

    plot_consecutive_overlays(
        unaligned,
        unaligned_output,
        "Unaligned consecutive-slice overlays",
        args.overwrite,
    )
    plot_consecutive_overlays(
        final,
        final_output,
        "Final Spateo consecutive-slice overlays",
        args.overwrite,
    )
    plot_metric_comparison(comparison, metric_output, args.overwrite)

    celltype_outputs = []
    for selected_celltypes, filename in CELLTYPE_FIGURES:
        output_path = output_dir / filename
        plot_unaligned_selected_celltypes(
            unaligned,
            unaligned_celltypes,
            selected_celltypes,
            output_path,
            args.overwrite,
        )
        celltype_outputs.append(output_path)

        is_drd1_drd2_compact = filename == (
            "Br6660_unaligned_D1_Island_A_B_DRD1_MSN_DRD2_MSN.png"
        )
        compact_slice_numbers = (
            (1, 3, 7, 8, 9, 11)
            if is_drd1_drd2_compact
            else SELECTED_SLICE_NUMBERS
        )
        compact_slice_tag = "_".join(map(str, compact_slice_numbers))
        compact_filename = (
            filename.removesuffix(".png")
            + f"_slices_{compact_slice_tag}.png"
        )
        compact_output_path = output_dir / compact_filename
        plot_unaligned_selected_celltypes(
            unaligned,
            unaligned_celltypes,
            selected_celltypes,
            compact_output_path,
            args.overwrite,
            slice_numbers=compact_slice_numbers,
            show_box=False,
            show_slice_number=not is_drd1_drd2_compact,
            vertical_spacing=-0.10 if is_drd1_drd2_compact else 0.02,
            height_per_slice=3.2 if is_drd1_drd2_compact else 3.8,
            extra_gap_after_first=0.02 if is_drd1_drd2_compact else 0.0,
        )
        celltype_outputs.append(compact_output_path)

    print(f"Saved: {unaligned_output}")
    print(f"Saved: {final_output}")
    print(f"Saved: {metric_output}")
    for output_path in celltype_outputs:
        print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/05_Xenium_alignment
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
# python supp_figure_Br6660.py --overwrite
