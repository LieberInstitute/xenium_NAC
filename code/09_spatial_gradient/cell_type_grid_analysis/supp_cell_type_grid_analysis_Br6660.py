"""Generate the Br6660 supplementary cell-type grid plots.

This plotting-only script reproduces two existing figures with slices labelled
1--11 in increasing anatomical depth. Outputs are written to the ``supp_fig``
subfolder and therefore do not overwrite the original figures.

Run from anywhere inside the repository with:

    python code/09_spatial_gradient/cell_type_grid_analysis/\
        supp_cell_type_grid_analysis_Br6660.py
"""

from pathlib import Path
import subprocess

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import numpy as np
import pandas as pd


DONOR = "Br6660"
GRID_SIZE = 5
POINT_SIZE = 0.5
POINT_ALPHA = 0.6
ANNOTATION_THRESHOLD = 0.05

CELLTYPE_PALETTE = {
    "Excitatory": "#DB1C5F",
    "Microglia_A": "#0D87E4",
    "D1_Island_B": "#00FF0D",
    "Astro_A": "#FD00FD",
    "Fibroblast_A": "#FFA7E2",
    "Fibroblast_B": "#2AFECA",
    "DRD2_MSN": "#7AAA16",
    "DRD1_MSN": "#9400FF",
    "Astro_B": "#823526",
    "MSN_Oligo": "#F5DEC0",
    "Inh_PVALB": "#93E5FF",
    "OPC": "#7A1699",
    "Microglia_B": "#FF0DBC",
    "Inh_SST": "#C4B3FB",
    "Astrocyte_Oligo": "#F80D2A",
    "D1_Island_A": "#224B82",
    "CHAT": "#FBA475",
    "Ependymal": "#73690D",
    "Microglia_Oligo": "#B62A7A",
    "WM": "orange",
}


def find_git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def load_plot_data(git_root: Path) -> tuple[pd.DataFrame, np.ndarray]:
    """Load Br6660 annotations and final second-pass XY coordinates."""
    celltype_path = (
        git_root
        / "processed-data"
        / "05_Clustering"
        / "Banksy_cell_types.csv"
    )
    coordinate_path = (
        git_root
        / "processed-data"
        / "05_Xenium_alignment"
        / "module02_Spateo_second_pass_Br6660"
        / "Br6660_second_pass_coordinates.csv"
    )

    celltypes = pd.read_csv(celltype_path, usecols=["cell_id", "CellType"])
    celltypes["cell_id"] = celltypes["cell_id"].astype(str)
    if celltypes["cell_id"].duplicated().any():
        raise ValueError("Cell-type table contains duplicate cell IDs")

    coordinates = pd.read_csv(
        coordinate_path,
        usecols=[
            "donor",
            "sample",
            "cell_id",
            "x_second_pass",
            "y_second_pass",
            "z_height",
        ],
    )
    coordinates["cell_id"] = coordinates["cell_id"].astype(str)
    if coordinates["cell_id"].duplicated().any():
        raise ValueError("Second-pass coordinate table contains duplicate cell IDs")
    if set(coordinates["donor"].dropna().astype(str)) != {DONOR}:
        raise ValueError("Second-pass coordinate table contains an unexpected donor")

    celltype_map = celltypes.set_index("cell_id")["CellType"]
    coordinates["CellType"] = coordinates["cell_id"].map(celltype_map)
    obs = coordinates.rename(columns={"sample": "Sample"})
    required = [
        "Sample", "CellType", "x_second_pass", "y_second_pass", "z_height"
    ]
    if obs[required].isna().any().any():
        missing_counts = obs[required].isna().sum()
        raise ValueError(
            "Missing annotations or coordinates:\n"
            + missing_counts[missing_counts.gt(0)].to_string()
        )

    sample_depth = (
        obs.groupby("Sample", sort=False)["z_height"].mean().sort_values()
    )
    if len(sample_depth) != 11:
        raise ValueError(f"Expected 11 Br6660 slices; found {len(sample_depth)}")
    slice_number = {
        sample: number
        for number, sample in enumerate(sample_depth.index.astype(str), start=1)
    }
    obs["slice_number"] = obs["Sample"].astype(str).map(slice_number).astype(int)
    xy = obs[["x_second_pass", "y_second_pass"]].to_numpy(dtype=float)
    return obs, xy


def make_slice_grid_overlay(
    obs: pd.DataFrame, xy: np.ndarray, output_path: Path
) -> None:
    """Plot the 11 depth-ordered slices in a fixed 2-by-6 layout."""
    x, y = xy[:, 0], xy[:, 1]
    xbins = np.linspace(x.min(), x.max(), GRID_SIZE + 1)
    ybins = np.linspace(y.min(), y.max(), GRID_SIZE + 1)
    colors = obs["CellType"].astype(str).map(CELLTYPE_PALETTE)
    if colors.isna().any():
        unknown = sorted(obs.loc[colors.isna(), "CellType"].astype(str).unique())
        raise ValueError(f"Missing colors for cell types: {unknown}")

    fig, axes = plt.subplots(2, 6, figsize=(23, 8), squeeze=False)
    axes_flat = axes.ravel()
    for slice_number in range(1, 12):
        ax = axes_flat[slice_number - 1]
        mask = obs["slice_number"].to_numpy() == slice_number
        ax.scatter(
            x[mask],
            y[mask],
            s=POINT_SIZE,
            c=colors.to_numpy()[mask],
            alpha=POINT_ALPHA,
            linewidths=0,
        )
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_aspect("equal")
        ax.set_title(f"Slice {slice_number}", fontsize=12)
        ax.set_xticks([])
        ax.set_yticks([])
        for boundary in xbins:
            ax.axvline(boundary, color="black", linewidth=0.6, alpha=0.6)
        for boundary in ybins:
            ax.axhline(boundary, color="black", linewidth=0.6, alpha=0.6)

    axes_flat[-1].axis("off")
    handles = [
        Patch(facecolor=color, label=celltype)
        for celltype, color in CELLTYPE_PALETTE.items()
    ]
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.875, 0.5),
        frameon=False,
        fontsize=9,
    )
    fig.subplots_adjust(right=0.86, wspace=0.12, hspace=0.20)
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def make_xy_proportion_grid(
    obs: pd.DataFrame, xy: np.ndarray, output_path: Path
) -> None:
    """Plot cell-type proportions with depth-ordered slice numbers on x axes."""
    celltype_order = list(CELLTYPE_PALETTE)
    frame = pd.DataFrame(
        {
            "X": xy[:, 0],
            "Y": xy[:, 1],
            "Slice": obs["slice_number"].to_numpy(dtype=int),
            "CellType": obs["CellType"].astype(str).to_numpy(),
        }
    )
    xbins = np.linspace(frame["X"].min(), frame["X"].max(), GRID_SIZE + 1)
    ybins = np.linspace(frame["Y"].min(), frame["Y"].max(), GRID_SIZE + 1)
    frame["Xbin"] = np.clip(
        np.digitize(frame["X"], xbins, right=False) - 1, 0, GRID_SIZE - 1
    )
    frame["Ybin"] = np.clip(
        np.digitize(frame["Y"], ybins, right=False) - 1, 0, GRID_SIZE - 1
    )

    fig, axes = plt.subplots(GRID_SIZE, GRID_SIZE, figsize=(16, 14))
    fig.suptitle(
        "Cell-type Proportions (ML/DV grid; bars = slices)",
        y=0.98,
        fontsize=14,
    )
    handles = [
        Patch(facecolor=color, label=celltype)
        for celltype, color in CELLTYPE_PALETTE.items()
    ]
    fig.legend(
        handles=handles,
        loc="center left",
        bbox_to_anchor=(0.88, 0.5),
        frameon=False,
        fontsize=8,
    )

    for i in range(GRID_SIZE):
        for j in range(GRID_SIZE):
            ax = axes[GRID_SIZE - 1 - j, i]
            subset = frame[(frame["Xbin"] == i) & (frame["Ybin"] == j)]
            if subset.empty:
                ax.axis("off")
                continue

            counts = (
                subset.groupby(["Slice", "CellType"], observed=False)
                .size()
                .unstack(fill_value=0)
                .reindex(index=np.arange(1, 12), fill_value=0)
                .reindex(columns=celltype_order, fill_value=0)
            )
            proportions = counts.div(counts.sum(axis=1).replace(0, np.nan), axis=0)
            proportions = proportions.fillna(0.0)
            xlocs = np.arange(1, 12)
            bottom = np.zeros(11, dtype=float)
            for celltype in celltype_order:
                values = proportions[celltype].to_numpy(dtype=float)
                ax.bar(
                    xlocs,
                    values,
                    bottom=bottom,
                    color=CELLTYPE_PALETTE[celltype],
                )
                bottom += values

            # Match the original plot: label segments only when proportion > 0.05.
            for bar_index in range(len(proportions)):
                cumulative = 0.0
                for celltype in celltype_order:
                    value = float(proportions.iloc[bar_index][celltype])
                    if value > ANNOTATION_THRESHOLD:
                        ax.text(
                            xlocs[bar_index],
                            cumulative + value / 2,
                            f"{value:.2f}",
                            ha="center",
                            va="center",
                            fontsize=6,
                        )
                    cumulative += value

            ax.set_ylim(0, 1)
            ax.set_xticks(xlocs)
            ax.set_xticklabels([str(number) for number in xlocs], fontsize=8)
            ax.set_title(
                f"ML/DV Grid ({i + 1},{j + 1}) | n={len(subset)}",
                fontsize=9,
            )
            ax.tick_params(axis="y", labelsize=7)

    fig.supxlabel("Slice number (ordered by depth)", fontsize=11)
    fig.subplots_adjust(
        left=0.06,
        right=0.86,
        top=0.92,
        bottom=0.08,
        hspace=0.35,
        wspace=0.20,
    )
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    git_root = find_git_root()
    output_dir = (
        git_root
        / "plots"
        / "09_spatial_gradient"
        / "cell_type_grid_analysis"
        / "supp_fig"
    )
    output_dir.mkdir(parents=True, exist_ok=True)
    obs, xy = load_plot_data(git_root)

    overlay_path = output_dir / "Br6660_slice_grid_overlay.png"
    proportion_path = output_dir / "Br6660_cell_type_proportion_grid_xy.png"
    make_slice_grid_overlay(obs, xy, overlay_path)
    make_xy_proportion_grid(obs, xy, proportion_path)
    print(f"Saved: {overlay_path}")
    print(f"Saved: {proportion_path}")


if __name__ == "__main__":
    main()
