#!/usr/bin/env python3
"""QC plots and 3D reconstructions for final Xenium–VisiumHD alignments.

Inputs
------
Xenium
    Final donor-level coordinates from
    ``processed-data/05_Xenium_alignment/module02_Spateo_second_pass_*``.
VisiumHD
    Final pair-level coordinates from
    ``processed-data/11_Xenium_VisiumHD_alignment/
    module_02_Spateo_second_pass``.

Outputs
-------
1. Four 2-by-4 pairwise overlay figures: all observations, D1_Island_A,
   D1_Island_B, and WM.
2. Six interactive HTML files: two donors x three highlighted labels,
   displayed with reversed physical z depth.
3. Six static three-view PNG files with the same donor/label combinations.

All observations remain visible in the 3D plots. Xenium highlighted labels use
bright colors; VHD uses darker counterparts. Other observations are grey.
Xenium and VHD use the same numeric point size: Xenium is rendered as flat
points, whereas VHD is rendered as spheres and with a darker grey background
so its section planes remain recognizable without implying different
observation areas.
"""

from __future__ import annotations

import argparse
import html
import os
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyvista as pv


REFERENCE_COLOR = "#4C78A8"
MOVING_COLOR = "#F58518"
BACKGROUND_XENIUM_COLOR = "#8F8F8F"
BACKGROUND_VHD_COLOR = "#585858"
HIGHLIGHT_PALETTE = {
    "D1_Island_A": "#00BFFF",
    "D1_Island_B": "#00FF0D",
    "WM": "#FFA500",
}
VHD_HIGHLIGHT_PALETTE = {
    "D1_Island_A": "#0B3C5D",
    "D1_Island_B": "#006400",
    "WM": "#8B4513",
}

PAIRWISE_ALL_POINT_SIZE = 0.20
PAIRWISE_ALL_ALPHA = 0.45
PAIRWISE_CELLTYPE_POINT_SIZE = 0.60
PAIRWISE_XENIUM_ALPHA = 0.60
PAIRWISE_VHD_ALPHA = 0.35

THREE_D_POINT_SIZE = 3.0
XENIUM_BACKGROUND_OPACITY = 0.05
VHD_BACKGROUND_OPACITY = 0.10
XENIUM_HIGHLIGHT_OPACITY = 0.55
VHD_HIGHLIGHT_OPACITY = 1.00


@dataclass(frozen=True)
class PairSpec:
    donor: str
    xenium_sample: str
    vhd_sample: str
    xenium_depth: float
    vhd_depth: float

    @property
    def tag(self) -> str:
        return f"{self.xenium_sample}-{self.vhd_sample}"


# Pairwise panel order: Br6660 first, then Br6436. Within each donor, sort by
# Xenium depth and use VHD depth as the tie-breaker.
PAIR_SPECS = (
    PairSpec(
        "Br6660",
        "Br6660_NAc4_2080",
        "VHD_H1_8MTH2TQ_A1",
        2080,
        2040,
    ),
    PairSpec(
        "Br6660",
        "Br6660_NAc6_3080",
        "VHD_H1_M3TCP9V_A1",
        3080,
        3040,
    ),
    PairSpec(
        "Br6660",
        "Br6660_NAc7_3580",
        "VHD_H1_XKYDCP3_A1",
        3580,
        3540,
    ),
    PairSpec(
        "Br6660",
        "Br6660_NAc7_3580",
        "VHD_H1_8MTH2TQ_D1",
        3580,
        3630,
    ),
    PairSpec(
        "Br6660",
        "Br6660_Nac10_4080",
        "VHD_H1_XKYDCP3_D1",
        4080,
        4120,
    ),
    PairSpec(
        "Br6660",
        "Br6660_NAc8_4580",
        "VHD_H1_M3TCP9V_D1",
        4580,
        4550,
    ),
    PairSpec(
        "Br6436",
        "Br6436_Nac1_650",
        "VHD_H1_XNQ4F2B_A1",
        650,
        630,
    ),
    PairSpec(
        "Br6436",
        "Br6436_Nac_9_4650",
        "VHD_H1_XNQ4F2B_D1",
        4650,
        4630,
    ),
)
DONORS = ("Br6660", "Br6436")


@dataclass
class DonorCloud:
    donor: str
    xenium_xyz: np.ndarray
    xenium_labels: np.ndarray
    vhd_xyz: np.ndarray
    vhd_labels: np.ndarray
    xenium_depths: tuple[float, ...]
    vhd_depths: tuple[float, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot final Xenium–VisiumHD pairwise alignments and donor-level "
            "D1 Island and WM 3D reconstructions."
        )
    )
    parser.add_argument(
        "--products",
        nargs="+",
        choices=("pairwise", "interactive", "static"),
        default=("pairwise", "interactive", "static"),
        help="Plot groups to generate (default: all three).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Override plots/11_Xenium_VisiumHD_alignment/"
            "module_03_Xenium_VHD_aligned_plots."
        ),
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of existing plot files.",
    )
    parser.add_argument(
        "--validate-only",
        action="store_true",
        help="Validate every input, pair, label, and depth without plotting.",
    )
    return parser.parse_args()


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def xenium_coordinate_path(root: Path, donor: str) -> Path:
    return (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / f"module02_Spateo_second_pass_{donor}"
        / f"{donor}_second_pass_coordinates.csv"
    )


def vhd_coordinate_path(root: Path, spec: PairSpec) -> Path:
    return (
        root
        / "processed-data"
        / "11_Xenium_VisiumHD_alignment"
        / "module_02_Spateo_second_pass"
        / spec.tag
        / f"{spec.tag}_second_pass_coordinates.csv"
    )


def require_file(path: Path) -> None:
    if not path.is_file():
        raise FileNotFoundError(path)


def require_output(path: Path, overwrite: bool) -> None:
    if path.exists() and not overwrite:
        raise FileExistsError(
            f"Output already exists; rerun with --overwrite: {path}"
        )
    path.parent.mkdir(parents=True, exist_ok=True)


def unique_numeric(
    values: pd.Series,
    description: str,
) -> float:
    observed = pd.to_numeric(values, errors="raise").dropna().unique()
    if len(observed) != 1:
        raise ValueError(
            f"Expected exactly one {description}; observed {observed.tolist()}"
        )
    return float(observed[0])


def check_depth(
    observed: float,
    expected: float,
    description: str,
) -> None:
    if not np.isclose(observed, expected):
        raise ValueError(
            f"{description} depth mismatch: expected {expected:g}, "
            f"observed {observed:g}"
        )


def load_xenium_donor(root: Path, donor: str) -> pd.DataFrame:
    path = xenium_coordinate_path(root, donor)
    require_file(path)
    required = [
        "donor",
        "sample",
        "cell_id",
        "cell_type",
        "x_second_pass",
        "y_second_pass",
        "z_height",
    ]
    data = pd.read_csv(
        path,
        usecols=required,
        dtype={
            "donor": "string",
            "sample": "string",
            "cell_id": "string",
            "cell_type": "string",
        },
    )
    if set(data["donor"].dropna().astype(str)) != {donor}:
        raise ValueError(f"Xenium coordinate donor mismatch in {path}")
    if data["cell_id"].isna().any() or data["cell_id"].duplicated().any():
        raise ValueError(f"Xenium cell IDs must be present and unique: {path}")
    coordinate_columns = ["x_second_pass", "y_second_pass", "z_height"]
    if data[coordinate_columns].isna().any().any():
        raise ValueError(f"Xenium coordinates contain missing values: {path}")
    if data["cell_type"].isna().any():
        raise ValueError(f"Xenium cell types contain missing values: {path}")
    data["cell_type"] = data["cell_type"].astype("category")
    return data


def load_vhd_pair(root: Path, spec: PairSpec) -> pd.DataFrame:
    path = vhd_coordinate_path(root, spec)
    require_file(path)
    required = [
        "donor",
        "xenium_sample",
        "vhd_sample",
        "vhd_obs_id",
        "x_second_pass",
        "y_second_pass",
        "xenium_z_height",
        "vhd_z_height",
        "Spatial_Domain",
    ]
    data = pd.read_csv(
        path,
        usecols=required,
        dtype={
            "donor": "string",
            "xenium_sample": "string",
            "vhd_sample": "string",
            "vhd_obs_id": "string",
            "Spatial_Domain": "string",
        },
    )
    expected_values = {
        "donor": spec.donor,
        "xenium_sample": spec.xenium_sample,
        "vhd_sample": spec.vhd_sample,
    }
    for column, expected in expected_values.items():
        observed = set(data[column].dropna().astype(str))
        if observed != {expected}:
            raise ValueError(
                f"{path}: expected {column}={expected!r}; "
                f"observed {sorted(observed)}"
            )
    if data["vhd_obs_id"].isna().any() or data["vhd_obs_id"].duplicated().any():
        raise ValueError(f"VHD observation IDs must be present and unique: {path}")
    coordinate_columns = ["x_second_pass", "y_second_pass", "vhd_z_height"]
    if data[coordinate_columns].isna().any().any():
        raise ValueError(f"VHD coordinates contain missing values: {path}")
    if data["Spatial_Domain"].isna().any():
        raise ValueError(f"VHD spatial domains contain missing values: {path}")
    check_depth(
        unique_numeric(data["xenium_z_height"], "Xenium depth"),
        spec.xenium_depth,
        f"{spec.tag} Xenium",
    )
    check_depth(
        unique_numeric(data["vhd_z_height"], "VHD depth"),
        spec.vhd_depth,
        f"{spec.tag} VHD",
    )
    data["Spatial_Domain"] = data["Spatial_Domain"].astype("category")
    return data


def validate_xenium_pair(
    xenium: pd.DataFrame,
    spec: PairSpec,
) -> pd.DataFrame:
    selected = xenium.loc[xenium["sample"].eq(spec.xenium_sample)]
    if selected.empty:
        available = sorted(xenium["sample"].dropna().astype(str).unique())
        raise ValueError(
            f"Missing Xenium sample {spec.xenium_sample}; "
            f"available={available}"
        )
    check_depth(
        unique_numeric(selected["z_height"], "Xenium slice depth"),
        spec.xenium_depth,
        spec.xenium_sample,
    )
    return selected


def short_sample_name(sample: str, donor: str) -> str:
    prefix = f"{donor}_"
    return sample[len(prefix) :] if sample.startswith(prefix) else sample


def plot_pairwise_grid(
    xenium_by_donor: dict[str, pd.DataFrame],
    vhd_by_pair: dict[str, pd.DataFrame],
    output_path: Path,
    overwrite: bool,
    cell_type: str | None,
) -> None:
    require_output(output_path, overwrite)
    figure, axes = plt.subplots(
        2,
        4,
        figsize=(20.0, 9.0),
        squeeze=False,
    )
    for axis, spec in zip(axes.ravel(), PAIR_SPECS):
        xenium = validate_xenium_pair(xenium_by_donor[spec.donor], spec)
        vhd = vhd_by_pair[spec.tag]
        if cell_type is not None:
            xenium = xenium.loc[xenium["cell_type"].eq(cell_type)]
            vhd = vhd.loc[vhd["Spatial_Domain"].eq(cell_type)]
            point_size = PAIRWISE_CELLTYPE_POINT_SIZE
            xenium_alpha = PAIRWISE_XENIUM_ALPHA
            vhd_alpha = PAIRWISE_VHD_ALPHA
        else:
            point_size = PAIRWISE_ALL_POINT_SIZE
            xenium_alpha = PAIRWISE_ALL_ALPHA
            vhd_alpha = PAIRWISE_ALL_ALPHA

        # Equal point sizes; opacity and color distinguish the technologies.
        axis.scatter(
            vhd["x_second_pass"],
            vhd["y_second_pass"],
            s=point_size,
            alpha=vhd_alpha,
            color=MOVING_COLOR,
            linewidths=0,
            label="VHD",
            zorder=1,
            rasterized=True,
        )
        axis.scatter(
            xenium["x_second_pass"],
            xenium["y_second_pass"],
            s=point_size,
            alpha=xenium_alpha,
            color=REFERENCE_COLOR,
            linewidths=0,
            label="Xenium",
            zorder=2,
            rasterized=True,
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (µm)")
        axis.set_ylabel("y (µm)")
        axis.set_title(
            f"{spec.donor}: "
            f"{short_sample_name(spec.xenium_sample, spec.donor)} ↔ "
            f"{spec.vhd_sample.removeprefix('VHD_')}\n"
            f"Xenium depth: {spec.xenium_depth:g} µm | "
            f"VHD depth: {spec.vhd_depth:g} µm",
            fontsize=8,
        )
        axis.legend(
            loc="upper right",
            fontsize=6,
            markerscale=7,
            frameon=False,
        )

    label = "all observations" if cell_type is None else cell_type
    figure.suptitle(
        f"Final Xenium–VisiumHD pairwise alignments: {label}",
        y=1.002,
        fontsize=13,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def build_donor_cloud(
    donor: str,
    xenium: pd.DataFrame,
    vhd_frames: Sequence[pd.DataFrame],
) -> DonorCloud:
    xenium_sorted = xenium.sort_values(
        ["z_height", "sample", "cell_id"],
        kind="stable",
    )
    vhd = pd.concat(vhd_frames, ignore_index=True)
    vhd_sorted = vhd.sort_values(
        ["vhd_z_height", "vhd_sample", "vhd_obs_id"],
        kind="stable",
    )
    xenium_xyz = xenium_sorted[
        ["x_second_pass", "y_second_pass", "z_height"]
    ].to_numpy(dtype=np.float32)
    vhd_xyz = vhd_sorted[
        ["x_second_pass", "y_second_pass", "vhd_z_height"]
    ].to_numpy(dtype=np.float32)
    xenium_labels = xenium_sorted["cell_type"].astype(str).to_numpy()
    vhd_labels = vhd_sorted["Spatial_Domain"].astype(str).to_numpy()
    xenium_depths = tuple(
        sorted(map(float, xenium_sorted["z_height"].unique()))
    )
    vhd_depths = tuple(
        sorted(map(float, vhd_sorted["vhd_z_height"].unique()))
    )
    if not xenium_depths or not vhd_depths:
        raise ValueError(f"{donor}: empty Xenium or VHD reconstruction")
    return DonorCloud(
        donor=donor,
        xenium_xyz=xenium_xyz,
        xenium_labels=xenium_labels,
        vhd_xyz=vhd_xyz,
        vhd_labels=vhd_labels,
        xenium_depths=xenium_depths,
        vhd_depths=vhd_depths,
    )


def points_for_orientation(
    points: np.ndarray,
    flipped_z: bool,
) -> np.ndarray:
    if not flipped_z:
        return points
    flipped = points.copy()
    flipped[:, 2] *= -1
    return flipped


def add_point_layer(
    plotter: pv.Plotter,
    points: np.ndarray,
    color: str,
    opacity: float,
    point_size: float,
    label: str,
    render_points_as_spheres: bool,
) -> None:
    if len(points) == 0:
        return
    model = pv.PolyData(points)
    plotter.add_mesh(
        model,
        color=color,
        opacity=opacity,
        style="points",
        point_size=point_size,
        render_points_as_spheres=render_points_as_spheres,
        ambient=0.35,
        smooth_shading=False,
        show_scalar_bar=False,
        label=label,
    )


def add_highlight_layers(
    plotter: pv.Plotter,
    cloud: DonorCloud,
    cell_type: str,
    flipped_z: bool,
    show_legend: bool,
    legend_size: tuple[float, float] = (0.30, 0.20),
    legend_font_size: int | None = None,
    legend_location: str = "upper right",
) -> None:
    xenium_highlight_color = HIGHLIGHT_PALETTE[cell_type]
    vhd_highlight_color = VHD_HIGHLIGHT_PALETTE[cell_type]
    short_cell_type = cell_type.replace("D1_Island_", "D1 ")
    xenium_points = points_for_orientation(cloud.xenium_xyz, flipped_z)
    vhd_points = points_for_orientation(cloud.vhd_xyz, flipped_z)
    xenium_highlight = cloud.xenium_labels == cell_type
    vhd_highlight = cloud.vhd_labels == cell_type

    # All layers use the same point size. VHD is shown as spherical points and
    # with darker grey so its physical slice planes remain distinguishable.
    add_point_layer(
        plotter,
        xenium_points[~xenium_highlight],
        BACKGROUND_XENIUM_COLOR,
        XENIUM_BACKGROUND_OPACITY,
        THREE_D_POINT_SIZE,
        "Xenium other",
        False,
    )
    add_point_layer(
        plotter,
        vhd_points[~vhd_highlight],
        BACKGROUND_VHD_COLOR,
        VHD_BACKGROUND_OPACITY,
        THREE_D_POINT_SIZE,
        "VHD other",
        True,
    )
    add_point_layer(
        plotter,
        xenium_points[xenium_highlight],
        xenium_highlight_color,
        XENIUM_HIGHLIGHT_OPACITY,
        THREE_D_POINT_SIZE,
        f"Xenium {short_cell_type}",
        False,
    )
    add_point_layer(
        plotter,
        vhd_points[vhd_highlight],
        vhd_highlight_color,
        VHD_HIGHLIGHT_OPACITY,
        THREE_D_POINT_SIZE,
        f"VHD {short_cell_type}",
        True,
    )
    if show_legend:
        legend = plotter.add_legend(
            loc=legend_location,
            bcolor="white",
            border=True,
            face="circle",
            size=legend_size,
        )
        if legend_font_size is not None:
            legend.GetEntryTextProperty().SetFontSize(legend_font_size)


def add_brain_reference_axes(plotter: pv.Plotter) -> None:
    """Add the same rotating anatomical axis widget as Xenium Module 02."""
    plotter.add_axes(
        interactive=False,
        line_width=3,
        color="black",
        x_color="#D95F5F",
        y_color="#5FA35F",
        z_color="#5F7FD9",
        xlabel="ML (L-M)",
        ylabel="DV (D-V)",
        zlabel="AP (A-P)",
        viewport=(0.0, 0.0, 0.22, 0.22),
        label_size=(0.34, 0.12),
    )


def inject_html_title(
    output_path: Path,
    donor: str,
    cell_type: str,
) -> None:
    overlay = f"""
<div style="position:fixed;left:18px;top:18px;z-index:1000;
            padding:8px 12px;background:rgba(255,255,255,0.90);
            border:1px solid #777;border-radius:4px;color:#111;
            font:600 18px Arial,sans-serif;">
  {html.escape(donor)} — {html.escape(cell_type)}
</div>"""
    html_text = output_path.read_text(encoding="utf-8")
    if "</body>" not in html_text:
        raise ValueError(f"Exported HTML has no closing body tag: {output_path}")
    output_path.write_text(
        html_text.replace("</body>", f"{overlay}\n</body>", 1),
        encoding="utf-8",
    )


def export_interactive(
    cloud: DonorCloud,
    cell_type: str,
    flipped_z: bool,
    output_path: Path,
    overwrite: bool,
) -> None:
    require_output(output_path, overwrite)
    plotter = pv.Plotter(off_screen=True, window_size=(1500, 1500))
    try:
        plotter.set_background("white")
        add_highlight_layers(
            plotter,
            cloud,
            cell_type,
            flipped_z,
            show_legend=True,
        )
        plotter.view_isometric()
        plotter.reset_camera()
        add_brain_reference_axes(plotter)
        plotter.export_html(str(output_path))
        inject_html_title(
            output_path,
            cloud.donor,
            cell_type,
        )
    finally:
        plotter.close()


def set_static_view(plotter: pv.Plotter, view: str) -> None:
    if view == "xy":
        plotter.view_xy()
    elif view == "xz":
        plotter.view_xz()
    elif view == "yz":
        plotter.view_yz()
    else:
        raise ValueError(f"Unknown camera view: {view}")
    plotter.reset_camera()


def export_static(
    cloud: DonorCloud,
    cell_type: str,
    flipped_z: bool,
    output_path: Path,
    overwrite: bool,
) -> None:
    require_output(output_path, overwrite)
    views = ("xy", "xz", "yz")
    plotter = pv.Plotter(
        shape=(1, 3),
        off_screen=True,
        window_size=(1800, 650),
        border=False,
    )
    try:
        for index, view in enumerate(views):
            plotter.subplot(0, index)
            plotter.set_background("white")
            add_highlight_layers(
                plotter,
                cloud,
                cell_type,
                flipped_z,
                show_legend=index == len(views) - 1,
                legend_size=(0.34, 0.15),
                legend_font_size=8,
                legend_location="upper left",
            )
            set_static_view(plotter, view)
            plotter.add_text(
                f"{cloud.donor} {cell_type} {view}",
                position="upper_edge",
                font_size=10,
                color="black",
            )
        plotter.show(
            screenshot=str(output_path),
            auto_close=False,
            interactive=False,
        )
    finally:
        plotter.close()


def initialize_headless_pyvista() -> None:
    try:
        pv.start_xvfb()
    except Exception:
        # EGL/OSMesa builds do not need an X virtual framebuffer.
        pass


def print_validation_summary(
    xenium_by_donor: dict[str, pd.DataFrame],
    vhd_by_pair: dict[str, pd.DataFrame],
) -> None:
    print("\nValidated pair inputs:")
    for spec in PAIR_SPECS:
        xenium = validate_xenium_pair(xenium_by_donor[spec.donor], spec)
        vhd = vhd_by_pair[spec.tag]
        counts = []
        for cell_type in HIGHLIGHT_PALETTE:
            counts.append(
                f"{cell_type}:X={int(xenium['cell_type'].eq(cell_type).sum()):,},"
                f"V={int(vhd['Spatial_Domain'].eq(cell_type).sum()):,}"
            )
        print(
            f"  {spec.tag}: X={len(xenium):,} at {spec.xenium_depth:g} µm; "
            f"V={len(vhd):,} at {spec.vhd_depth:g} µm; "
            + "; ".join(counts)
        )


def main() -> None:
    args = parse_args()
    root = git_root()
    output_dir = args.output_dir or (
        root
        / "plots"
        / "11_Xenium_VisiumHD_alignment"
        / "module_03_Xenium_VHD_aligned_plots"
    )

    xenium_by_donor = {
        donor: load_xenium_donor(root, donor) for donor in DONORS
    }
    vhd_by_pair = {
        spec.tag: load_vhd_pair(root, spec) for spec in PAIR_SPECS
    }
    for spec in PAIR_SPECS:
        validate_xenium_pair(xenium_by_donor[spec.donor], spec)
    print_validation_summary(xenium_by_donor, vhd_by_pair)
    if args.validate_only:
        print("\nValidation completed; no plots were generated.")
        return

    products = set(args.products)
    if "pairwise" in products:
        pairwise_dir = output_dir / "pairwise_2d"
        plot_pairwise_grid(
            xenium_by_donor,
            vhd_by_pair,
            pairwise_dir
            / "Xenium_VHD_pairwise_all_observations.png",
            args.overwrite,
            cell_type=None,
        )
        for cell_type in HIGHLIGHT_PALETTE:
            plot_pairwise_grid(
                xenium_by_donor,
                vhd_by_pair,
                pairwise_dir
                / f"Xenium_VHD_pairwise_{cell_type}.png",
                args.overwrite,
                cell_type=cell_type,
            )

    if {"interactive", "static"} & products:
        initialize_headless_pyvista()
        for donor in DONORS:
            donor_specs = [
                spec for spec in PAIR_SPECS if spec.donor == donor
            ]
            cloud = build_donor_cloud(
                donor,
                xenium_by_donor[donor],
                [vhd_by_pair[spec.tag] for spec in donor_specs],
            )
            print(
                f"\n{donor} 3D order — Xenium depths: "
                f"{list(cloud.xenium_depths)}; "
                f"VHD depths: {list(cloud.vhd_depths)}"
            )
            for cell_type in HIGHLIGHT_PALETTE:
                # Display the anatomical stack with reversed physical z depth.
                # Output names and plot annotations intentionally do not expose
                # this implementation detail.
                flipped_z = True
                if "interactive" in products:
                    export_interactive(
                        cloud,
                        cell_type,
                        flipped_z,
                        output_dir
                        / "by_celltype_3d_interactive"
                        / donor
                        / f"{cell_type}_reconstruction_3d.html",
                        args.overwrite,
                    )
                if "static" in products:
                    export_static(
                        cloud,
                        cell_type,
                        flipped_z,
                        output_dir
                        / "by_celltype_3d_static"
                        / donor
                        / f"{cell_type}_reconstruction_3d.png",
                        args.overwrite,
                    )

    print("\nSaved plots:", output_dir)
    print("Done.")


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/11_Xenium_VisiumHD_alignment
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# Validate all eight pairs and physical depths without creating plots:
# python 06_Xenium_VHD_aligned_plot.py --validate-only
#
# Generate all requested plots:
# python 06_Xenium_VHD_aligned_plot.py --overwrite
#
# Generate only the four 2-by-4 pairwise figures:
# python 06_Xenium_VHD_aligned_plot.py --products pairwise --overwrite
