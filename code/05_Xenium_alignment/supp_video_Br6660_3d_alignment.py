#!/usr/bin/env python3
"""Render reproducible Br6660 final-alignment 3D rotation videos.

The scene reproduces the all-cell-type point cloud used by
``Br6660_reconstruction_3d_second_spateo_flipped_z.html`` directly from the
final second-pass coordinate CSV. Physical z depth is multiplied by -1 to
match that HTML. Either all cells or only D1_Island_A and D1_Island_B can be
rendered. Each combined movie contains two consecutive camera orbits:

1. one full rotation about the X/medial-lateral (ML) axis;
2. one full rotation about the Y/dorsal-ventral (DV) axis.

Both rotations orbit around the point-cloud center. The first frame of each
segment is the frontal XY view, and a fixed upper-left caption identifies the
rotation axis. The lower-left anatomical axes rotate with the camera.
"""

from __future__ import annotations

import argparse
import os
import subprocess
import time
from pathlib import Path

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ.setdefault("VTK_DEFAULT_OPENGL_WINDOW", "vtkEGLRenderWindow")
# Ignore stale SSH-forwarded X displays: this movie is rendered through EGL.
os.environ.pop("DISPLAY", None)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import cv2
import numpy as np
import pandas as pd
import pyvista as pv
from matplotlib.colors import to_rgba


DONOR = "Br6660"
D1_ISLAND_CELLTYPES = ("D1_Island_A", "D1_Island_B")

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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cells",
        choices=("all", "d1-islands"),
        default="all",
        help=(
            "Cells to render: all cell types, or only D1_Island_A and "
            "D1_Island_B (default: all)."
        ),
    )
    parser.add_argument(
        "--fps",
        type=int,
        default=30,
        help="Movie frames per second (default: 30).",
    )
    parser.add_argument(
        "--seconds-per-rotation",
        type=float,
        default=8.0,
        help="Duration of each 360-degree orbit (default: 8 seconds).",
    )
    parser.add_argument(
        "--window-size",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(1440, 1440),
        help="Output frame dimensions (default: 1440 1440).",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=3.0,
        help="Rendered point size, matching the HTML default (default: 3).",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Optional output MP4 path; defaults to the supp_video folder.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow replacement of an existing movie.",
    )
    args = parser.parse_args()

    if args.fps < 1:
        parser.error("--fps must be positive")
    if not np.isfinite(args.seconds_per_rotation) or args.seconds_per_rotation <= 0:
        parser.error("--seconds-per-rotation must be positive and finite")
    if any(value < 2 or value % 2 for value in args.window_size):
        parser.error("--window-size values must be positive even integers")
    if not np.isfinite(args.point_size) or args.point_size <= 0:
        parser.error("--point-size must be positive and finite")
    return args


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def input_coordinate_path(root: Path) -> Path:
    return (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / "module02_Spateo_second_pass_Br6660"
        / "Br6660_second_pass_coordinates.csv"
    )


def default_output_path(root: Path, cells: str) -> Path:
    filename = (
        "Br6660_reconstruction_3d_second_spateo_all_cells.mp4"
        if cells == "all"
        else "Br6660_reconstruction_3d_second_spateo_D1_Island_A_B.mp4"
    )
    return (
        root
        / "plots"
        / "05_Xenium_alignment"
        / "supp_video"
        / filename
    )


def prepare_output(path: Path, overwrite: bool) -> None:
    if path.suffix.lower() != ".mp4":
        raise ValueError(f"Output must have an .mp4 extension: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {path}")
    if path.exists():
        path.unlink()


def load_flipped_point_cloud(
    root: Path,
    cells: str,
) -> tuple[pv.PolyData, pd.DataFrame]:
    """Load requested final cells and reproduce the HTML's z *= -1 operation."""
    path = input_coordinate_path(root)
    required = {
        "donor",
        "sample",
        "cell_id",
        "cell_type",
        "x_second_pass",
        "y_second_pass",
        "z_height",
    }
    frame = pd.read_csv(path)
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"Coordinate CSV is missing columns: {sorted(missing)}")
    if set(frame["donor"].dropna().astype(str)) != {DONOR}:
        raise ValueError("Coordinate CSV contains an unexpected donor")
    if frame["cell_id"].astype(str).duplicated().any():
        raise ValueError("Coordinate CSV contains duplicate cell IDs")

    if cells == "d1-islands":
        frame = frame.loc[
            frame["cell_type"].astype(str).isin(D1_ISLAND_CELLTYPES)
        ].copy()
        observed = set(frame["cell_type"].astype(str))
        missing_d1 = set(D1_ISLAND_CELLTYPES) - observed
        if missing_d1:
            raise ValueError(
                f"Requested D1 Island cell types are absent: {sorted(missing_d1)}"
            )
    elif cells != "all":
        raise ValueError(f"Unsupported cell selection: {cells}")
    if frame.empty:
        raise ValueError(f"No cells remain for selection {cells!r}")

    coordinate_columns = ["x_second_pass", "y_second_pass", "z_height"]
    coordinates = frame[coordinate_columns].to_numpy(dtype=np.float64)
    if not np.isfinite(coordinates).all():
        raise ValueError("Coordinate CSV contains non-finite spatial values")
    coordinates[:, 2] *= -1.0

    labels = frame["cell_type"].astype(str)
    unknown = sorted(set(labels) - set(CELLTYPE_PALETTE))
    if unknown:
        raise ValueError(f"Cell types missing from palette: {unknown}")

    palette_labels = list(CELLTYPE_PALETTE)
    color_table = np.vstack(
        [
            np.rint(np.asarray(to_rgba(CELLTYPE_PALETTE[label])) * 255)
            .astype(np.uint8)
            for label in palette_labels
        ]
    )
    label_codes = pd.Categorical(labels, categories=palette_labels).codes
    if np.any(label_codes < 0):
        raise ValueError("Failed to encode one or more cell-type colors")
    rgba = color_table[label_codes]

    cloud = pv.PolyData(coordinates)
    cloud["annotation_rgba"] = rgba
    return cloud, frame


def add_brain_reference_axes(plotter: pv.Plotter) -> None:
    """Add the same rotating anatomical axes used by the flipped-Z HTML."""
    plotter.add_axes(
        interactive=False,
        line_width=3,
        color="black",
        x_color="#D95F5F",
        y_color="#5FA35F",
        z_color="#5F7FD9",
        xlabel="ML",
        ylabel="DV",
        zlabel="AP",
        viewport=(0.0, 0.0, 0.22, 0.22),
        label_size=(0.34, 0.12),
    )


def add_d1_island_legend(plotter: pv.Plotter) -> None:
    """Add a compact, inset legend used only by the D1-Island video."""
    legend = plotter.add_legend(
        labels=[
            ["D1 Island A", CELLTYPE_PALETTE["D1_Island_A"], "circle"],
            ["D1 Island B", CELLTYPE_PALETTE["D1_Island_B"], "circle"],
        ],
        bcolor="white",
        border=False,
        size=(0.15, 0.055),
        loc=None,
        background_opacity=0.78,
        font_family="arial",
        name="d1_island_legend",
    )
    # Normalized viewport coordinates: upper right with a small inset.
    legend.SetPosition(0.79, 0.89)
    legend.SetPosition2(0.15, 0.055)
    text_property = legend.GetEntryTextProperty()
    text_property.SetFontSize(8)
    text_property.SetBold(False)


def rotation_x(angle: float) -> np.ndarray:
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array(
        [[1.0, 0.0, 0.0], [0.0, cosine, -sine], [0.0, sine, cosine]]
    )


def rotation_y(angle: float) -> np.ndarray:
    cosine, sine = np.cos(angle), np.sin(angle)
    return np.array(
        [[cosine, 0.0, sine], [0.0, 1.0, 0.0], [-sine, 0.0, cosine]]
    )


def set_camera(
    plotter: pv.Plotter,
    center: np.ndarray,
    initial_offset: np.ndarray,
    initial_up: np.ndarray,
    rotation: np.ndarray,
) -> None:
    """Orbit the complete camera rig around the fixed point-cloud center."""
    offset = rotation @ initial_offset
    view_up = rotation @ initial_up
    plotter.camera.position = tuple(center + offset)
    plotter.camera.focal_point = tuple(center)
    plotter.camera.up = tuple(view_up)
    plotter.reset_camera_clipping_range()
    plotter.render()


def render_orbit(
    plotter: pv.Plotter,
    writer: cv2.VideoWriter,
    caption: str,
    center: np.ndarray,
    initial_offset: np.ndarray,
    initial_up: np.ndarray,
    rotation_function,
    frame_count: int,
    frame_offset: int,
    total_frames: int,
    started: float,
) -> None:
    """Write one complete 360-degree orbit to the open movie writer."""
    plotter.add_text(
        caption,
        position=(0.025, 0.90),
        viewport=True,
        font_size=14,
        color="black",
        shadow=False,
        name="rotation_caption",
    )
    angles = np.linspace(0.0, 2.0 * np.pi, frame_count, endpoint=True)
    for local_index, angle in enumerate(angles):
        set_camera(
            plotter,
            center,
            initial_offset,
            initial_up,
            rotation_function(float(angle)),
        )
        rgb_frame = np.asarray(plotter.image)
        expected_shape = (
            int(plotter.window_size[1]),
            int(plotter.window_size[0]),
        )
        if rgb_frame.shape[:2] != expected_shape or rgb_frame.shape[2] != 3:
            raise ValueError(
                f"Unexpected rendered frame shape {rgb_frame.shape}; "
                f"expected {expected_shape + (3,)}"
            )
        writer.write(cv2.cvtColor(rgb_frame, cv2.COLOR_RGB2BGR))

        completed = frame_offset + local_index + 1
        if completed == 1 or completed % 10 == 0 or completed == total_frames:
            elapsed = time.perf_counter() - started
            render_fps = completed / elapsed
            eta_seconds = (total_frames - completed) / render_fps
            print(
                f"Frames {completed}/{total_frames} | "
                f"render rate {render_fps:.2f} frame/s | "
                f"ETA {eta_seconds / 60:.1f} min",
                flush=True,
            )


def report_renderer(plotter: pv.Plotter) -> None:
    """Report the actual VTK/OpenGL backend rather than assuming GPU use."""
    render_window = plotter.render_window
    window_class = render_window.GetClassName()
    capabilities = str(render_window.ReportCapabilities())
    print(f"VTK render window: {window_class}", flush=True)
    informative_lines = [
        line.strip()
        for line in capabilities.splitlines()
        if any(
            token in line.lower()
            for token in ("vendor string", "renderer string", "version string")
        )
    ]
    for line in informative_lines:
        print(line, flush=True)
    lowered = capabilities.lower()
    if window_class != "vtkEGLRenderWindow":
        raise RuntimeError(
            "EGL was requested but VTK created "
            f"{window_class}. Do not continue with this 480-frame render."
        )
    if "llvmpipe" in lowered or "softpipe" in lowered:
        raise RuntimeError(
            "VTK selected a software OpenGL renderer instead of the GPU."
        )


def render_video(
    cloud: pv.PolyData,
    output_path: Path,
    *,
    cells: str,
    fps: int,
    seconds_per_rotation: float,
    window_size: tuple[int, int],
    point_size: float,
) -> None:
    """Render the ML orbit followed immediately by the DV orbit."""
    frame_count = max(2, int(round(fps * seconds_per_rotation)))
    temporary_path = output_path.with_name(
        f"{output_path.stem}.tmp{output_path.suffix}"
    )
    if temporary_path.exists():
        temporary_path.unlink()
    writer = None
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    try:
        plotter.set_background("white")
        plotter.add_mesh(
            cloud,
            scalars="annotation_rgba",
            rgba=True,
            style="points",
            render_points_as_spheres=True,
            point_size=point_size,
            ambient=0.2,
            smooth_shading=True,
            show_scalar_bar=False,
        )
        add_brain_reference_axes(plotter)
        if cells == "d1-islands":
            add_d1_island_legend(plotter)

        # Standard positive-Z frontal view of the displayed (already flipped-Z)
        # point cloud. This is the positive XY plane requested for frame one.
        plotter.view_xy()
        plotter.reset_camera()
        # Initialize the off-screen framebuffer before accessing plotter.image.
        # PyVista's open_movie() used to do this implicitly; our explicit
        # OpenCV/FFmpeg writer requires this one-time initialization here.
        plotter.show(auto_close=False, interactive=False)
        report_renderer(plotter)
        center = np.asarray(plotter.camera.focal_point, dtype=float)
        initial_position = np.asarray(plotter.camera.position, dtype=float)
        initial_offset = initial_position - center
        initial_up = np.asarray(plotter.camera.up, dtype=float)

        fourcc = cv2.VideoWriter_fourcc(*"mp4v")
        writer = cv2.VideoWriter(
            str(temporary_path),
            fourcc,
            float(fps),
            tuple(window_size),
        )
        if not writer.isOpened():
            raise RuntimeError(
                "OpenCV/FFmpeg could not open the MP4 writer for "
                f"{temporary_path}"
            )
        total_frames = 2 * frame_count
        started = time.perf_counter()
        render_orbit(
            plotter,
            writer,
            "Rotation about the ML axis",
            center,
            initial_offset,
            initial_up,
            rotation_x,
            frame_count,
            0,
            total_frames,
            started,
        )
        render_orbit(
            plotter,
            writer,
            "Rotation about the DV axis",
            center,
            initial_offset,
            initial_up,
            rotation_y,
            frame_count,
            frame_count,
            total_frames,
            started,
        )
        writer.release()
        writer = None
        temporary_path.replace(output_path)
    finally:
        if writer is not None:
            writer.release()
        plotter.close()


def main() -> None:
    args = parse_args()
    root = git_root()
    output_path = args.output or default_output_path(root, args.cells)
    if not output_path.is_absolute():
        output_path = Path.cwd() / output_path
    prepare_output(output_path, args.overwrite)

    cloud, metadata = load_flipped_point_cloud(root, args.cells)
    frame_count = max(2, int(round(args.fps * args.seconds_per_rotation)))
    print(f"Cell selection: {args.cells}")
    print(f"Loaded {cloud.n_points:,} cells across {metadata['sample'].nunique()} slices")
    print(f"Cell types: {metadata['cell_type'].nunique()}")
    print(
        f"Rendering {2 * frame_count:,} frames at {args.fps} fps "
        f"({2 * frame_count / args.fps:.1f} seconds total)"
    )
    render_video(
        cloud,
        output_path,
        cells=args.cells,
        fps=args.fps,
        seconds_per_rotation=args.seconds_per_rotation,
        window_size=tuple(args.window_size),
        point_size=args.point_size,
    )
    print(f"Saved: {output_path}")


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/05_Xenium_alignment
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# Best/smoothest default version:
#   1440 x 1440, 30 fps, 8 seconds per rotation, 16 seconds total.
#   This is the preferred full-quality version; all defaults remain unchanged.
# python supp_video_Br6660_3d_alignment.py --cells all --overwrite
#
# D1 Island A and B only, using the same full-quality defaults:
# python supp_video_Br6660_3d_alignment.py --cells d1-islands --overwrite
#
# Smaller, broadly compatible version (same resolution and rotation speed):
#   1440 x 1440, 15 fps, 8 seconds per rotation, 16 seconds total.
#   Expected size is approximately 115--125 MB with the current mp4v encoder.
# python supp_video_Br6660_3d_alignment.py \
#     --cells all \
#     --fps 15 \
#     --seconds-per-rotation 8 \
#     --window-size 1440 1440 \
#     --overwrite
#
# Lower-resolution preview at the same 8-second-per-rotation playback speed:
# python supp_video_Br6660_3d_alignment.py \
#     --cells all \
#     --fps 10 \
#     --seconds-per-rotation 8 \
#     --window-size 720 720 \
#     --output ../../plots/05_Xenium_alignment/supp_video/preview.mp4 \
#     --overwrite
