#!/usr/bin/env python3
"""Render supplementary rotation videos from final 3D heatmap HTML files.

The Plotly HTML files contain the final regular-grid coordinates, fitted GAM
proportions, display thresholds, opacity, surface count, and Viridis colorscale.
This script reads those saved values directly and reconstructs the heatmaps as
PyVista isosurfaces; it does not rerun voxelization or GAM fitting.

Each movie starts in the frontal XY view, rotates 90 degrees about the ML (X)
axis, and then continues from that camera pose through a 360-degree rotation
about the AP (Z) axis. The proportion colorbar and Plotly background grid are
intentionally omitted. A lower-left ML/DV/AP reference-axis widget rotates
with the model, and an upper-left caption identifies both the cell type and
the current orbit.
"""

from __future__ import annotations

import argparse
import base64
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ.setdefault("VTK_DEFAULT_OPENGL_WINDOW", "vtkEGLRenderWindow")
os.environ.pop("DISPLAY", None)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

import cv2
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pyvista as pv


@dataclass(frozen=True)
class VideoSpec:
    key: str
    donor: str
    cell_type: str

    @property
    def html_name(self) -> str:
        return f"{self.donor}_3d_heatmap_figure_{self.cell_type}.html"

    @property
    def video_name(self) -> str:
        return f"{self.donor}_3d_heatmap_{self.cell_type}.mp4"


VIDEO_SPECS = (
    VideoSpec("Br6660_D1_Island_A", "Br6660", "D1_Island_A"),
    VideoSpec("Br6660_D1_Island_B", "Br6660", "D1_Island_B"),
    VideoSpec("Br6660_Excitatory", "Br6660", "Excitatory"),
    VideoSpec("Br6660_WM", "Br6660", "WM"),
)
VIDEO_SPEC_BY_KEY = {spec.key: spec for spec in VIDEO_SPECS}


@dataclass
class HeatmapVolume:
    grid: pv.StructuredGrid
    isomin: float
    isomax: float
    opacity: float
    surface_count: int
    colormap: LinearSegmentedColormap


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--videos",
        nargs="+",
        choices=tuple(VIDEO_SPEC_BY_KEY),
        default=list(VIDEO_SPEC_BY_KEY),
        help="Heatmap videos to render (default: all four).",
    )
    parser.add_argument("--fps", type=int, default=30)
    parser.add_argument(
        "--seconds-per-90-degrees",
        type=float,
        default=4.0,
        help=(
            "Rotation speed expressed as seconds per 90 degrees (default: 4). "
            "Thus ML lasts 4 seconds and AP lasts 16 seconds."
        ),
    )
    parser.add_argument(
        "--window-size",
        type=int,
        nargs=2,
        metavar=("WIDTH", "HEIGHT"),
        default=(1440, 1440),
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.fps < 1:
        parser.error("--fps must be positive")
    if (
        not np.isfinite(args.seconds_per_90_degrees)
        or args.seconds_per_90_degrees <= 0
    ):
        parser.error("--seconds-per-90-degrees must be positive and finite")
    if any(value < 2 or value % 2 for value in args.window_size):
        parser.error("--window-size values must be positive even integers")
    if len(args.videos) != len(set(args.videos)):
        parser.error("--videos contains duplicates")
    return args


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def html_path(root: Path, spec: VideoSpec) -> Path:
    return (
        root
        / "plots"
        / "09_spatial_gradient"
        / "cell_type_pseudotime"
        / "3d_heatmap_figure"
        / spec.html_name
    )


def output_path(root: Path, spec: VideoSpec) -> Path:
    return (
        root
        / "plots"
        / "09_spatial_gradient"
        / "cell_type_pseudotime"
        / "supp_video"
        / spec.video_name
    )


def prepare_output(path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {path}")
    if path.exists():
        path.unlink()


def decode_plotly_array(encoded: object, name: str) -> np.ndarray:
    """Decode Plotly 3 typed-array JSON or accept a plain numeric list."""
    if isinstance(encoded, list):
        result = np.asarray(encoded, dtype=float)
    elif isinstance(encoded, dict) and {"dtype", "bdata"} <= set(encoded):
        dtype = np.dtype(str(encoded["dtype"]))
        result = np.frombuffer(base64.b64decode(encoded["bdata"]), dtype=dtype)
    else:
        raise TypeError(f"Unsupported Plotly representation for {name}")
    if result.ndim != 1 or not len(result):
        raise ValueError(f"Plotly array {name} must be a nonempty vector")
    return result.astype(np.float64, copy=False)


def parse_plotly_volume(path: Path) -> HeatmapVolume:
    """Extract the single go.Volume trace and create a structured grid."""
    text = path.read_text(encoding="utf-8")
    marker = "Plotly.newPlot("
    position = text.find(marker)
    if position < 0:
        raise ValueError(f"No Plotly.newPlot call found in {path}")
    position += len(marker)
    decoder = json.JSONDecoder()

    def decode_next(start: int) -> tuple[object, int]:
        while start < len(text) and (text[start].isspace() or text[start] == ","):
            start += 1
        return decoder.raw_decode(text, start)

    _, position = decode_next(position)  # div identifier
    data, position = decode_next(position)
    _layout, _ = decode_next(position)
    if not isinstance(data, list) or len(data) != 1:
        raise ValueError(f"Expected exactly one Plotly trace in {path}")
    trace = data[0]
    if not isinstance(trace, dict) or trace.get("type") != "volume":
        raise ValueError(f"Expected one Plotly volume trace in {path}")

    x = decode_plotly_array(trace["x"], "x")
    y = decode_plotly_array(trace["y"], "y")
    z = decode_plotly_array(trace["z"], "z")
    values = decode_plotly_array(trace["value"], "value")
    if not (len(x) == len(y) == len(z) == len(values)):
        raise ValueError(f"Plotly volume arrays have different lengths in {path}")

    nx, ny, nz = len(np.unique(x)), len(np.unique(y)), len(np.unique(z))
    if nx * ny * nz != len(values):
        raise ValueError(
            f"Volume is not a complete regular grid in {path}: "
            f"{nx}*{ny}*{nz} != {len(values)}"
        )
    shape = (nx, ny, nz)
    x_grid = x.reshape(shape, order="C")
    y_grid = y.reshape(shape, order="C")
    z_grid = z.reshape(shape, order="C")
    value_grid = values.reshape(shape, order="C")
    grid = pv.StructuredGrid(x_grid, y_grid, z_grid)
    grid["fitted_proportion"] = value_grid.ravel(order="F")

    isomin = float(trace["isomin"])
    isomax = float(trace["isomax"])
    opacity = float(trace["opacity"])
    surface_count = int(trace["surface"]["count"])
    if not np.isfinite([isomin, isomax, opacity]).all() or isomin >= isomax:
        raise ValueError(f"Invalid volume display parameters in {path}")
    colorscale = trace["colorscale"]
    if not isinstance(colorscale, list) or len(colorscale) < 2:
        raise ValueError(f"Invalid colorscale in {path}")
    colormap = LinearSegmentedColormap.from_list(
        "saved_plotly_colorscale",
        [(float(position), str(color)) for position, color in colorscale],
    )
    print(
        f"Loaded {path.name}: grid={shape}, isomin={isomin:.5g}, "
        f"isomax={isomax:.5g}, opacity={opacity:g}, surfaces={surface_count}",
        flush=True,
    )
    return HeatmapVolume(
        grid=grid,
        isomin=isomin,
        isomax=isomax,
        opacity=opacity,
        surface_count=surface_count,
        colormap=colormap,
    )


def add_heatmap(plotter: pv.Plotter, volume: HeatmapVolume) -> None:
    levels = np.linspace(volume.isomin, volume.isomax, volume.surface_count)
    surfaces = volume.grid.contour(
        isosurfaces=levels,
        scalars="fitted_proportion",
        preference="point",
    )
    if surfaces.n_points == 0:
        raise ValueError("No isosurfaces were generated from the saved volume")
    plotter.add_mesh(
        surfaces,
        scalars="fitted_proportion",
        cmap=volume.colormap,
        clim=(volume.isomin, volume.isomax),
        opacity=volume.opacity,
        lighting=False,
        show_scalar_bar=False,
    )

    # Add denser surfaces only across the upper 50% of the displayed range.
    # Keeping this overlay away from the low-value purple region avoids a hard
    # green/purple boundary while retaining the emphasized yellow-green tail.
    high_value_start = volume.isomin + 0.50 * (volume.isomax - volume.isomin)
    high_value_levels = np.linspace(high_value_start, volume.isomax, 50)
    high_value_colormap = LinearSegmentedColormap.from_list(
        "high_value_viridis",
        volume.colormap(np.linspace(0.55, 1.0, 256)),
    )
    high_value_surfaces = volume.grid.contour(
        isosurfaces=high_value_levels,
        scalars="fitted_proportion",
        preference="point",
    )
    if high_value_surfaces.n_points:
        plotter.add_mesh(
            high_value_surfaces,
            scalars="fitted_proportion",
            cmap=high_value_colormap,
            clim=(high_value_start, volume.isomax),
            opacity=min(1.0, volume.opacity * 1.15),
            lighting=False,
            show_scalar_bar=False,
        )

    # Keep the previously accepted bright-yellow treatment for the extreme
    # high-value tail.
    peak_start = volume.isomin + 0.88 * (volume.isomax - volume.isomin)
    peak_levels = np.linspace(peak_start, volume.isomax, 12)
    peak_colormap = LinearSegmentedColormap.from_list(
        "peak_viridis",
        volume.colormap(np.linspace(0.90, 1.0, 256)),
    )
    peak_surfaces = volume.grid.contour(
        isosurfaces=peak_levels,
        scalars="fitted_proportion",
        preference="point",
    )
    if peak_surfaces.n_points:
        plotter.add_mesh(
            peak_surfaces,
            scalars="fitted_proportion",
            cmap=peak_colormap,
            clim=(peak_start, volume.isomax),
            opacity=0.45,
            lighting=False,
            show_scalar_bar=False,
        )

def add_brain_reference_axes(plotter: pv.Plotter) -> None:
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


def add_spatial_bounding_box(
    plotter: pv.Plotter,
    bounds: tuple[float, float, float, float, float, float],
) -> None:
    """Add a minimal shared spatial frame and anatomical edge directions."""
    xmin, xmax, ymin, ymax, zmin, zmax = map(float, bounds)
    box = pv.Box(bounds=bounds)
    plotter.add_mesh(
        box,
        style="wireframe",
        color="#B8B8B8",
        opacity=0.55,
        line_width=1.2,
        lighting=False,
        show_scalar_bar=False,
    )

    # Three mutually non-overlapping box edges encode the coordinate polarity.
    directional_edges = (
        ((xmin, ymin, zmin), (xmax, ymin, zmin)),  # X: L -> M
        ((xmin, ymin, zmax), (xmin, ymax, zmax)),  # Y: V -> D
        ((xmax, ymax, zmin), (xmax, ymax, zmax)),  # Z: P -> A
    )
    for start, end in directional_edges:
        plotter.add_mesh(
            pv.Line(start, end),
            color="#707070",
            line_width=2.2,
            lighting=False,
            show_scalar_bar=False,
        )

    dx, dy, dz = xmax - xmin, ymax - ymin, zmax - zmin
    label_points = np.array(
        [
            (xmin, ymin - 0.025 * dy, zmin - 0.025 * dz),  # L
            (xmax, ymin - 0.025 * dy, zmin - 0.025 * dz),  # M
            (xmin - 0.025 * dx, ymax, zmax + 0.025 * dz),  # D
            (xmin - 0.025 * dx, ymin, zmax + 0.025 * dz),  # V
            (xmax + 0.025 * dx, ymax + 0.025 * dy, zmax),  # A
            (xmax + 0.025 * dx, ymax + 0.025 * dy, zmin),  # P
        ],
        dtype=float,
    )
    plotter.add_point_labels(
        label_points,
        ["L", "M", "D", "V", "A", "P"],
        show_points=False,
        shape=None,
        font_size=18,
        bold=True,
        text_color="#505050",
        always_visible=False,
        name="anatomical_box_labels",
    )


def rotation_x(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]], dtype=float)


def rotation_z(angle: float) -> np.ndarray:
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]], dtype=float)


def set_camera(
    plotter: pv.Plotter,
    center: np.ndarray,
    initial_offset: np.ndarray,
    initial_up: np.ndarray,
    rotation: np.ndarray,
) -> None:
    plotter.camera.position = tuple(center + rotation @ initial_offset)
    plotter.camera.focal_point = tuple(center)
    plotter.camera.up = tuple(rotation @ initial_up)
    plotter.reset_camera_clipping_range()
    plotter.render()


def report_renderer(plotter: pv.Plotter) -> None:
    window = plotter.render_window
    window_class = window.GetClassName()
    capabilities = str(window.ReportCapabilities())
    print(f"VTK render window: {window_class}", flush=True)
    for line in capabilities.splitlines():
        if any(
            token in line.lower()
            for token in ("vendor string", "renderer string", "version string")
        ):
            print(line.strip(), flush=True)
    lowered = capabilities.lower()
    if window_class != "vtkEGLRenderWindow":
        raise RuntimeError(f"Expected vtkEGLRenderWindow; got {window_class}")
    if "llvmpipe" in lowered or "softpipe" in lowered:
        raise RuntimeError("VTK selected software rendering instead of the GPU")


def write_orbit(
    plotter: pv.Plotter,
    writer: cv2.VideoWriter,
    caption: str,
    center: np.ndarray,
    initial_offset: np.ndarray,
    initial_up: np.ndarray,
    rotation_function: Callable[[float], np.ndarray],
    starting_rotation: np.ndarray,
    sweep_radians: float,
    frame_count: int,
    frame_offset: int,
    total_frames: int,
    started: float,
    include_start: bool,
) -> np.ndarray:
    plotter.add_text(
        caption,
        position=(0.025, 0.90),
        viewport=True,
        font_size=14,
        color="black",
        name="rotation_caption",
    )
    if include_start:
        angles = np.linspace(0.0, sweep_radians, frame_count, endpoint=True)
    else:
        angles = np.linspace(
            sweep_radians / frame_count,
            sweep_radians,
            frame_count,
            endpoint=True,
        )
    for local_index, angle in enumerate(angles):
        rotation = rotation_function(float(angle)) @ starting_rotation
        set_camera(
            plotter, center, initial_offset, initial_up,
            rotation,
        )
        rgb = np.asarray(plotter.image)
        expected = (int(plotter.window_size[1]), int(plotter.window_size[0]), 3)
        if rgb.shape != expected:
            raise ValueError(f"Unexpected frame shape {rgb.shape}; expected {expected}")
        writer.write(cv2.cvtColor(rgb, cv2.COLOR_RGB2BGR))
        completed = frame_offset + local_index + 1
        if completed == 1 or completed % 10 == 0 or completed == total_frames:
            elapsed = time.perf_counter() - started
            rate = completed / elapsed
            eta = (total_frames - completed) / rate / 60
            print(
                f"Frames {completed}/{total_frames} | "
                f"{rate:.2f} frame/s | ETA {eta:.1f} min",
                flush=True,
            )
    return rotation_function(sweep_radians) @ starting_rotation


def render_video(
    volume: HeatmapVolume,
    path: Path,
    *,
    cell_type: str,
    fps: int,
    seconds_per_90_degrees: float,
    window_size: tuple[int, int],
    overwrite: bool,
) -> None:
    prepare_output(path, overwrite)
    temporary = path.with_name(f"{path.stem}.tmp{path.suffix}")
    if temporary.exists():
        temporary.unlink()
    ml_frame_count = max(2, int(round(fps * seconds_per_90_degrees)))
    ap_frame_count = max(2, int(round(fps * 4 * seconds_per_90_degrees)))
    total_frames = ml_frame_count + ap_frame_count
    writer = None
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    try:
        plotter.set_background("white")
        add_heatmap(plotter, volume)
        add_spatial_bounding_box(plotter, volume.grid.bounds)
        add_brain_reference_axes(plotter)
        plotter.view_xy()
        plotter.reset_camera()
        plotter.show(auto_close=False, interactive=False)
        report_renderer(plotter)
        center = np.asarray(plotter.camera.focal_point, dtype=float)
        initial_offset = np.asarray(plotter.camera.position, dtype=float) - center
        initial_up = np.asarray(plotter.camera.up, dtype=float)

        writer = cv2.VideoWriter(
            str(temporary), cv2.VideoWriter_fourcc(*"mp4v"),
            float(fps), tuple(window_size),
        )
        if not writer.isOpened():
            raise RuntimeError(f"Could not open MP4 writer: {temporary}")
        started = time.perf_counter()
        current_rotation = np.eye(3)
        cell_type_label = cell_type.replace("_", " ")
        current_rotation = write_orbit(
            plotter, writer,
            f"{cell_type_label}: Rotation about the ML axis",
            center, initial_offset, initial_up, rotation_x,
            current_rotation, np.pi / 2, ml_frame_count, 0,
            total_frames, started, True,
        )
        write_orbit(
            plotter, writer,
            f"{cell_type_label}: Rotation about the AP axis",
            center, initial_offset, initial_up, rotation_z,
            current_rotation, 2 * np.pi, ap_frame_count, ml_frame_count,
            total_frames, started, False,
        )
        writer.release()
        writer = None
        temporary.replace(path)
    finally:
        if writer is not None:
            writer.release()
        plotter.close()


def main() -> None:
    args = parse_args()
    root = git_root()
    for key in args.videos:
        spec = VIDEO_SPEC_BY_KEY[key]
        source = html_path(root, spec)
        destination = output_path(root, spec)
        print(f"\nPreparing {spec.key} from {source.name}", flush=True)
        volume = parse_plotly_volume(source)
        render_video(
            volume,
            destination,
            cell_type=spec.cell_type,
            fps=args.fps,
            seconds_per_90_degrees=args.seconds_per_90_degrees,
            window_size=tuple(args.window_size),
            overwrite=args.overwrite,
        )
        print(f"Saved: {destination}", flush=True)


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/09_spatial_gradient/cell_type_pseudotime
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# Quick first-look preview using the Br6660 WM heatmap:
# python supp_video_3d_heatmap.py \
#     --videos Br6660_WM \
#     --fps 10 \
#     --seconds-per-90-degrees 4 \
#     --window-size 720 720 \
#     --overwrite
#
# All four full-quality videos (30 fps, 1440 x 1440, 20 seconds each):
# python supp_video_3d_heatmap.py --overwrite
#
# All four smaller-file videos at the same rotation speed and resolution:
# python supp_video_3d_heatmap.py \
#     --fps 15 \
#     --seconds-per-90-degrees 4 \
#     --window-size 1440 1440 \
#     --overwrite
