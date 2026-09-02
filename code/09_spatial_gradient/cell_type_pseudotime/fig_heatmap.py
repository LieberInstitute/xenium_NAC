#!/usr/bin/env python3
"""Render frontal XY screenshots of the four final Br6660 3D heatmaps.

The HTML files are used only as containers for the saved fitted-grid data.
Each screenshot independently reproduces the current supplementary video's
PyVista isosurface rendering, including its color, surface-density, and opacity
settings; Plotly/HTML rendering is not used. Each PNG has a white background,
an unlabeled fitted-proportion colorbar, and a minimal 3D bounding box. Titles,
reference-axis widgets, background grids, and direction letters are omitted.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ.setdefault("VTK_DEFAULT_OPENGL_WINDOW", "vtkEGLRenderWindow")
os.environ.pop("DISPLAY", None)
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pyvista as pv

from supp_video_3d_heatmap import (
    VIDEO_SPECS,
    HeatmapVolume,
    VideoSpec,
    git_root,
    html_path,
    parse_plotly_volume,
    report_renderer,
)


SPEC_BY_CELLTYPE = {spec.cell_type: spec for spec in VIDEO_SPECS}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cell-types",
        nargs="+",
        choices=tuple(SPEC_BY_CELLTYPE),
        default=list(SPEC_BY_CELLTYPE),
        help="Cell types to render (default: all four).",
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
    if any(value < 2 for value in args.window_size):
        parser.error("--window-size values must be at least 2")
    if len(args.cell_types) != len(set(args.cell_types)):
        parser.error("--cell-types contains duplicates")
    return args


def output_path(root: Path, spec: VideoSpec) -> Path:
    return (
        root
        / "plots"
        / "09_spatial_gradient"
        / "cell_type_pseudotime"
        / "3d_heatmap_figure"
        / f"{spec.donor}_3d_heatmap_front_xy_{spec.cell_type}.png"
    )


def prepare_output(path: Path, overwrite: bool) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists() and not overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {path}")


def add_clean_bounding_box(
    plotter: pv.Plotter,
    bounds: tuple[float, float, float, float, float, float],
) -> None:
    """Draw the shared box and three thicker edges without direction labels."""
    xmin, xmax, ymin, ymax, zmin, zmax = map(float, bounds)
    plotter.add_mesh(
        pv.Box(bounds=bounds),
        style="wireframe",
        color="#B8B8B8",
        opacity=0.55,
        line_width=6.0,
        lighting=False,
        show_scalar_bar=False,
    )

    emphasized_edges = (
        ((xmin, ymin, zmin), (xmax, ymin, zmin)),
        ((xmin, ymin, zmax), (xmin, ymax, zmax)),
        ((xmax, ymax, zmin), (xmax, ymax, zmax)),
    )
    for start, end in emphasized_edges:
        plotter.add_mesh(
            pv.Line(start, end),
            color="#606060",
            line_width=10.0,
            lighting=False,
            show_scalar_bar=False,
        )


def add_video_style_heatmap(
    plotter: pv.Plotter,
    volume: HeatmapVolume,
) -> pv.Actor:
    """Render the heatmap using the exact current supplementary-video style."""
    levels = np.linspace(volume.isomin, volume.isomax, volume.surface_count)
    surfaces = volume.grid.contour(
        isosurfaces=levels,
        scalars="fitted_proportion",
        preference="point",
    )
    if surfaces.n_points == 0:
        raise ValueError("No isosurfaces were generated from the saved volume")
    base_actor = plotter.add_mesh(
        surfaces,
        scalars="fitted_proportion",
        cmap=volume.colormap,
        clim=(volume.isomin, volume.isomax),
        opacity=volume.opacity,
        lighting=False,
        show_scalar_bar=False,
    )

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

    return base_actor


def render_front_xy(
    source: Path,
    destination: Path,
    *,
    window_size: tuple[int, int],
    overwrite: bool,
) -> None:
    prepare_output(destination, overwrite)
    volume = parse_plotly_volume(source)
    plotter = pv.Plotter(off_screen=True, window_size=window_size)
    try:
        plotter.set_background("white")
        base_actor = add_video_style_heatmap(plotter, volume)
        add_clean_bounding_box(plotter, volume.grid.bounds)
        plotter.add_scalar_bar(
            title="",
            mapper=base_actor.mapper,
            vertical=True,
            position_x=0.92,
            position_y=0.18,
            width=0.025,
            height=0.64,
            n_labels=5,
            label_font_size=40,
            color="black",
            fmt="%.2f",
            outline=False,
        )
        plotter.view_xy()
        plotter.reset_camera()
        plotter.show(auto_close=False, interactive=False)
        report_renderer(plotter)
        plotter.screenshot(str(destination), return_img=False)
    finally:
        plotter.close()


def main() -> None:
    args = parse_args()
    root = git_root()
    for cell_type in args.cell_types:
        spec = SPEC_BY_CELLTYPE[cell_type]
        source = html_path(root, spec)
        destination = output_path(root, spec)
        print(f"Rendering {cell_type} from {source.name}", flush=True)
        render_front_xy(
            source,
            destination,
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
# Preview WM only:
# python fig_heatmap.py \
#     --cell-types WM \
#     --window-size 720 720 \
#     --overwrite
#
# Render all four final 1440 x 1440 PNGs:
# python fig_heatmap.py --overwrite
