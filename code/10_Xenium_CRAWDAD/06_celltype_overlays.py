#!/usr/bin/env python3
"""Plot selected cell-type overlays within one CRAWDAD ROI.

All cells in the selected sample/ROI are retained as spatial context. Cell
types assigned to an overlay group use the established CRAWDAD color
palette; every other cell is drawn in gray beneath them.
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-celltype-overlays")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


CELLTYPE_COLORS = {
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
    "WM": "#FFA500",
}

DEFAULT_GROUPS = (
    (
        "fibroblast_microglia",
        ("Fibroblast_A", "Fibroblast_B", "Microglia_B"),
    ),
    (
        "msn_island_astro",
        ("DRD1_MSN", "DRD2_MSN", "MSN_Oligo", "D1_Island_A", "Astro_B"),
    ),
    (
        "astro_wm_island_oligo",
        ("Astro_A", "WM", "D1_Island_B", "Astrocyte_Oligo"),
    ),
)


def parse_group(value: str) -> tuple[str, tuple[str, ...]]:
    if "=" not in value:
        raise argparse.ArgumentTypeError(
            "--group must use NAME=CELLTYPE[,CELLTYPE...]"
        )
    name, raw_types = value.split("=", 1)
    name = name.strip()
    celltypes = tuple(item.strip() for item in raw_types.split(",") if item.strip())
    if not name or not celltypes:
        raise argparse.ArgumentTypeError(
            "--group requires a nonempty name and at least one cell type"
        )
    if len(set(celltypes)) != len(celltypes):
        raise argparse.ArgumentTypeError(f"--group {name!r} contains duplicates")
    return name, celltypes


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input-root",
        required=True,
        type=Path,
        help="Module 01 directory containing crawdad_inputs/SAMPLE/ROI.csv.gz.",
    )
    parser.add_argument(
        "--samples", nargs="+", default=["Br6660_Nac10_4080"]
    )
    parser.add_argument("--regions", nargs="+", default=["dorsomedial"])
    parser.add_argument(
        "--group",
        action="append",
        type=parse_group,
        help=(
            "Overlay definition NAME=CELLTYPE[,CELLTYPE...]. Repeat for multiple "
            "figures. The three default groups are used when omitted."
        ),
    )
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--plot-dir", required=True, type=Path)
    parser.add_argument("--background-color", default="#AFAFAF")
    parser.add_argument("--background-size", type=float, default=0.55)
    parser.add_argument("--highlight-size", type=float, default=2.2)
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    args.groups = tuple(args.group) if args.group else DEFAULT_GROUPS
    if args.background_size <= 0 or args.highlight_size <= 0:
        parser.error("Point sizes must be positive")
    if args.dpi < 72:
        parser.error("--dpi must be >= 72")
    group_names = [name for name, _ in args.groups]
    if len(set(group_names)) != len(group_names):
        parser.error("Overlay group names must be unique")
    if len(set(args.samples)) != len(args.samples):
        parser.error("--samples contains duplicates")
    if len(set(args.regions)) != len(args.regions):
        parser.error("--regions contains duplicates")
    return args


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    if not token:
        raise ValueError(f"Cannot make a filename token from {value!r}")
    return token


def read_cells(
    input_root: Path, samples: list[str], regions: list[str]
) -> tuple[pd.DataFrame, Path]:
    path = input_root / "cell_assignments" / "all_cells_roi_assignment.parquet"
    if not path.exists():
        raise FileNotFoundError(f"Missing Module 01 cell assignment parquet: {path}")
    columns = [
        "cell_id", "sample", "cell_type", "x_aligned", "y_aligned",
        "primary_roi",
    ]
    cells = pd.read_parquet(
        path,
        columns=columns,
        filters=[("sample", "in", samples), ("primary_roi", "in", regions)],
    ).rename(
        columns={
            "cell_type": "celltype",
            "x_aligned": "x",
            "y_aligned": "y",
            "primary_roi": "roi",
        }
    )
    required = {"cell_id", "x", "y", "celltype", "sample", "roi"}
    missing = sorted(required - set(cells.columns))
    if missing:
        raise ValueError(f"{path} is missing columns: {missing}")
    cells = cells.copy()
    cells["sample"] = cells["sample"].astype(str)
    cells["roi"] = cells["roi"].astype(str)
    cells["celltype"] = cells["celltype"].astype(str)
    cells["x"] = pd.to_numeric(cells["x"], errors="coerce")
    cells["y"] = pd.to_numeric(cells["y"], errors="coerce")
    missing_samples = sorted(set(samples) - set(cells["sample"].unique()))
    missing_regions = sorted(set(regions) - set(cells["roi"].unique()))
    if missing_samples:
        raise ValueError(f"Samples absent from {path}: {missing_samples}")
    if missing_regions:
        raise ValueError(f"Regions absent from {path}: {missing_regions}")
    if not np.isfinite(cells[["x", "y"]].to_numpy(dtype=float)).all():
        raise ValueError(f"Non-finite aligned coordinates in {path}")
    if cells["cell_id"].duplicated().any():
        raise ValueError(f"Duplicate cell IDs in {path}")
    observed_tasks = set(zip(cells["sample"], cells["roi"]))
    missing_tasks = [
        f"{sample}/{region}"
        for sample in samples
        for region in regions
        if (sample, region) not in observed_tasks
    ]
    if missing_tasks:
        raise ValueError(f"Missing sample-region assignments: {missing_tasks}")
    return cells, path


def atomic_save(fig: plt.Figure, target: Path, dpi: int) -> None:
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(f".{target.stem}.tmp.{os.getpid()}.png")
    fig.savefig(
        temporary,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "06_celltype_overlays.py"},
    )
    plt.close(fig)
    temporary.replace(target)


def render_group(
    cells: pd.DataFrame,
    group_index: int,
    group_name: str,
    selected_types: tuple[str, ...],
    args: argparse.Namespace,
    sample: str,
    region: str,
    plot_root: Path,
    limits: tuple[float, float, float, float],
) -> dict[str, Any]:
    unknown_palette = sorted(set(selected_types) - set(CELLTYPE_COLORS))
    if unknown_palette:
        raise ValueError(f"Cell types missing from the CRAWDAD palette: {unknown_palette}")
    absent = sorted(set(selected_types) - set(cells["celltype"].unique()))
    if absent:
        raise ValueError(
            f"Requested cell types are absent from {sample}/{region}: {absent}"
        )

    token = safe_token(group_name)
    target = plot_root / f"overlay{group_index:02d}_{token}.png"
    if target.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists; use --overwrite: {target}")

    highlighted = cells["celltype"].isin(selected_types)
    background = cells.loc[~highlighted]
    fig, axis = plt.subplots(figsize=(10, 9))
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
    legend_handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", color=args.background_color,
            markersize=5, label=f"Other ROI cells (n={len(background):,})",
        )
    ]
    selected_counts: dict[str, int] = {}
    for celltype in selected_types:
        subset = cells.loc[cells["celltype"].eq(celltype)]
        selected_counts[celltype] = len(subset)
        axis.scatter(
            subset["x"],
            subset["y"],
            s=args.highlight_size,
            c=CELLTYPE_COLORS[celltype],
            alpha=0.9,
            linewidths=0,
            rasterized=True,
            zorder=2,
        )
        legend_handles.append(
            Line2D(
                [0], [0], marker="o", linestyle="",
                color=CELLTYPE_COLORS[celltype], markersize=6,
                label=f"{celltype} (n={len(subset):,})",
            )
        )

    min_x, max_x, min_y, max_y = limits
    axis.set_xlim(min_x, max_x)
    axis.set_ylim(min_y, max_y)
    axis.set_aspect("equal")
    axis.set_xlabel("Aligned x (µm)")
    axis.set_ylabel("Aligned y (µm)")
    axis.set_title(
        f"{sample} | {region}\nSelected cell-type overlay {group_index}"
    )
    axis.legend(
        handles=legend_handles,
        loc="center left",
        bbox_to_anchor=(1.02, 0.5),
        frameon=False,
    )
    axis.grid(False)
    atomic_save(fig, target, args.dpi)
    return {
        "group_index": group_index,
        "group_name": group_name,
        "celltypes": ";".join(selected_types),
        "n_roi_cells": len(cells),
        "n_background_cells": len(background),
        "n_highlighted_cells": int(highlighted.sum()),
        "highlighted_counts_json": json.dumps(selected_counts, sort_keys=True),
        "plot": str(target.resolve()),
        "status": "written",
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    all_cells, input_path = read_cells(args.input_root, args.samples, args.regions)
    all_rows: list[dict[str, Any]] = []
    all_counts: list[pd.DataFrame] = []
    n_tasks = len(args.samples) * len(args.regions)
    completed_tasks = 0
    for sample in args.samples:
        for region in args.regions:
            cells = all_cells.loc[
                all_cells["sample"].eq(sample) & all_cells["roi"].eq(region)
            ].copy()
            output_root = args.output_dir / sample / region
            plot_root = args.plot_dir / sample / region
            output_root.mkdir(parents=True, exist_ok=True)
            plot_root.mkdir(parents=True, exist_ok=True)

            x_span = max(float(cells["x"].max() - cells["x"].min()), 1.0)
            y_span = max(float(cells["y"].max() - cells["y"].min()), 1.0)
            limits = (
                float(cells["x"].min() - 0.02 * x_span),
                float(cells["x"].max() + 0.02 * x_span),
                float(cells["y"].min() - 0.02 * y_span),
                float(cells["y"].max() + 0.02 * y_span),
            )

            rows = []
            for index, (name, celltypes) in enumerate(args.groups, start=1):
                row = render_group(
                    cells, index, name, celltypes, args, sample, region,
                    plot_root, limits
                )
                row = {"sample": sample, "region": region, **row}
                rows.append(row)
                all_rows.append(row)
            pd.DataFrame(rows).to_csv(
                output_root / "plot_manifest.csv", index=False
            )
            counts = (
                cells["celltype"]
                .value_counts()
                .rename_axis("celltype")
                .reset_index(name="n_cells")
                .sort_values("celltype")
            )
            counts.to_csv(output_root / "celltype_counts.csv", index=False)
            all_counts.append(
                counts.assign(sample=sample, region=region)[
                    ["sample", "region", "celltype", "n_cells"]
                ]
            )
            metadata = {
                "sample": sample,
                "region": region,
                "groups": [
                    {"name": name, "celltypes": list(celltypes)}
                    for name, celltypes in args.groups
                ],
                "input_file": str(input_path.resolve()),
                "coordinate_definition": "Module 01 aligned Xenium x/y",
                "n_roi_cells": len(cells),
                "plot_limits": {
                    "x_min": limits[0], "x_max": limits[1],
                    "y_min": limits[2], "y_max": limits[3],
                },
            }
            with (output_root / "resolved_arguments.json").open(
                "w", encoding="utf-8"
            ) as handle:
                json.dump(metadata, handle, indent=2)
                handle.write("\n")
            completed_tasks += 1
            print(
                f"Completed {completed_tasks}/{n_tasks}: {sample}/{region} "
                f"({len(cells):,} cells)",
                flush=True,
            )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_rows).to_csv(
        args.output_dir / "plot_manifest.csv", index=False
    )
    pd.concat(all_counts, ignore_index=True).to_csv(
        args.output_dir / "celltype_counts_all_tasks.csv", index=False
    )
    print(f"Input ROI cells across tasks: {len(all_cells):,}")
    print(f"Plots written: {len(all_rows)}")
    print(f"Plot directory: {args.plot_dir.resolve()}")
    print(f"Manifest: {(args.output_dir / 'plot_manifest.csv').resolve()}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
