#!/usr/bin/env python3
"""Draw six supplemental CRAWDAD relationship summaries from Module 02."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path


SCRIPT_DIR = Path(__file__).resolve().parent
CODE_ROOT = SCRIPT_DIR.parents[1]
PROJECT_ROOT = CODE_ROOT.parent
DEFAULT_INPUT_ROOT = (
    PROJECT_ROOT
    / "processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth"
    / "neighdist_50"
)
DEFAULT_PLOT_DIR = (
    PROJECT_ROOT / "plots/10_Xenium_CRAWDAD/figure_plots/supp_fig"
)
DEFAULT_RSCRIPT = Path(
    "/jhpce/shared/community/core/conda_R/4.5/R/bin/Rscript"
)
RENDERER = SCRIPT_DIR / "03_relationship_summary_overlay_groups.R"

DEPTHS = (580, 1090, 1580, 2080, 2580, 3080, 3580, 4080, 4580, 5080, 5580)
REQUESTED_PLOTS = (
    (1580, "dorsomedial"),
    (4080, "dorsomedial"),
    (4080, "lateral"),
    (4080, "ventromedial"),
    (4080, "outside_global_roi"),
    (5580, "dorsomedial"),
)
REGION_LABELS = {
    "lateral": "Lateral",
    "dorsomedial": "Dorsomedial",
    "ventromedial": "Ventromedial",
    "outside_global_roi": "Outside global ROI",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=DEFAULT_INPUT_ROOT)
    parser.add_argument("--plot-dir", type=Path, default=DEFAULT_PLOT_DIR)
    parser.add_argument("--rscript", type=Path, default=DEFAULT_RSCRIPT)
    parser.add_argument(
        "--scale-range-um",
        nargs=2,
        type=float,
        default=(200.0, 1000.0),
        metavar=("MIN", "MAX"),
    )
    parser.add_argument("--scale-interval-um", type=float, default=100.0)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()
    if args.scale_range_um[0] > args.scale_range_um[1]:
        parser.error("--scale-range-um MIN must be <= MAX")
    if args.scale_interval_um <= 0:
        parser.error("--scale-interval-um must be positive")
    return args


def find_task(input_root: Path, depth: int, region: str) -> Path:
    task_root = input_root / region / "tasks"
    if not task_root.is_dir():
        raise FileNotFoundError(task_root)
    matches = []
    for metadata_path in task_root.glob("*/task_metadata.json"):
        metadata = json.loads(metadata_path.read_text())
        if (
            int(metadata.get("depth", -1)) == depth
            and metadata.get("analysis_region") == region
        ):
            matches.append(metadata_path.parent)
    if len(matches) != 1:
        raise ValueError(
            f"Expected one Module 02 task for depth={depth}, region={region}; "
            f"found {len(matches)}"
        )
    task_dir = matches[0]
    for filename in ("trends_permutation.csv.gz", "task_identity.json"):
        if not (task_dir / filename).is_file():
            raise FileNotFoundError(task_dir / filename)
    return task_dir


def main() -> int:
    args = parse_args()
    if not args.rscript.is_file():
        raise FileNotFoundError(args.rscript)
    if not RENDERER.is_file():
        raise FileNotFoundError(RENDERER)

    jobs = []
    for depth, region in REQUESTED_PLOTS:
        slice_number = DEPTHS.index(depth) + 1
        task_dir = find_task(args.input_root, depth, region)
        target = args.plot_dir / f"relationship_summary_{depth}_{region}.png"
        title = f"Slice {slice_number} | {REGION_LABELS[region]}"
        jobs.append((task_dir, target, title))

    existing = [target for _, target, _ in jobs if target.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            f"{len(existing)} output plots already exist; use --overwrite"
        )
    args.plot_dir.mkdir(parents=True, exist_ok=True)

    environment = os.environ.copy()
    environment.setdefault("R_LIBS_USER", "/users/jyao/R/4.5")
    environment.setdefault(
        "XDG_CACHE_HOME", str(Path(os.environ.get("TMPDIR", "/tmp")) / "crawdad-supp-fig-cache")
    )
    for task_dir, target, title in jobs:
        command = [
            str(args.rscript),
            str(RENDERER),
            "--task-dir",
            str(task_dir),
            "--output",
            str(target),
            "--scale-range-um",
            f"{args.scale_range_um[0]:g}",
            f"{args.scale_range_um[1]:g}",
            "--scale-interval-um",
            f"{args.scale_interval_um:g}",
            "--title",
            title,
            "--no-boxes",
        ]
        if args.overwrite:
            command.append("--overwrite")
        subprocess.run(command, check=True, env=environment)

    print(f"Saved {len(jobs)} plots under {args.plot_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


# Example usage from the code directory:
# /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/bin/python \
#   10_Xenium_CRAWDAD/figure_plots/03_relationship_summary.py
# Add --overwrite when replacing existing plots.
