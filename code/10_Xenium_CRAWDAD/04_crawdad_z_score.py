#!/usr/bin/env python3
"""Plot CRAWDAD Z-score trends across AP depth and three anatomical ROIs.

Each reference-neighbor cell-type pair receives one faceted plot. The four
anatomical regions are shown in separate panels, AP depth is encoded by color,
and every AP trend is drawn with the same solid line style.

Example
-------
python 10_Xenium_CRAWDAD/04_crawdad_z_score.py \
  --input-root ../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth \
  --regions lateral dorsomedial ventromedial outside_global_roi \
  --neighborhood-distance-um 50 \
  --scale-range-um 200 1000 \
  --scale-interval-um 100 \
  --workers 4 \
  --output-dir ../processed-data/10_Xenium_CRAWDAD/module04_crawdad_z_score \
  --plot-dir ../plots/10_Xenium_CRAWDAD/module04_crawdad_z_score \
  --overwrite
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import math
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-z-score")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colormaps
from matplotlib.colors import Normalize


DEFAULT_REGIONS = (
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
REQUIRED_COLUMNS = {
    "sample",
    "depth",
    "analysis_region",
    "neighborhood_distance_um",
    "reference",
    "neighbor",
    "scale_um",
    "mean_z",
    "threshold_used",
    "reference_eligible",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Make one CRAWDAD mean-Z trend plot per directional reference-neighbor "
            "pair across AP depth and ROI."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Example\n-------\n", 1)[1],
    )
    parser.add_argument(
        "--input-root",
        type=Path,
        required=True,
        help=(
            "Module 02 root containing REGION/combined/all_trends_summary.parquet "
            "(or .csv.gz)."
        ),
    )
    parser.add_argument(
        "--regions",
        nargs="+",
        default=list(DEFAULT_REGIONS),
        help="Regions to compare [lateral dorsomedial ventromedial].",
    )
    parser.add_argument(
        "--neighborhood-distance-um",
        type=float,
        default=50,
        help="Neighborhood distance to plot [50].",
    )
    parser.add_argument(
        "--scale-range-um",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(200.0, 1000.0),
        help="Inclusive shuffle-scale range to plot [200 1000].",
    )
    parser.add_argument(
        "--scale-interval-um",
        type=float,
        default=100.0,
        help="Expected interval within --scale-range-um [100].",
    )
    parser.add_argument(
        "--references",
        nargs="+",
        default=["all"],
        help="Reference cell types to plot, or all [all].",
    )
    parser.add_argument(
        "--neighbors",
        nargs="+",
        default=["all"],
        help="Neighbor cell types to plot, or all [all].",
    )
    parser.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Independent plotting processes [1].",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=250,
        help="PNG resolution [250].",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        required=True,
        help="Directory for manifests and QC tables.",
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        required=True,
        help="Directory for PNG figures.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Replace existing plots for selected pairs.",
    )
    args = parser.parse_args(argv)
    if args.workers < 1:
        parser.error("--workers must be >= 1")
    if args.dpi < 72:
        parser.error("--dpi must be >= 72")
    if len(set(args.regions)) != len(args.regions):
        parser.error("--regions contains duplicates")
    scale_min, scale_max = args.scale_range_um
    if not np.isfinite(scale_min) or not np.isfinite(scale_max):
        parser.error("--scale-range-um values must be finite")
    if scale_min > scale_max:
        parser.error("--scale-range-um MIN must be <= MAX")
    if not np.isfinite(args.scale_interval_um) or args.scale_interval_um <= 0:
        parser.error("--scale-interval-um must be positive and finite")
    n_intervals = (scale_max - scale_min) / args.scale_interval_um
    if not np.isclose(n_intervals, round(n_intervals), atol=1e-8, rtol=0):
        parser.error(
            "--scale-range-um endpoints must align to --scale-interval-um"
        )
    args.selected_scales_um = [
        float(scale_min + index * args.scale_interval_um)
        for index in range(int(round(n_intervals)) + 1)
    ]
    return args


def _read_table(path: Path) -> pd.DataFrame:
    if path.suffix == ".parquet":
        return pd.read_parquet(path)
    if path.name.endswith(".csv.gz"):
        return pd.read_csv(path, low_memory=False)
    if path.suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported input table: {path}")


def _resolve_region_input(
    input_root: Path, region: str, neighborhood_distance_um: float
) -> Path:
    combined = (
        input_root
        / f"neighdist_{neighborhood_distance_um:g}"
        / region
        / "combined"
    )
    candidates = (
        combined / "all_trends_summary.parquet",
        combined / "all_trends_summary.csv.gz",
        combined / "all_trends_summary.csv",
    )
    for candidate in candidates:
        if candidate.exists():
            return candidate
    tried = "\n  ".join(str(x) for x in candidates)
    raise FileNotFoundError(
        f"No Module 02 combined trend summary found for region {region!r}. Tried:\n  {tried}"
    )


def _as_bool(series: pd.Series) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    normalized = series.astype(str).str.strip().str.lower()
    valid = {"true", "false", "1", "0", "yes", "no"}
    unexpected = sorted(set(normalized.dropna()) - valid - {"nan", ""})
    if unexpected:
        raise ValueError(
            "Unexpected reference_eligible values: "
            + ", ".join(unexpected[:10])
        )
    return normalized.isin({"true", "1", "yes"})


def load_inputs(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    input_qc: list[dict[str, Any]] = []
    for region in args.regions:
        path = _resolve_region_input(
            args.input_root, region, args.neighborhood_distance_um
        )
        frame = _read_table(path)
        missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing required columns: {missing}")

        frame = frame.copy()
        frame["analysis_region"] = frame["analysis_region"].astype(str)
        observed_regions = sorted(frame["analysis_region"].dropna().unique())
        if observed_regions != [region]:
            raise ValueError(
                f"{path} should contain only analysis_region={region!r}; "
                f"found {observed_regions}"
            )
        frame["reference_eligible"] = _as_bool(frame["reference_eligible"])
        for column in (
            "depth",
            "neighborhood_distance_um",
            "scale_um",
            "mean_z",
            "threshold_used",
        ):
            frame[column] = pd.to_numeric(frame[column], errors="coerce")

        selected = frame.loc[
            np.isclose(
                frame["neighborhood_distance_um"],
                args.neighborhood_distance_um,
                equal_nan=False,
            )
        ].copy()
        if selected.empty:
            available = sorted(frame["neighborhood_distance_um"].dropna().unique())
            raise ValueError(
                f"Neighborhood distance {args.neighborhood_distance_um:g} is absent "
                f"from {path}; available values: {available}"
            )

        available_scales = np.sort(selected["scale_um"].dropna().unique())
        missing_scales = [
            scale
            for scale in args.selected_scales_um
            if not np.isclose(available_scales, scale, atol=1e-8, rtol=0).any()
        ]
        if missing_scales:
            raise ValueError(
                f"Requested scales are absent from {path}: {missing_scales}; "
                f"available scales: {available_scales.tolist()}"
            )
        scale_values = selected["scale_um"].to_numpy(dtype=float)
        scale_mask = np.isclose(
            scale_values[:, np.newaxis],
            np.asarray(args.selected_scales_um)[np.newaxis, :],
            atol=1e-8,
            rtol=0,
        ).any(axis=1)
        rows_before_scale_filter = len(selected)
        selected = selected.loc[scale_mask].copy()
        excluded_scales = [
            float(scale)
            for scale in available_scales
            if not np.isclose(
                np.asarray(args.selected_scales_um), scale, atol=1e-8, rtol=0
            ).any()
        ]

        input_qc.append(
            {
                "analysis_region": region,
                "input_path": str(path.resolve()),
                "rows_selected": len(selected),
                "rows_before_scale_filter": rows_before_scale_filter,
                "rows_excluded_by_scale_filter": (
                    rows_before_scale_filter - len(selected)
                ),
                "available_scales_um": ";".join(
                    f"{value:g}" for value in available_scales
                ),
                "selected_scales_um": ";".join(
                    f"{value:g}" for value in args.selected_scales_um
                ),
                "excluded_scales_um": ";".join(
                    f"{value:g}" for value in excluded_scales
                ),
                "n_samples": selected["sample"].nunique(),
                "n_depths": selected["depth"].nunique(),
                "n_references": selected["reference"].nunique(),
                "n_neighbors": selected["neighbor"].nunique(),
                "n_scales": selected["scale_um"].nunique(),
                "n_reference_eligible_rows": int(selected["reference_eligible"].sum()),
            }
        )
        frames.append(selected)

    data = pd.concat(frames, ignore_index=True)
    duplicate_key = [
        "analysis_region",
        "sample",
        "reference",
        "neighbor",
        "scale_um",
        "neighborhood_distance_um",
    ]
    duplicates = data.duplicated(duplicate_key, keep=False)
    if duplicates.any():
        example = data.loc[duplicates, duplicate_key].head(5).to_dict("records")
        raise ValueError(f"Duplicate trend-summary rows found, for example: {example}")

    finite_thresholds = np.sort(data["threshold_used"].dropna().unique())
    if len(finite_thresholds) != 1:
        raise ValueError(
            "Expected one shared Z threshold across selected inputs; found "
            f"{finite_thresholds.tolist()}"
        )
    return data, pd.DataFrame(input_qc)


def select_levels(
    requested: Iterable[str], available: list[str], option_name: str
) -> list[str]:
    values = list(requested)
    if values == ["all"]:
        return available
    if "all" in values:
        raise ValueError(f"{option_name}: use either all or explicit cell types")
    unknown = sorted(set(values) - set(available))
    if unknown:
        raise ValueError(f"{option_name}: unknown cell types: {unknown}")
    return [value for value in available if value in values]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    if not token:
        token = "unnamed"
    return token[:80]


def validate_safe_tokens(values: list[str], label: str) -> None:
    token_to_values: dict[str, list[str]] = {}
    for value in values:
        token_to_values.setdefault(safe_token(value), []).append(value)
    collisions = {
        token: original
        for token, original in token_to_values.items()
        if len(original) > 1
    }
    if collisions:
        raise ValueError(
            f"{label} names collide after filename cleaning: {collisions}"
        )


def _region_label(region: str) -> str:
    return REGION_LABELS.get(region, region.replace("_", " ").title())


def _add_guides(axis: plt.Axes, z_threshold: float) -> None:
    axis.axhline(0, color="#303030", linewidth=0.8, alpha=0.8, zorder=0)
    axis.axhline(
        z_threshold, color="#6f6f6f", linewidth=0.8, linestyle=(0, (2, 2)), zorder=0
    )
    axis.axhline(
        -z_threshold, color="#6f6f6f", linewidth=0.8, linestyle=(0, (2, 2)), zorder=0
    )
    axis.grid(axis="x", color="#dddddd", linewidth=0.5, alpha=0.6)
    axis.set_axisbelow(True)


def _plot_curves(
    axis: plt.Axes,
    data: pd.DataFrame,
    cmap: Any,
    norm: Normalize,
) -> int:
    n_curves = 0
    grouped = data.sort_values(["depth", "scale_um"]).groupby(
        ["sample", "depth"], sort=True, dropna=False
    )
    for (_, depth), curve in grouped:
        valid = curve.loc[
            np.isfinite(curve["scale_um"]) & np.isfinite(curve["mean_z"])
        ].sort_values("scale_um")
        if valid.empty or not np.isfinite(depth):
            continue
        axis.plot(
            valid["scale_um"],
            valid["mean_z"],
            color=cmap(norm(float(depth))),
            linestyle="-",
            linewidth=1.15,
            alpha=0.86,
            solid_capstyle="round",
        )
        n_curves += 1
    return n_curves


def _save_figure(fig: plt.Figure, output_path: Path, dpi: int) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "04_crawdad_z_score.py"},
    )
    plt.close(fig)


def render_plot(
    pair_data: pd.DataFrame,
    reference: str,
    neighbor: str,
    regions: list[str],
    z_threshold: float,
    depth_min: float,
    depth_max: float,
    cmap_name: str,
    output_path: Path,
    dpi: int,
) -> int:
    cmap = colormaps[cmap_name]
    norm = Normalize(vmin=depth_min, vmax=depth_max)
    fig, axes = plt.subplots(
        1,
        len(regions),
        figsize=(5.1 * len(regions), 5.3),
        sharex=True,
        sharey=True,
        squeeze=False,
    )
    curve_count = 0
    for index, region in enumerate(regions):
        axis = axes[0, index]
        region_data = pair_data.loc[pair_data["analysis_region"] == region]
        n_region_curves = _plot_curves(axis, region_data, cmap, norm)
        curve_count += n_region_curves
        _add_guides(axis, z_threshold)
        axis.set_title(f"{_region_label(region)}  ·  n={n_region_curves} depths")
        axis.set_xlabel("Shuffle scale (µm)")
        if index == 0:
            axis.set_ylabel("Mean CRAWDAD Z score")

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    fig.subplots_adjust(left=0.055, right=0.91, bottom=0.12, top=0.82, wspace=0.08)
    colorbar_axis = fig.add_axes((0.93, 0.18, 0.012, 0.57))
    colorbar = fig.colorbar(scalar, cax=colorbar_axis)
    colorbar.set_label("AP depth (µm)")
    fig.suptitle(
        f"Reference: {reference}  →  Neighbor: {neighbor}",
        fontsize=13,
        fontweight="semibold",
        y=1.015,
    )
    fig.text(
        0.5,
        0.965,
        f"Eligible reference curves only; horizontal guides = ±{z_threshold:g}",
        ha="center",
        va="top",
        fontsize=9,
        color="#555555",
    )
    _save_figure(fig, output_path, dpi)
    return curve_count


def render_pair(task: dict[str, Any]) -> dict[str, Any]:
    data: pd.DataFrame = task["data"]
    reference = task["reference"]
    neighbor = task["neighbor"]
    regions = task["regions"]
    output_path = Path(task["output_path"])
    overwrite = task["overwrite"]

    pair = data.loc[
        (data["reference"] == reference)
        & (data["neighbor"] == neighbor)
        & data["reference_eligible"]
        & np.isfinite(data["mean_z"])
    ].copy()

    expected_curves = int(task["expected_curves"])
    eligible_reference_curves = int(task["eligible_reference_curves"])
    plotted_groups = pair[["analysis_region", "sample", "depth"]].drop_duplicates()
    plotted_curves = len(plotted_groups)

    result: dict[str, Any] = {
        "reference": reference,
        "neighbor": neighbor,
        "expected_curves": expected_curves,
        "eligible_reference_curves": eligible_reference_curves,
        "plotted_curves": plotted_curves,
        "reference_ineligible_curves": expected_curves - eligible_reference_curves,
        "eligible_but_missing_result_curves": max(
            0, eligible_reference_curves - plotted_curves
        ),
        "n_scales_present": int(pair["scale_um"].nunique()),
        "max_abs_mean_z": (
            float(pair["mean_z"].abs().max()) if not pair.empty else math.nan
        ),
        "plot": str(output_path.resolve()),
    }

    if output_path.exists() and not overwrite:
        result["status"] = "existing_not_overwritten"
        return result
    if pair.empty:
        fig, axis = plt.subplots(figsize=(9, 5))
        axis.axis("off")
        axis.text(
            0.5,
            0.55,
            f"Reference: {reference}  →  Neighbor: {neighbor}",
            ha="center",
            va="center",
            fontsize=13,
            fontweight="semibold",
        )
        axis.text(
            0.5,
            0.43,
            "No eligible reference curve with a finite mean Z score.",
            ha="center",
            va="center",
            fontsize=11,
            color="#555555",
        )
        _save_figure(fig, output_path, task["dpi"])
        result["status"] = "placeholder_no_eligible_data"
        return result

    render_plot(
        pair_data=pair,
        reference=reference,
        neighbor=neighbor,
        regions=regions,
        z_threshold=task["z_threshold"],
        depth_min=task["depth_min"],
        depth_max=task["depth_max"],
        cmap_name=task["cmap_name"],
        output_path=output_path,
        dpi=task["dpi"],
    )
    result["status"] = "written"
    return result


def _jsonable_args(args: argparse.Namespace) -> dict[str, Any]:
    converted: dict[str, Any] = {}
    for key, value in vars(args).items():
        if isinstance(value, Path):
            converted[key] = str(value.resolve())
        else:
            converted[key] = value
    return converted


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data, input_qc = load_inputs(args)

    available_references = sorted(data["reference"].dropna().astype(str).unique())
    available_neighbors = sorted(data["neighbor"].dropna().astype(str).unique())
    references = select_levels(args.references, available_references, "--references")
    neighbors = select_levels(args.neighbors, available_neighbors, "--neighbors")
    if not references or not neighbors:
        raise ValueError("No reference-neighbor pairs selected")
    validate_safe_tokens(references, "Reference")
    validate_safe_tokens(neighbors, "Neighbor")

    depths = np.sort(data["depth"].dropna().unique())
    if len(depths) < 2:
        raise ValueError("At least two AP depths are required for the color scale")
    z_threshold = float(data["threshold_used"].dropna().iloc[0])

    distance_token = f"neighdist_{args.neighborhood_distance_um:g}"
    output_root = args.output_dir / distance_token
    plot_root = args.plot_dir / distance_token
    output_root.mkdir(parents=True, exist_ok=True)
    plot_root.mkdir(parents=True, exist_ok=True)

    tasks: list[dict[str, Any]] = []
    expected_curves = int(
        data[["analysis_region", "sample"]].drop_duplicates().shape[0]
    )
    eligibility_by_reference = (
        data.groupby(["reference", "analysis_region", "sample"], as_index=False)[
            "reference_eligible"
        ]
        .max()
        .groupby("reference")["reference_eligible"]
        .sum()
        .astype(int)
        .to_dict()
    )
    for reference in references:
        reference_token = safe_token(reference)
        for neighbor in neighbors:
            neighbor_token = safe_token(neighbor)
            output_path = (
                plot_root
                / "by_reference"
                / reference_token
                / f"{reference_token}__to__{neighbor_token}.png"
            )
            tasks.append(
                {
                    "data": data.loc[
                        (data["reference"] == reference)
                        & (data["neighbor"] == neighbor)
                    ].copy(),
                    "reference": reference,
                    "neighbor": neighbor,
                    "regions": list(args.regions),
                    "expected_curves": expected_curves,
                    "eligible_reference_curves": eligibility_by_reference.get(
                        reference, 0
                    ),
                    "output_path": str(output_path),
                    "overwrite": args.overwrite,
                    "dpi": args.dpi,
                    "z_threshold": z_threshold,
                    "depth_min": float(depths.min()),
                    "depth_max": float(depths.max()),
                    "cmap_name": "viridis",
                }
            )

    print(
        f"Plotting {len(references)} references × {len(neighbors)} neighbors "
        f"= {len(tasks)} directional pairs; workers={args.workers}",
        flush=True,
    )
    results: list[dict[str, Any]] = []
    if args.workers == 1:
        for index, task in enumerate(tasks, start=1):
            results.append(render_pair(task))
            if index == 1 or index % 25 == 0 or index == len(tasks):
                print(f"Completed {index}/{len(tasks)} pairs", flush=True)
    else:
        with concurrent.futures.ProcessPoolExecutor(
            max_workers=args.workers
        ) as executor:
            futures = [executor.submit(render_pair, task) for task in tasks]
            for index, future in enumerate(
                concurrent.futures.as_completed(futures), start=1
            ):
                results.append(future.result())
                if index == 1 or index % 25 == 0 or index == len(tasks):
                    print(f"Completed {index}/{len(tasks)} pairs", flush=True)

    manifest = pd.DataFrame(results).sort_values(
        ["reference", "neighbor"], kind="stable"
    )
    manifest_path = output_root / "plot_manifest.csv"
    input_qc_path = output_root / "input_qc.csv"
    args_path = output_root / "resolved_arguments.json"
    manifest.to_csv(manifest_path, index=False)
    input_qc.to_csv(input_qc_path, index=False)
    with args_path.open("w", encoding="utf-8") as handle:
        json.dump(
            {
                "arguments": _jsonable_args(args),
                "resolved_references": references,
                "resolved_neighbors": neighbors,
                "ap_depths_um": [float(value) for value in depths],
                "observed_scales_um": [
                    float(value) for value in np.sort(data["scale_um"].unique())
                ],
                "z_threshold": z_threshold,
                "n_directional_pairs": len(tasks),
                "plot_structure": "one faceted PNG per directional pair",
            },
            handle,
            indent=2,
        )
        handle.write("\n")

    print(f"Plots: {plot_root.resolve()}", flush=True)
    print(f"Manifest: {manifest_path.resolve()}", flush=True)
    print(f"Input QC: {input_qc_path.resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
