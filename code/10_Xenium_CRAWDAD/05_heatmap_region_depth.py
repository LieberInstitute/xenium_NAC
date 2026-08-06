#!/usr/bin/env python3
"""Plot Module 02 CRAWDAD relationships across AP depth and region.

One PNG is created for every directional reference-to-neighbor pair. The
x-axis is Br6660 AP depth and the y-axis is anatomical region. As in Module 02
``relationship_summary.png``, a point is drawn only when the relationship is
significant at at least one scale. Point color is the mean Z score at the first
(smallest) significant scale, and point size is inversely related to that scale.

Example
-------
python 10_Xenium_CRAWDAD/05_heatmap_region_depth.py \
  --input-root ../processed-data/10_Xenium_CRAWDAD/module02_crawdad_within_depth \
  --regions lateral dorsomedial ventromedial outside_global_roi \
  --neighborhood-distance-um 50 \
  --scale-range-um 200 1000 \
  --scale-interval-um 100 \
  --workers 4 \
  --output-dir ../processed-data/10_Xenium_CRAWDAD/module05_heatmap_region_depth \
  --plot-dir ../plots/10_Xenium_CRAWDAD/module05_heatmap_region_depth \
  --overwrite
"""

from __future__ import annotations

import argparse
import concurrent.futures
import json
import os
import re
import sys
from pathlib import Path
from typing import Any, Iterable

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-crawdad-region-depth")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.lines import Line2D


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
    "donor",
    "sample",
    "depth",
    "analysis_region",
    "neighborhood_distance_um",
    "reference",
    "neighbor",
    "scale_um",
    "mean_z",
    "reference_eligible",
    "threshold_used",
}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Plot one AP-depth-by-region CRAWDAD relationship dot heatmap per "
            "directional cell-type pair."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__.split("Example\n-------\n", 1)[1],
    )
    parser.add_argument("--input-root", required=True, type=Path)
    parser.add_argument("--regions", nargs="+", default=list(DEFAULT_REGIONS))
    parser.add_argument("--neighborhood-distance-um", type=float, default=50)
    parser.add_argument("--references", nargs="+", default=["all"])
    parser.add_argument("--neighbors", nargs="+", default=["all"])
    parser.add_argument("--z-score-limit", type=float, default=10)
    parser.add_argument(
        "--scale-range-um",
        nargs=2,
        type=float,
        metavar=("MIN", "MAX"),
        default=(200.0, 1000.0),
    )
    parser.add_argument("--scale-interval-um", type=float, default=100)
    parser.add_argument("--workers", type=int, default=1)
    parser.add_argument("--dpi", type=int, default=250)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--plot-dir", required=True, type=Path)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args(argv)
    if len(set(args.regions)) != len(args.regions):
        parser.error("--regions contains duplicates")
    if args.workers < 1:
        parser.error("--workers must be >= 1")
    if args.dpi < 72:
        parser.error("--dpi must be >= 72")
    if not np.isfinite(args.neighborhood_distance_um) or args.neighborhood_distance_um <= 0:
        parser.error("--neighborhood-distance-um must be positive")
    if not np.isfinite(args.z_score_limit) or args.z_score_limit <= 0:
        parser.error("--z-score-limit must be positive")
    if not np.isfinite(args.scale_interval_um) or args.scale_interval_um <= 0:
        parser.error("--scale-interval-um must be positive")
    scale_min, scale_max = args.scale_range_um
    if not np.isfinite(scale_min) or not np.isfinite(scale_max):
        parser.error("--scale-range-um values must be finite")
    if scale_min > scale_max:
        parser.error("--scale-range-um MIN must be <= MAX")
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
    if path.name.endswith(".csv.gz") or path.suffix == ".csv":
        return pd.read_csv(path, low_memory=False)
    raise ValueError(f"Unsupported input table: {path}")


def _resolve_input(
    input_root: Path, region: str, neighborhood_distance_um: float
) -> Path:
    root = (
        input_root
        / f"neighdist_{neighborhood_distance_um:g}"
        / region
        / "combined"
    )
    candidates = (
        root / "all_trends_summary.parquet",
        root / "all_trends_summary.csv.gz",
        root / "all_trends_summary.csv",
    )
    existing = [path for path in candidates if path.exists()]
    if not existing:
        raise FileNotFoundError(
            f"No Module 02 trend-summary table found for {region}: "
            + ", ".join(str(path) for path in candidates)
        )
    return existing[0]


def _resolve_eligibility_input(
    input_root: Path, region: str, neighborhood_distance_um: float
) -> Path:
    path = (
        input_root
        / f"neighdist_{neighborhood_distance_um:g}"
        / region
        / "qc"
        / "reference_eligibility.csv"
    )
    if not path.exists():
        raise FileNotFoundError(f"Missing Module 02 reference eligibility table: {path}")
    return path


def _as_bool(series: pd.Series, column: str) -> pd.Series:
    if pd.api.types.is_bool_dtype(series):
        return series.fillna(False)
    values = series.astype(str).str.strip().str.lower()
    allowed = {"true", "false", "1", "0", "yes", "no", "nan", ""}
    unexpected = sorted(set(values) - allowed)
    if unexpected:
        raise ValueError(f"Unexpected {column} values: {unexpected[:10]}")
    return values.isin({"true", "1", "yes"})


def load_inputs(
    args: argparse.Namespace,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    frames: list[pd.DataFrame] = []
    eligibility_frames: list[pd.DataFrame] = []
    qc_rows: list[dict[str, Any]] = []
    for region in args.regions:
        path = _resolve_input(
            args.input_root, region, args.neighborhood_distance_um
        )
        eligibility_path = _resolve_eligibility_input(
            args.input_root, region, args.neighborhood_distance_um
        )
        frame = _read_table(path)
        missing = sorted(REQUIRED_COLUMNS - set(frame.columns))
        if missing:
            raise ValueError(f"{path} is missing columns: {missing}")
        frame = frame.copy()
        if sorted(frame["analysis_region"].dropna().astype(str).unique()) != [region]:
            raise ValueError(f"Input region mismatch in {path}")
        frame["analysis_region"] = frame["analysis_region"].astype(str)
        frame["reference"] = frame["reference"].astype(str)
        frame["neighbor"] = frame["neighbor"].astype(str)
        frame["sample"] = frame["sample"].astype(str)
        frame["reference_eligible"] = _as_bool(
            frame["reference_eligible"], "reference_eligible"
        )
        for column in (
            "depth",
            "neighborhood_distance_um",
            "scale_um",
            "mean_z",
            "threshold_used",
        ):
            frame[column] = pd.to_numeric(frame[column], errors="coerce")
        frame = frame.loc[
            np.isclose(
                frame["neighborhood_distance_um"],
                args.neighborhood_distance_um,
                equal_nan=False,
            )
        ].copy()
        if frame.empty:
            raise ValueError(
                f"Neighborhood distance {args.neighborhood_distance_um:g} is absent in {path}"
            )
        available_scales = np.sort(frame["scale_um"].dropna().unique())
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
        scale_values = frame["scale_um"].to_numpy(dtype=float)
        scale_mask = np.isclose(
            scale_values[:, np.newaxis],
            np.asarray(args.selected_scales_um)[np.newaxis, :],
            atol=1e-8,
            rtol=0,
        ).any(axis=1)
        rows_before_scale_filter = len(frame)
        frame = frame.loc[scale_mask].copy()
        excluded_scales = [
            float(scale)
            for scale in available_scales
            if not np.isclose(
                np.asarray(args.selected_scales_um), scale, atol=1e-8, rtol=0
            ).any()
        ]
        pair_key = [
            "donor", "sample", "depth", "analysis_region",
            "neighborhood_distance_um", "reference", "neighbor",
        ]
        scale_key = pair_key + ["scale_um"]
        if frame.duplicated(scale_key, keep=False).any():
            example = frame.loc[
                frame.duplicated(scale_key, keep=False), scale_key
            ].head()
            raise ValueError(
                f"Duplicate directional-pair scale rows in {path}: "
                f"{example.to_dict('records')}"
            )
        consistency = frame.groupby(pair_key, dropna=False).agg(
            n_reference_eligible_values=("reference_eligible", "nunique"),
            n_threshold_values=("threshold_used", "nunique"),
        )
        if (
            (consistency["n_reference_eligible_values"] != 1).any()
            or (consistency["n_threshold_values"] != 1).any()
        ):
            raise ValueError(
                f"Reference eligibility or threshold changes across scales in {path}"
            )
        pair_summary = (
            frame.sort_values("scale_um")
            .groupby(pair_key, as_index=False, dropna=False)
            .agg(
                reference_eligible=("reference_eligible", "first"),
                threshold_used=("threshold_used", "first"),
            )
        )
        significant = frame.loc[
            frame["reference_eligible"]
            & np.isfinite(frame["mean_z"])
            & (frame["mean_z"].abs() >= frame["threshold_used"])
        ].sort_values("scale_um")
        first_hits = significant.drop_duplicates(pair_key, keep="first")[
            pair_key + ["scale_um", "mean_z"]
        ].rename(
            columns={
                "scale_um": "first_significant_scale_um",
                "mean_z": "z_at_first_significant_scale",
            }
        )
        pair_summary = pair_summary.merge(
            first_hits, on=pair_key, how="left", validate="one_to_one"
        )
        pair_summary["ever_significant"] = pair_summary[
            "first_significant_scale_um"
        ].notna()
        pair_summary["relationship_at_first_significant_scale"] = np.select(
            [
                pair_summary["z_at_first_significant_scale"] > 0,
                pair_summary["z_at_first_significant_scale"] < 0,
            ],
            ["enrichment", "depletion"],
            default=None,
        )
        pair_summary = pair_summary.loc[
            pair_summary["reference_eligible"]
        ].copy()

        eligibility = pd.read_csv(eligibility_path, low_memory=False)
        eligibility_required = {
            "sample", "depth", "analysis_region", "cell_type",
            "present", "reference_eligible",
        }
        eligibility_missing = sorted(eligibility_required - set(eligibility.columns))
        if eligibility_missing:
            raise ValueError(
                f"{eligibility_path} is missing columns: {eligibility_missing}"
            )
        eligibility = eligibility.copy()
        eligibility["sample"] = eligibility["sample"].astype(str)
        eligibility["analysis_region"] = eligibility["analysis_region"].astype(str)
        eligibility["reference"] = eligibility["cell_type"].astype(str)
        eligibility["depth"] = pd.to_numeric(eligibility["depth"], errors="coerce")
        eligibility["present"] = _as_bool(eligibility["present"], "present")
        eligibility["reference_eligible"] = _as_bool(
            eligibility["reference_eligible"], "reference_eligible"
        )
        if sorted(eligibility["analysis_region"].dropna().unique()) != [region]:
            raise ValueError(f"Eligibility region mismatch in {eligibility_path}")
        eligibility_key = ["analysis_region", "sample", "reference"]
        if eligibility.duplicated(eligibility_key, keep=False).any():
            raise ValueError(f"Duplicate reference eligibility rows in {eligibility_path}")
        eligibility_frames.append(
            eligibility[
                eligibility_key + ["depth", "present", "reference_eligible"]
            ]
        )
        qc_rows.append(
            {
                "analysis_region": region,
                "input_path": str(path.resolve()),
                "eligibility_path": str(eligibility_path.resolve()),
                "n_summary_rows_before_scale_filter": rows_before_scale_filter,
                "n_summary_rows_after_scale_filter": len(frame),
                "n_summary_rows_excluded_by_scale_filter": (
                    rows_before_scale_filter - len(frame)
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
                "n_pair_rows": len(pair_summary),
                "n_samples": pair_summary["sample"].nunique(),
                "n_depths": pair_summary["depth"].nunique(),
                "n_references": pair_summary["reference"].nunique(),
                "n_neighbors": pair_summary["neighbor"].nunique(),
                "n_significant_rows": int(pair_summary["ever_significant"].sum()),
            }
        )
        frames.append(pair_summary)

    data = pd.concat(frames, ignore_index=True)
    thresholds = np.sort(data["threshold_used"].dropna().unique())
    if len(thresholds) != 1:
        raise ValueError(f"Expected one shared threshold; found {thresholds.tolist()}")
    mappings = data[["sample", "depth"]].drop_duplicates()
    if mappings["sample"].duplicated().any():
        raise ValueError("Sample-to-depth mapping is not one-to-one")
    region_depths = data.groupby("analysis_region")["depth"].apply(
        lambda x: tuple(sorted(x.dropna().unique()))
    )
    if len(set(region_depths)) != 1:
        raise ValueError("AP-depth coverage differs across regions")
    eligibility_data = pd.concat(eligibility_frames, ignore_index=True)
    return data, pd.DataFrame(qc_rows), eligibility_data


def select_levels(
    requested: Iterable[str], available: list[str], option_name: str
) -> list[str]:
    values = list(requested)
    if values == ["all"]:
        return available
    if "all" in values:
        raise ValueError(f"{option_name}: use either all or explicit values")
    unknown = sorted(set(values) - set(available))
    if unknown:
        raise ValueError(f"{option_name}: unknown cell types: {unknown}")
    return [value for value in available if value in values]


def safe_token(value: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._-]+", "_", value.strip()).strip("._")
    return (token or "unnamed")[:80]


def validate_safe_tokens(values: list[str], label: str) -> None:
    tokens: dict[str, list[str]] = {}
    for value in values:
        tokens.setdefault(safe_token(value), []).append(value)
    collisions = {key: value for key, value in tokens.items() if len(value) > 1}
    if collisions:
        raise ValueError(f"{label} filename-token collisions: {collisions}")


def _radius(scale: np.ndarray | float, levels: np.ndarray) -> np.ndarray:
    values = np.asarray(scale, dtype=float)
    if len(levels) == 1:
        return np.full_like(values, 6.0)
    indices = np.array(
        [int(np.argmin(np.abs(levels - value))) for value in values.ravel()]
    ).reshape(values.shape)
    return 10.0 - indices * (8.0 / (len(levels) - 1))


def _save_figure(fig: plt.Figure, output: Path, dpi: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(f".{output.stem}.tmp.{os.getpid()}.png")
    fig.savefig(
        temporary,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        metadata={"Software": "05_heatmap_region_depth.py"},
    )
    plt.close(fig)
    temporary.replace(output)


def render_pair(task: dict[str, Any]) -> dict[str, Any]:
    pair = task["data"].copy()
    output = Path(task["output"])
    if output.exists() and not task["overwrite"]:
        return {
            "reference": task["reference"],
            "neighbor": task["neighbor"],
            "status": "existing_not_overwritten",
            "n_expected_positions": task["n_expected_positions"],
            "n_eligible_positions": int(pair["reference_eligible"].sum()),
            "n_significant_points": int(pair["ever_significant"].sum()),
            "plot": str(output.resolve()),
        }

    regions = task["regions"]
    significant = pair.loc[pair["ever_significant"]].copy()
    cmap = LinearSegmentedColormap.from_list(
        "crawdad_relationship", ["blue", "white", "red"]
    )
    norm = Normalize(vmin=-task["z_score_limit"], vmax=task["z_score_limit"])
    fig, axis = plt.subplots(figsize=(13.5, 5.2))
    y_positions = {region: index for index, region in enumerate(regions)}

    if not significant.empty:
        colors = np.clip(
            significant["z_at_first_significant_scale"].to_numpy(float),
            -task["z_score_limit"],
            task["z_score_limit"],
        )
        radii = _radius(
            significant["first_significant_scale_um"].to_numpy(float),
            task["scale_levels"],
        )
        axis.scatter(
            significant["depth"],
            significant["analysis_region"].map(y_positions),
            c=colors,
            s=radii**2,
            cmap=cmap,
            norm=norm,
            edgecolors="#505050",
            linewidths=0.35,
            zorder=3,
        )

    depths = task["depths"]
    axis.set_xticks(depths)
    axis.set_xticklabels([f"{value:g}" for value in depths], rotation=45, ha="right")
    axis.set_yticks(range(len(regions)))
    axis.set_yticklabels([REGION_LABELS.get(x, x.replace("_", " ").title()) for x in regions])
    axis.set_xlim(min(depths) - 250, max(depths) + 250)
    axis.set_ylim(-0.6, len(regions) - 0.4)
    axis.invert_yaxis()
    axis.grid(color="#dddddd", linewidth=0.7, zorder=0)
    axis.set_axisbelow(True)
    axis.set_xlabel("AP depth (µm)")
    axis.set_ylabel("Region")
    fig.suptitle(
        f"Reference: {task['reference']}  →  Neighbor: {task['neighbor']}",
        fontweight="semibold",
        fontsize=13,
        y=0.97,
    )
    fig.text(
        0.5,
        0.915,
        (
            f"First significant CRAWDAD relationship; |mean Z| ≥ "
            f"{task['z_threshold']:g}; no point = no significant result (reason saved in point table)"
        ),
        ha="center",
        va="center",
        fontsize=8.5,
        color="#555555",
    )

    scalar = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    scalar.set_array([])
    colorbar = fig.colorbar(scalar, ax=axis, fraction=0.028, pad=0.025)
    colorbar.set_label("Mean Z at first significant scale")
    scale_breaks = task["scale_levels"]
    handles = [
        Line2D(
            [0], [0], marker="o", linestyle="", markerfacecolor="#bdbdbd",
            markeredgecolor="#505050", markeredgewidth=0.35,
            markersize=float(_radius(value, task["scale_levels"])),
            label=f"{value:g} µm",
        )
        for value in scale_breaks
    ]
    axis.legend(
        handles=handles,
        title=(
            f"First significant scale ({task['scale_interval_um']:g} µm levels)\n"
            "smaller scale = larger point"
        ),
        loc="upper left",
        bbox_to_anchor=(1.25, 0.66),
        frameon=False,
    )
    fig.subplots_adjust(right=0.70, bottom=0.2, top=0.84)
    _save_figure(fig, output, task["dpi"])
    return {
        "reference": task["reference"],
        "neighbor": task["neighbor"],
        "status": "written" if not significant.empty else "written_no_significant_points",
        "n_expected_positions": task["n_expected_positions"],
        "n_eligible_positions": int(pair["reference_eligible"].sum()),
        "n_significant_points": len(significant),
        "plot": str(output.resolve()),
    }


def _write_dataframe(data: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f".{path.name}.tmp.{os.getpid()}")
    if path.suffix == ".parquet":
        data.to_parquet(temporary, index=False)
    elif path.name.endswith(".csv.gz"):
        data.to_csv(temporary, index=False, compression="gzip")
    else:
        data.to_csv(temporary, index=False)
    temporary.replace(path)


def _jsonable_args(args: argparse.Namespace) -> dict[str, Any]:
    return {
        key: str(value.resolve()) if isinstance(value, Path) else value
        for key, value in vars(args).items()
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    data, input_qc, eligibility = load_inputs(args)
    references = select_levels(
        args.references, sorted(data["reference"].unique()), "--references"
    )
    neighbors = select_levels(
        args.neighbors, sorted(data["neighbor"].unique()), "--neighbors"
    )
    if not references or not neighbors:
        raise ValueError("No directional pairs selected")
    validate_safe_tokens(references, "Reference")
    validate_safe_tokens(neighbors, "Neighbor")

    depth_map = data[["sample", "depth"]].drop_duplicates().sort_values("depth")
    depths = depth_map["depth"].to_numpy(float)
    z_threshold = float(data["threshold_used"].dropna().iloc[0])
    observed_scales = data["first_significant_scale_um"].dropna().to_numpy(float)
    scale_levels = np.asarray(args.selected_scales_um, dtype=float)
    scale_min = float(scale_levels.min())
    scale_max = float(scale_levels.max())
    distances_to_levels = np.min(
        np.abs(observed_scales[:, np.newaxis] - scale_levels[np.newaxis, :]),
        axis=1,
    )
    if np.any(distances_to_levels > 1e-8):
        invalid = sorted(np.unique(observed_scales[distances_to_levels > 1e-8]))
        raise ValueError(
            "First-significant scales do not align to the requested "
            f"{args.scale_interval_um:g} µm interval: {invalid}"
        )

    expected = pd.MultiIndex.from_product(
        [args.regions, depth_map["sample"], references, neighbors],
        names=["analysis_region", "sample", "reference", "neighbor"],
    ).to_frame(index=False)
    expected = expected.merge(depth_map, on="sample", how="left", validate="many_to_one")
    selected_columns = [
        "analysis_region", "sample", "reference", "neighbor", "donor",
        "depth", "neighborhood_distance_um", "reference_eligible",
        "first_significant_scale_um", "z_at_first_significant_scale",
        "relationship_at_first_significant_scale", "ever_significant",
        "threshold_used",
    ]
    point_summary = expected.merge(
        eligibility,
        on=["analysis_region", "sample", "depth", "reference"],
        how="left",
        validate="many_to_one",
    )
    if point_summary["reference_eligible"].isna().any():
        raise ValueError("Reference eligibility is missing for expected region/sample/reference rows")
    point_summary = point_summary.merge(
        data[selected_columns],
        on=["analysis_region", "sample", "depth", "reference", "neighbor"],
        how="left",
        validate="one_to_one",
        suffixes=("", "_trend"),
    )
    eligibility_disagreement = point_summary["reference_eligible_trend"].notna() & (
        point_summary["reference_eligible_trend"]
        != point_summary["reference_eligible"]
    )
    if eligibility_disagreement.any():
        raise ValueError("Trend-summary and QC reference eligibility disagree")
    point_summary["ever_significant"] = point_summary["ever_significant"].eq(True)
    point_summary["pair_result_present"] = point_summary[
        "reference_eligible_trend"
    ].notna()
    point_summary["point_status"] = np.select(
        [
            ~point_summary["reference_eligible"],
            ~point_summary["pair_result_present"],
            point_summary["ever_significant"],
        ],
        ["reference_ineligible", "eligible_pair_absent", "significant_point"],
        default="eligible_not_significant",
    )
    point_summary = point_summary.sort_values(
        ["reference", "neighbor", "depth", "analysis_region"], kind="stable"
    )

    distance_token = f"neighdist_{args.neighborhood_distance_um:g}"
    output_root = args.output_dir / distance_token
    plot_root = args.plot_dir / distance_token
    output_root.mkdir(parents=True, exist_ok=True)
    plot_root.mkdir(parents=True, exist_ok=True)
    _write_dataframe(point_summary, output_root / "relationship_points.parquet")
    _write_dataframe(point_summary, output_root / "relationship_points.csv.gz")
    _write_dataframe(input_qc, output_root / "input_qc.csv")

    tasks: list[dict[str, Any]] = []
    for reference in references:
        for neighbor in neighbors:
            reference_token = safe_token(reference)
            neighbor_token = safe_token(neighbor)
            tasks.append(
                {
                    "data": point_summary.loc[
                        (point_summary["reference"] == reference)
                        & (point_summary["neighbor"] == neighbor)
                    ].copy(),
                    "reference": reference,
                    "neighbor": neighbor,
                    "regions": list(args.regions),
                    "depths": depths,
                    "n_expected_positions": len(args.regions) * len(depths),
                    "z_threshold": z_threshold,
                    "z_score_limit": args.z_score_limit,
                    "scale_min": scale_min,
                    "scale_max": scale_max,
                    "scale_levels": scale_levels,
                    "scale_interval_um": args.scale_interval_um,
                    "dpi": args.dpi,
                    "overwrite": args.overwrite,
                    "output": str(
                        plot_root / "by_reference" / reference_token
                        / f"{reference_token}__to__{neighbor_token}.png"
                    ),
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
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = [executor.submit(render_pair, task) for task in tasks]
            for index, future in enumerate(concurrent.futures.as_completed(futures), start=1):
                results.append(future.result())
                if index == 1 or index % 25 == 0 or index == len(tasks):
                    print(f"Completed {index}/{len(tasks)} pairs", flush=True)

    manifest = pd.DataFrame(results).sort_values(["reference", "neighbor"], kind="stable")
    _write_dataframe(manifest, output_root / "plot_manifest.csv")
    metadata = {
        "arguments": _jsonable_args(args),
        "input_definition": (
            "Module 02 all_trends_summary filtered to selected scales; "
            "first-significant relationships recomputed by Module 05"
        ),
        "point_definition": {
            "presence": "relationship first reaches abs(mean Z) >= threshold at any tested scale",
            "color": "mean Z at the first (smallest) significant scale",
            "size": "first significant scale, reverse mapped so smaller scale is a larger point",
            "no_point": "eligible but not significant, reference ineligible, or eligible pair absent; distinguish using relationship_points tables",
        },
        "resolved_references": references,
        "resolved_neighbors": neighbors,
        "regions": list(args.regions),
        "ap_depths_um": depths.tolist(),
        "z_threshold": z_threshold,
        "z_score_color_limit": args.z_score_limit,
        "first_significant_scale_range_um": [scale_min, scale_max],
        "first_significant_scale_interval_um": args.scale_interval_um,
        "first_significant_scale_levels_um": scale_levels.tolist(),
        "n_directional_pairs": len(tasks),
    }
    metadata_path = output_root / "resolved_arguments.json"
    temporary = metadata_path.with_name(f".{metadata_path.name}.tmp.{os.getpid()}")
    with temporary.open("w", encoding="utf-8") as handle:
        json.dump(metadata, handle, indent=2)
        handle.write("\n")
    temporary.replace(metadata_path)
    print(f"Plots: {plot_root.resolve()}", flush=True)
    print(f"Point table: {(output_root / 'relationship_points.parquet').resolve()}", flush=True)
    print(f"Manifest: {(output_root / 'plot_manifest.csv').resolve()}", flush=True)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        raise
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)
