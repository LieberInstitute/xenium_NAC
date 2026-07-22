#!/usr/bin/env python3
"""Compare unaligned, first-pass, and second-pass pair-level similarity metrics.

The default comparison uses the reliability-filtered multiscale median from
the Module 00 and two Module 03 pair-summary tables. Each metric remains
separate; no cross-metric composite score is calculated.
"""

from __future__ import annotations

import argparse
import json
import logging
import platform
import re
import sys
import tempfile
from pathlib import Path
from typing import Any, Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


LOGGER = logging.getLogger("module04_feature_similarity_comparison")
SCRIPT_VERSION = "2.0"
METRICS = (
    ("pcc_all", "Pearson correlation"),
    ("cos_sim_all", "Cosine similarity"),
    ("ssim_all", "SSIM"),
    ("mi_all", "Mutual information"),
)
SUMMARY_STATISTICS = ("multiscale_median", "gene_median")


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unaligned-pair-summary-csv", type=Path)
    parser.add_argument("--first-pass-pair-summary-csv", type=Path)
    parser.add_argument("--second-pass-pair-summary-csv", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--plot-dir", type=Path)
    parser.add_argument("--donor")
    parser.add_argument(
        "--summary-statistic", choices=SUMMARY_STATISTICS,
        default="multiscale_median",
        help=(
            "multiscale_median is reliability-filtered; gene_median is the "
            "unfiltered finite multiscale-gene median."
        ),
    )
    parser.add_argument("--plot-dpi", type=int, default=250)
    parser.add_argument("--overwrite", action="store_true")
    parser.add_argument("--run-self-tests", action="store_true")
    args = parser.parse_args(argv)
    required = (
        "unaligned_pair_summary_csv", "first_pass_pair_summary_csv",
        "second_pass_pair_summary_csv", "output_dir", "plot_dir", "donor",
    )
    if not args.run_self_tests and any(getattr(args, name) is None for name in required):
        parser.error(
            "--unaligned-pair-summary-csv, --first-pass-pair-summary-csv, "
            "--second-pass-pair-summary-csv, --output-dir, --plot-dir, and "
            "--donor are required unless --run-self-tests is used"
        )
    if args.plot_dpi < 1:
        parser.error("--plot-dpi must be positive")
    if not args.run_self_tests and args.output_dir.resolve() == args.plot_dir.resolve():
        parser.error("--output-dir and --plot-dir must differ")
    return args


def safe_filename(value: str) -> str:
    """Convert a value into a conservative filename component."""
    cleaned = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    return cleaned or "unnamed"


def display_sample_name(sample: str, donor: str) -> str:
    """Remove an exact donor prefix for plot display only."""
    prefix = f"{donor}_"
    return str(sample)[len(prefix):] if str(sample).startswith(prefix) else str(sample)


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    """Return all generated output paths."""
    donor = safe_filename(args.donor)
    statistic = safe_filename(args.summary_statistic)
    return {
        "pair_comparison": args.output_dir / f"{donor}_unaligned_first_second_{statistic}_by_pair.csv",
        "metric_summary": args.output_dir / f"{donor}_unaligned_first_second_{statistic}_metric_summary.csv",
        "metadata": args.output_dir / f"{donor}_unaligned_first_second_{statistic}_metadata.json",
        "plot": args.plot_dir / f"{donor}_unaligned_first_second_{statistic}.png",
        "final_pair_comparison": args.output_dir / f"{donor}_unaligned_vs_final_spateo_{statistic}_by_pair.csv",
        "final_metric_summary": args.output_dir / f"{donor}_unaligned_vs_final_spateo_{statistic}_metric_summary.csv",
        "final_plot": args.plot_dir / f"{donor}_unaligned_vs_final_spateo_{statistic}.png",
    }


def prepare_output_dirs(args: argparse.Namespace) -> dict[str, Path]:
    """Create output directories and enforce overwrite protection."""
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)
    paths = output_paths(args)
    conflicts = [path for path in paths.values() if path.exists()]
    if conflicts and not args.overwrite:
        raise FileExistsError(
            "Outputs already exist; use --overwrite: "
            + ", ".join(map(str, conflicts))
        )
    if args.overwrite:
        for path in conflicts:
            if path.is_file():
                path.unlink()
    return paths


def atomic_csv_dump(frame: pd.DataFrame, path: Path) -> None:
    """Atomically write a CSV table."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False)
    temporary.replace(path)


def atomic_json_dump(payload: dict[str, Any], path: Path) -> None:
    """Atomically write JSON metadata."""
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def validate_summary_table(
    frame: pd.DataFrame, table_name: str, statistic: str
) -> None:
    """Validate pair identity and selected metric columns."""
    required = {
        "pair_number", "pair_id", "reference_slice", "moving_slice",
        *(f"{metric}_{statistic}" for metric, _ in METRICS),
    }
    missing = required - set(frame.columns)
    if missing:
        raise ValueError(f"{table_name} is missing columns: {sorted(missing)}")
    if frame["pair_id"].isna().any() or frame["pair_id"].duplicated().any():
        raise ValueError(f"{table_name} pair_id values must be nonmissing and unique")


def build_pair_comparison(
    unaligned: pd.DataFrame,
    first_pass: pd.DataFrame,
    second_pass: pd.DataFrame,
    statistic: str,
) -> pd.DataFrame:
    """Build a long table strictly matched across all three alignment stages."""
    validate_summary_table(unaligned, "Unaligned summary", statistic)
    validate_summary_table(first_pass, "First-pass summary", statistic)
    validate_summary_table(second_pass, "Second-pass summary", statistic)
    id_sets = {
        "unaligned": set(unaligned["pair_id"]),
        "first-pass": set(first_pass["pair_id"]),
        "second-pass": set(second_pass["pair_id"]),
    }
    all_ids = set.union(*id_sets.values())
    if any(ids != all_ids for ids in id_sets.values()):
        raise ValueError(
            "Pair IDs differ between summaries; missing IDs by stage="
            + json.dumps(
                {
                    stage: sorted(all_ids - ids)
                    for stage, ids in id_sets.items()
                    if ids != all_ids
                },
                sort_keys=True,
            )
        )
    identity = ["pair_id", "reference_slice", "moving_slice"]
    metric_columns = [f"{metric}_{statistic}" for metric, _ in METRICS]

    def prepare(frame: pd.DataFrame, stage: str) -> pd.DataFrame:
        renamed = {"pair_number": f"pair_number_{stage}"}
        renamed.update({column: f"{column}_{stage}" for column in metric_columns})
        return frame[["pair_number", *identity, *metric_columns]].rename(
            columns=renamed
        )

    merged = (
        prepare(unaligned, "unaligned")
        .merge(
            prepare(first_pass, "first_pass"),
            on=identity,
            validate="one_to_one",
        )
        .merge(
            prepare(second_pass, "second_pass"),
            on=identity,
            validate="one_to_one",
        )
        .sort_values("pair_number_unaligned")
    )
    pair_numbers = [
        merged[f"pair_number_{stage}"].to_numpy()
        for stage in ("unaligned", "first_pass", "second_pass")
    ]
    if not all(np.array_equal(pair_numbers[0], values) for values in pair_numbers[1:]):
        raise ValueError("pair_number differs among the three summaries")

    rows: list[dict[str, Any]] = []
    for pair in merged.itertuples(index=False):
        for metric, metric_label in METRICS:
            column = f"{metric}_{statistic}"
            unaligned_score = float(getattr(pair, f"{column}_unaligned"))
            first_pass_score = float(getattr(pair, f"{column}_first_pass"))
            second_pass_score = float(getattr(pair, f"{column}_second_pass"))
            all_scores_finite = bool(np.isfinite([
                unaligned_score, first_pass_score, second_pass_score,
            ]).all())
            rows.append(
                {
                    "pair_number": int(pair.pair_number_unaligned),
                    "pair_id": pair.pair_id,
                    "reference_slice": pair.reference_slice,
                    "moving_slice": pair.moving_slice,
                    "metric": metric,
                    "metric_label": metric_label,
                    "summary_statistic": statistic,
                    "unaligned_score": unaligned_score,
                    "first_pass_score": first_pass_score,
                    "second_pass_score": second_pass_score,
                    "delta_first_pass_minus_unaligned": (
                        first_pass_score - unaligned_score
                        if all_scores_finite else np.nan
                    ),
                    "delta_second_pass_minus_first_pass": (
                        second_pass_score - first_pass_score
                        if all_scores_finite else np.nan
                    ),
                    "delta_second_pass_minus_unaligned": (
                        second_pass_score - unaligned_score
                        if all_scores_finite else np.nan
                    ),
                    "all_scores_finite": all_scores_finite,
                    "first_pass_higher": (
                        bool(first_pass_score > unaligned_score)
                        if all_scores_finite else pd.NA
                    ),
                    "second_pass_higher_than_first": (
                        bool(second_pass_score > first_pass_score)
                        if all_scores_finite else pd.NA
                    ),
                    "second_pass_higher_than_unaligned": (
                        bool(second_pass_score > unaligned_score)
                        if all_scores_finite else pd.NA
                    ),
                }
            )
    return pd.DataFrame.from_records(rows)


def summarize_across_pairs(comparison: pd.DataFrame) -> pd.DataFrame:
    """Summarize each metric across finite matched pairs without combining metrics."""
    rows: list[dict[str, Any]] = []
    for metric, metric_label in METRICS:
        group = comparison[
            comparison["metric"].eq(metric) & comparison["all_scores_finite"]
        ]
        rows.append(
            {
                "metric": metric,
                "metric_label": metric_label,
                "summary_statistic": comparison["summary_statistic"].iloc[0],
                "n_pairs_total": int(comparison["pair_id"].nunique()),
                "n_pairs_finite": len(group),
                "median_unaligned_score": (
                    float(group["unaligned_score"].median()) if len(group) else np.nan
                ),
                "median_first_pass_score": (
                    float(group["first_pass_score"].median()) if len(group) else np.nan
                ),
                "median_second_pass_score": (
                    float(group["second_pass_score"].median()) if len(group) else np.nan
                ),
                "median_delta_first_pass_minus_unaligned": (
                    float(group["delta_first_pass_minus_unaligned"].median())
                    if len(group) else np.nan
                ),
                "median_delta_second_pass_minus_first_pass": (
                    float(group["delta_second_pass_minus_first_pass"].median())
                    if len(group) else np.nan
                ),
                "median_delta_second_pass_minus_unaligned": (
                    float(group["delta_second_pass_minus_unaligned"].median())
                    if len(group) else np.nan
                ),
                "n_pairs_first_pass_higher": (
                    int(group["first_pass_higher"].sum()) if len(group) else 0
                ),
                "fraction_pairs_first_pass_higher": (
                    float(group["first_pass_higher"].mean()) if len(group) else np.nan
                ),
                "n_pairs_second_pass_higher_than_first": (
                    int(group["second_pass_higher_than_first"].sum())
                    if len(group) else 0
                ),
                "fraction_pairs_second_pass_higher_than_first": (
                    float(group["second_pass_higher_than_first"].mean())
                    if len(group) else np.nan
                ),
                "n_pairs_second_pass_higher_than_unaligned": (
                    int(group["second_pass_higher_than_unaligned"].sum())
                    if len(group) else 0
                ),
                "fraction_pairs_second_pass_higher_than_unaligned": (
                    float(group["second_pass_higher_than_unaligned"].mean())
                    if len(group) else np.nan
                ),
            }
        )
    return pd.DataFrame.from_records(rows)


def build_final_comparison(comparison: pd.DataFrame) -> pd.DataFrame:
    """Extract an independent unaligned-versus-final comparison table."""
    columns = [
        "pair_number", "pair_id", "reference_slice", "moving_slice",
        "metric", "metric_label", "summary_statistic", "unaligned_score",
        "second_pass_score",
    ]
    final = comparison[columns].copy().rename(
        columns={"second_pass_score": "final_spateo_score"}
    )
    finite = np.isfinite(final["unaligned_score"].to_numpy(float)) & np.isfinite(
        final["final_spateo_score"].to_numpy(float)
    )
    final["both_scores_finite"] = finite
    final["delta_final_spateo_minus_unaligned"] = np.where(
        finite,
        final["final_spateo_score"] - final["unaligned_score"],
        np.nan,
    )
    final["final_spateo_higher"] = pd.array(
        np.where(
            finite,
            final["final_spateo_score"] > final["unaligned_score"],
            pd.NA,
        ),
        dtype="boolean",
    )
    return final


def summarize_final_across_pairs(final: pd.DataFrame) -> pd.DataFrame:
    """Summarize unaligned versus final Spateo without using first-pass scores."""
    rows: list[dict[str, Any]] = []
    for metric, metric_label in METRICS:
        group = final[
            final["metric"].eq(metric) & final["both_scores_finite"]
        ]
        rows.append(
            {
                "metric": metric,
                "metric_label": metric_label,
                "summary_statistic": final["summary_statistic"].iloc[0],
                "n_pairs_total": int(final["pair_id"].nunique()),
                "n_pairs_finite": len(group),
                "median_unaligned_score": (
                    float(group["unaligned_score"].median())
                    if len(group) else np.nan
                ),
                "median_final_spateo_score": (
                    float(group["final_spateo_score"].median())
                    if len(group) else np.nan
                ),
                "median_delta_final_spateo_minus_unaligned": (
                    float(group["delta_final_spateo_minus_unaligned"].median())
                    if len(group) else np.nan
                ),
                "n_pairs_final_spateo_higher": (
                    int(group["final_spateo_higher"].sum()) if len(group) else 0
                ),
                "fraction_pairs_final_spateo_higher": (
                    float(group["final_spateo_higher"].mean())
                    if len(group) else np.nan
                ),
            }
        )
    return pd.DataFrame.from_records(rows)


def plot_pair_comparison(
    comparison: pd.DataFrame,
    output_path: Path,
    donor: str,
    statistic: str,
    dpi: int,
) -> None:
    """Plot three-stage scores in four separate metric panels."""
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), squeeze=False)
    for axis, (metric, metric_label) in zip(axes.ravel(), METRICS):
        group = comparison[comparison["metric"].eq(metric)].sort_values("pair_number")
        positions = np.arange(len(group))
        width = 0.26
        axis.bar(
            positions - width, group["unaligned_score"], width,
            label="Unaligned", color="#4C78A8",
        )
        axis.bar(
            positions, group["first_pass_score"], width,
            label="First-pass Spateo", color="#F58518",
        )
        axis.bar(
            positions + width, group["second_pass_score"], width,
            label="Second-pass Spateo", color="#54A24B",
        )
        labels = [
            f"{display_sample_name(reference, donor)}\n→\n"
            f"{display_sample_name(moving, donor)}"
            for reference, moving in zip(group["reference_slice"], group["moving_slice"])
        ]
        first_delta = group["delta_first_pass_minus_unaligned"].dropna()
        second_delta = group["delta_second_pass_minus_first_pass"].dropna()
        first_delta_text = (
            f"{first_delta.median():.3g}" if len(first_delta) else "NA"
        )
        second_delta_text = (
            f"{second_delta.median():.3g}" if len(second_delta) else "NA"
        )
        delta_text = (
            f"Median Δ first−unaligned={first_delta_text}\n"
            f"Median Δ second−first={second_delta_text}"
        )
        axis.text(
            0.02, 0.97, delta_text, transform=axis.transAxes,
            ha="left", va="top", fontsize=9,
        )
        axis.set_xticks(positions, labels, rotation=45, ha="right", fontsize=7)
        axis.set_ylabel(f"{metric}_{statistic}")
        axis.set_title(metric_label)
        axis.legend(loc="upper right", frameon=False, fontsize=8)
        axis.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    figure.suptitle(
        "Unaligned vs first-pass vs second-pass Spateo: "
        f"{statistic.replace('_', ' ')}",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def plot_final_comparison(
    final: pd.DataFrame,
    output_path: Path,
    donor: str,
    statistic: str,
    dpi: int,
) -> None:
    """Plot paired unaligned and final Spateo scores for each metric."""
    figure, axes = plt.subplots(2, 2, figsize=(16, 10), squeeze=False)
    for axis, (metric, metric_label) in zip(axes.ravel(), METRICS):
        group = final[final["metric"].eq(metric)].sort_values("pair_number")
        positions = np.arange(len(group))
        width = 0.38
        axis.bar(
            positions - width / 2, group["unaligned_score"], width,
            label="Unaligned", color="#4C78A8",
        )
        axis.bar(
            positions + width / 2, group["final_spateo_score"], width,
            label="Final Spateo alignment", color="#54A24B",
        )
        labels = [
            f"{display_sample_name(reference, donor)}\n→\n"
            f"{display_sample_name(moving, donor)}"
            for reference, moving in zip(
                group["reference_slice"], group["moving_slice"]
            )
        ]
        finite_delta = group["delta_final_spateo_minus_unaligned"].dropna()
        delta_text = (
            f"Median Δ final−unaligned={finite_delta.median():.3g}"
            if len(finite_delta) else "Median Δ final−unaligned=NA"
        )
        axis.text(
            0.02, 0.97, delta_text, transform=axis.transAxes,
            ha="left", va="top", fontsize=9,
        )
        axis.set_xticks(positions, labels, rotation=45, ha="right", fontsize=7)
        axis.set_ylabel(f"{metric}_{statistic}")
        axis.set_title(metric_label)
        axis.legend(loc="upper right", frameon=False, fontsize=8)
        axis.axhline(0, color="black", linewidth=0.6, alpha=0.4)
    figure.suptitle(
        "Unaligned vs final Spateo alignment: "
        f"{statistic.replace('_', ' ')}",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)


def run_self_tests() -> None:
    """Test strict pair matching, deltas, aggregation, and plotting."""
    rows = []
    for pair_number, (reference, moving) in enumerate(
        (("D1_S1", "D1_S2"), ("D1_S2", "D1_S3")), 1
    ):
        row: dict[str, Any] = {
            "pair_number": pair_number,
            "pair_id": f"pair_{pair_number:02d}_{reference}__{moving}",
            "reference_slice": reference,
            "moving_slice": moving,
        }
        for metric_index, (metric, _) in enumerate(METRICS, 1):
            row[f"{metric}_multiscale_median"] = 0.1 * metric_index + pair_number
            row[f"{metric}_gene_median"] = 0.2 * metric_index + pair_number
        rows.append(row)
    unaligned = pd.DataFrame.from_records(rows)
    first_pass = unaligned.copy()
    second_pass = unaligned.copy()
    for metric, _ in METRICS:
        first_pass[f"{metric}_multiscale_median"] += 0.25
        first_pass[f"{metric}_gene_median"] += 0.10
        second_pass[f"{metric}_multiscale_median"] += 0.50
        second_pass[f"{metric}_gene_median"] += 0.20
    first_pass = first_pass.iloc[::-1].reset_index(drop=True)
    second_pass = second_pass.iloc[::-1].reset_index(drop=True)
    comparison = build_pair_comparison(
        unaligned, first_pass, second_pass, "multiscale_median"
    )
    assert len(comparison) == 2 * len(METRICS)
    assert np.allclose(comparison["delta_first_pass_minus_unaligned"], 0.25)
    assert np.allclose(comparison["delta_second_pass_minus_first_pass"], 0.25)
    assert np.allclose(comparison["delta_second_pass_minus_unaligned"], 0.50)
    assert comparison["pair_id"].drop_duplicates().tolist() == unaligned["pair_id"].tolist()
    metric_summary = summarize_across_pairs(comparison)
    assert len(metric_summary) == len(METRICS)
    assert np.allclose(
        metric_summary["median_delta_first_pass_minus_unaligned"], 0.25
    )
    assert np.allclose(
        metric_summary["median_delta_second_pass_minus_first_pass"], 0.25
    )
    final_input = comparison.copy()
    final_input.loc[0, "first_pass_score"] = np.nan
    final_comparison = build_final_comparison(final_input)
    assert final_comparison.loc[0, "both_scores_finite"]
    assert np.allclose(
        final_comparison["delta_final_spateo_minus_unaligned"], 0.50
    )
    final_summary = summarize_final_across_pairs(final_comparison)
    assert len(final_summary) == len(METRICS)
    assert np.allclose(
        final_summary["median_delta_final_spateo_minus_unaligned"], 0.50
    )
    with tempfile.TemporaryDirectory() as temporary_directory:
        plot_path = Path(temporary_directory) / "comparison.png"
        final_plot_path = Path(temporary_directory) / "final_comparison.png"
        plot_pair_comparison(
            comparison, plot_path, "D1", "multiscale_median", dpi=80
        )
        plot_final_comparison(
            final_comparison, final_plot_path, "D1", "multiscale_median",
            dpi=80,
        )
        assert plot_path.is_file() and final_plot_path.is_file()
    LOGGER.info("All self-tests passed")


def main(argv: Sequence[str] | None = None) -> int:
    """Run the unaligned/first-pass/second-pass comparison."""
    args = parse_args(argv)
    logging.basicConfig(
        level=logging.INFO, format="%(asctime)s | %(levelname)s | %(message)s"
    )
    if args.run_self_tests:
        run_self_tests()
        return 0
    for path, label in (
        (args.unaligned_pair_summary_csv, "Unaligned pair summary"),
        (args.first_pass_pair_summary_csv, "First-pass pair summary"),
        (args.second_pass_pair_summary_csv, "Second-pass pair summary"),
    ):
        if not path.is_file():
            raise FileNotFoundError(f"{label} does not exist: {path}")
    paths = prepare_output_dirs(args)
    unaligned = pd.read_csv(args.unaligned_pair_summary_csv)
    first_pass = pd.read_csv(args.first_pass_pair_summary_csv)
    second_pass = pd.read_csv(args.second_pass_pair_summary_csv)
    comparison = build_pair_comparison(
        unaligned, first_pass, second_pass, args.summary_statistic
    )
    metric_summary = summarize_across_pairs(comparison)
    final_comparison = build_final_comparison(comparison)
    final_metric_summary = summarize_final_across_pairs(final_comparison)
    atomic_csv_dump(comparison, paths["pair_comparison"])
    atomic_csv_dump(metric_summary, paths["metric_summary"])
    atomic_csv_dump(final_comparison, paths["final_pair_comparison"])
    atomic_csv_dump(final_metric_summary, paths["final_metric_summary"])
    plot_pair_comparison(
        comparison, paths["plot"], args.donor,
        args.summary_statistic, args.plot_dpi,
    )
    plot_final_comparison(
        final_comparison, paths["final_plot"], args.donor,
        args.summary_statistic, args.plot_dpi,
    )
    metadata = {
        "script_version": SCRIPT_VERSION,
        "unaligned_pair_summary_csv": str(args.unaligned_pair_summary_csv.resolve()),
        "first_pass_pair_summary_csv": str(args.first_pass_pair_summary_csv.resolve()),
        "second_pass_pair_summary_csv": str(args.second_pass_pair_summary_csv.resolve()),
        "donor": str(args.donor),
        "summary_statistic": args.summary_statistic,
        "metrics": [metric for metric, _ in METRICS],
        "cross_metric_composite_score": False,
        "delta_definitions": [
            "first_pass_score - unaligned_score",
            "second_pass_score - first_pass_score",
            "second_pass_score - unaligned_score",
        ],
        "n_pairs": int(comparison["pair_id"].nunique()),
        "outputs": {name: str(path) for name, path in paths.items()},
        "software_versions": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "matplotlib": matplotlib.__version__,
        },
    }
    atomic_json_dump(metadata, paths["metadata"])
    LOGGER.info("Comparison complete. Plot: %s", paths["plot"])
    return 0


if __name__ == "__main__":
    sys.exit(main())


# Example usage:
# python 04_feature_similarity_comparison.py \
#     --unaligned-pair-summary-csv /path/to/Br6660_unaligned_feature_similarity_pair_summary.csv \
#     --first-pass-pair-summary-csv /path/to/Br6660_first_pass_feature_similarity_pair_summary.csv \
#     --second-pass-pair-summary-csv /path/to/Br6660_second_pass_feature_similarity_pair_summary.csv \
#     --output-dir /path/to/module04_feature_similarity_comparison_Br6660 \
#     --plot-dir /path/to/module04_feature_similarity_comparison_Br6660 \
#     --donor Br6660 \
#     --summary-statistic multiscale_median \
#     --plot-dpi 250 \
#     --overwrite
