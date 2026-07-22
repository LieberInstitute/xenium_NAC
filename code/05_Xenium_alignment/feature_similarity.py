"""Reusable grid-based gene-similarity scores for Xenium alignment evaluation.

Coordinates are interpreted in micrometers. Both slices must already be placed
in the candidate alignment being evaluated; this module evaluates alignment and
does not perform alignment. Partial tissue overlap is allowed, and overlap area
is treated as a reliability diagnostic rather than an optimization objective.

The caller controls the gene list and expression preprocessing and should reuse
the same expression data and genes for every alignment candidate. The default
Xenium scales are 200, 300, and 400 µm to match Module 1. Results are retained
per scale and also summarized across scales. Multiple grid offsets reduce
sensitivity to arbitrary square-grid boundaries, but do not make the score
rotation invariant.

This module contains no pipeline file I/O, directory searching, Spateo calls,
candidate selection, or slice-stack orchestration.
"""

from __future__ import annotations

import warnings
from collections.abc import Mapping, Sequence
from typing import Any, Optional

import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.stats import pearsonr
from sklearn.feature_selection import mutual_info_regression


DEFAULT_GRID_SIZES_UM = (200.0, 300.0, 400.0)
SUPPORTED_METRICS = ("pcc", "cosine", "ssim", "mi")
METRIC_PREFIX = {
    "pcc": "pcc",
    "cosine": "cos_sim",
    "ssim": "ssim",
    "mi": "mi",
}

# Retained only for backward compatibility with the legacy cross-technology API.
GRID_SIZE_PX = {
    ("visium", "visium"): 500,
    ("visium", "visium_hd"): 500,
    ("visium", "xenium"): 500,
    ("visium_hd", "visium_hd"): 300,
    ("visium_hd", "xenium"): 300,
    ("xenium", "xenium"): 250,
    ("merfish", "merfish"): 250,
    ("bins", "bins"): 100,
    ("cell_obs_sectioned", "cell_obs_sectioned"): 300,
    ("spots", "spots"): 200,
    ("bins", "cell_obs_sectioned"): 300,
    ("bins", "spots"): 200,
    ("cell_obs_sectioned", "spots"): 300,
}


# ---------------------------------------------------------------------------
# Metric and expression helpers
# ---------------------------------------------------------------------------

def _pair_key(tech1: str, tech2: str) -> tuple[str, str]:
    """Return the canonical legacy technology-pair key."""
    return tuple(sorted((tech1, tech2)))  # type: ignore[return-value]


def _dense_column(matrix: Any, column_index: int) -> np.ndarray:
    """Extract one sparse or dense expression column as a float vector."""
    column = matrix[:, column_index]
    if sp.issparse(column):
        column = column.toarray()
    return np.asarray(column, dtype=np.float64).ravel()


def _grid_means(expression: np.ndarray, labels: np.ndarray, n_grids: int) -> np.ndarray:
    """Calculate the mean expression in each nonnegative grid label."""
    output = np.zeros(n_grids, dtype=np.float64)
    counts = np.zeros(n_grids, dtype=np.int64)
    mask = labels >= 0
    np.add.at(output, labels[mask], expression[mask])
    np.add.at(counts, labels[mask], 1)
    with np.errstate(divide="ignore", invalid="ignore"):
        output = np.where(counts > 0, output / counts, 0.0)
    return output


def _pcc(a: np.ndarray, b: np.ndarray) -> float:
    """Return Pearson correlation, or NaN when mathematically undefined."""
    if a.size < 2 or a.std() == 0 or b.std() == 0:
        return np.nan
    correlation, _ = pearsonr(a, b)
    return float(correlation)


def _cos_sim(a: np.ndarray, b: np.ndarray) -> float:
    """Return cosine similarity, or NaN for zero/insufficient vectors."""
    norm_a, norm_b = np.linalg.norm(a), np.linalg.norm(b)
    if a.size < 2 or norm_a == 0 or norm_b == 0:
        return np.nan
    return float(a @ b / (norm_a * norm_b))


def _ssim(a: np.ndarray, b: np.ndarray, dynamic_range: float = 1.0) -> float:
    """Return 1-D SSIM after independently max-normalizing both vectors."""
    if a.size < 2:
        return np.nan
    max_a, max_b = a.max(), b.max()
    if max_a == 0 or max_b == 0:
        return np.nan
    a = a / max_a
    b = b / max_b
    mean_a, mean_b = a.mean(), b.mean()
    sd_a = np.sqrt(((a - mean_a) ** 2).mean())
    sd_b = np.sqrt(((b - mean_b) ** 2).mean())
    covariance = ((a - mean_a) * (b - mean_b)).mean()
    c1 = (0.01 * dynamic_range) ** 2
    c2 = (0.03 * dynamic_range) ** 2
    c3 = c2 / 2
    luminance = (2 * mean_a * mean_b + c1) / (mean_a**2 + mean_b**2 + c1)
    contrast = (2 * sd_a * sd_b + c2) / (sd_a**2 + sd_b**2 + c2)
    denominator = sd_a * sd_b + c3
    structure = (covariance + c3) / denominator if denominator > 0 else np.nan
    return float(luminance * contrast * structure)


def _mi(a: np.ndarray, b: np.ndarray, random_state: int = 0) -> float:
    """Return mutual information, or NaN when estimation is undefined."""
    if a.size < 4 or a.std() == 0 or b.std() == 0:
        return np.nan
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            value = mutual_info_regression(
                a.reshape(-1, 1), b, random_state=random_state
            )
        except ValueError:
            return np.nan
    return float(value[0])


def _metric_value(
    metric: str, a: np.ndarray, b: np.ndarray, random_state: int
) -> float:
    if metric == "pcc":
        return _pcc(a, b)
    if metric == "cosine":
        return _cos_sim(a, b)
    if metric == "ssim":
        return _ssim(a, b)
    if metric == "mi":
        return _mi(a, b, random_state=random_state)
    raise ValueError(f"Unsupported metric: {metric!r}")


def _finite_summary(values: Sequence[float]) -> tuple[float, float, float, float]:
    """Return median, minimum, maximum, and population SD over finite values."""
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return np.nan, np.nan, np.nan, np.nan
    return (
        float(np.median(array)),
        float(np.min(array)),
        float(np.max(array)),
        float(np.std(array, ddof=0)),
    )


def _validate_coordinates(
    coordinates: Any, n_observations: int, name: str
) -> np.ndarray:
    array = np.asarray(coordinates, dtype=np.float64)
    if array.shape != (n_observations, 2):
        raise ValueError(
            f"{name} must have shape ({n_observations}, 2); got {array.shape}"
        )
    if not np.isfinite(array).all():
        raise ValueError(f"{name} contains nonfinite coordinates")
    return array.copy()


def _resolve_coordinates(
    adata: Any,
    explicit_coordinates: Any,
    coordinate_key: str,
    name: str,
) -> tuple[np.ndarray, str]:
    if explicit_coordinates is not None:
        return _validate_coordinates(explicit_coordinates, adata.n_obs, name), "explicit"
    if coordinate_key not in adata.obsm:
        raise KeyError(f"{name} was not supplied and adata.obsm[{coordinate_key!r}] is absent")
    return (
        _validate_coordinates(adata.obsm[coordinate_key], adata.n_obs, name),
        f"obsm:{coordinate_key}",
    )


def _resolve_expression_matrix(adata: Any, expression_layer: Optional[str]) -> Any:
    if expression_layer is None:
        return adata.X
    if expression_layer not in adata.layers:
        raise KeyError(f"AnnData layer {expression_layer!r} is absent")
    return adata.layers[expression_layer]


def _resolve_genes(
    reference_adata: Any,
    moving_adata: Any,
    gene_list: Optional[Sequence[str]],
) -> tuple[list[str], list[str], list[str]]:
    reference_genes = set(map(str, reference_adata.var_names))
    moving_genes = set(map(str, moving_adata.var_names))
    if gene_list is None:
        requested = [
            str(gene) for gene in reference_adata.var_names
            if str(gene) in moving_genes
        ]
    else:
        requested = [str(gene) for gene in gene_list]
        if len(requested) != len(set(requested)):
            raise ValueError("gene_list contains duplicate genes")
    scored = [
        gene for gene in requested
        if gene in reference_genes and gene in moving_genes
    ]
    missing = [gene for gene in requested if gene not in reference_genes or gene not in moving_genes]
    if missing:
        preview = ", ".join(missing[:10])
        warnings.warn(
            f"{len(missing)} requested genes are absent from one or both AnnData "
            f"objects and will not be scored. Examples: {preview}",
            UserWarning,
            stacklevel=2,
        )
    return requested, scored, missing


def _validate_grid_sizes(grid_sizes_um: Sequence[float]) -> tuple[float, ...]:
    values = tuple(float(value) for value in grid_sizes_um)
    if not values:
        raise ValueError("grid_sizes_um must contain at least one scale")
    if any(not np.isfinite(value) or value <= 0 for value in values):
        raise ValueError("Every grid size must be finite and greater than zero")
    if len(values) != len(set(values)):
        raise ValueError("grid_sizes_um contains duplicate scales")
    return values


def _validate_metrics(metrics: Sequence[str]) -> tuple[str, ...]:
    values = tuple(str(metric).lower() for metric in metrics)
    if not values:
        raise ValueError("metrics must contain at least one metric")
    invalid = [metric for metric in values if metric not in SUPPORTED_METRICS]
    if invalid:
        raise ValueError(
            f"Unsupported metrics {invalid}; choose from {SUPPORTED_METRICS}"
        )
    if len(values) != len(set(values)):
        raise ValueError("metrics contains duplicates")
    return values


def _validate_offset(offset: Sequence[float], grid_size_um: float) -> tuple[float, float]:
    if len(offset) != 2:
        raise ValueError(f"Each grid offset must contain two values; got {offset}")
    x, y = float(offset[0]), float(offset[1])
    if not np.isfinite([x, y]).all():
        raise ValueError(f"Grid offset must be finite; got {offset}")
    return x % grid_size_um, y % grid_size_um


def _offsets_for_scale(
    grid_offsets_um: Optional[
        Sequence[Sequence[float]] | Mapping[float, Sequence[Sequence[float]]]
    ],
    grid_size_um: float,
) -> list[tuple[float, float]]:
    if grid_offsets_um is None:
        half = grid_size_um / 2.0
        raw_offsets: Sequence[Sequence[float]] = (
            (0.0, 0.0), (half, 0.0), (0.0, half), (half, half)
        )
    elif isinstance(grid_offsets_um, Mapping):
        matching_key = next(
            (
                key for key in grid_offsets_um
                if np.isclose(float(key), grid_size_um)
            ),
            None,
        )
        if matching_key is None:
            raise KeyError(f"No grid offsets supplied for scale {grid_size_um:g} µm")
        raw_offsets = grid_offsets_um[matching_key]
    else:
        raw_offsets = grid_offsets_um
    offsets = [_validate_offset(offset, grid_size_um) for offset in raw_offsets]
    if not offsets:
        raise ValueError(f"No grid offsets supplied for scale {grid_size_um:g} µm")
    stable_offsets = [tuple(np.round(offset, decimals=12)) for offset in offsets]
    if len(stable_offsets) != len(set(stable_offsets)):
        raise ValueError(
            "Duplicate grid offsets remain after modulo normalization "
            f"for scale {grid_size_um:g} µm"
        )
    return offsets


# ---------------------------------------------------------------------------
# Grid construction
# ---------------------------------------------------------------------------

def assign_grid_labels(
    reference_coords: Any,
    moving_coords: Any,
    *,
    grid_size_um: float,
    offset_um: Sequence[float] = (0.0, 0.0),
    min_cells_reference: int = 5,
    min_cells_moving: int = 5,
) -> tuple[np.ndarray, np.ndarray, int, np.ndarray, np.ndarray, dict[str, Any]]:
    """Label mutually supported grids in the coordinate intersection.

    Returned count arrays contain counts for valid grids only and correspond to
    the renumbered labels ``0..n_valid_grids-1``.
    """
    reference = np.asarray(reference_coords, dtype=np.float64)
    moving = np.asarray(moving_coords, dtype=np.float64)
    if reference.ndim != 2 or reference.shape[1] != 2 or reference.shape[0] == 0:
        raise ValueError(f"reference_coords must have shape (n_cells, 2); got {reference.shape}")
    if moving.ndim != 2 or moving.shape[1] != 2 or moving.shape[0] == 0:
        raise ValueError(f"moving_coords must have shape (n_cells, 2); got {moving.shape}")
    if not np.isfinite(reference).all() or not np.isfinite(moving).all():
        raise ValueError("Coordinate arrays must contain only finite values")
    grid_size_um = float(grid_size_um)
    if not np.isfinite(grid_size_um) or grid_size_um <= 0:
        raise ValueError("grid_size_um must be finite and greater than zero")
    if min_cells_reference < 1 or min_cells_moving < 1:
        raise ValueError("Minimum cell thresholds must be positive integers")
    offset_x, offset_y = _validate_offset(offset_um, grid_size_um)

    x_min = max(reference[:, 0].min(), moving[:, 0].min())
    x_max = min(reference[:, 0].max(), moving[:, 0].max())
    y_min = max(reference[:, 1].min(), moving[:, 1].min())
    y_max = min(reference[:, 1].max(), moving[:, 1].max())
    base_metadata = {
        "grid_size_um": grid_size_um,
        "offset_x_um": offset_x,
        "offset_y_um": offset_y,
        "x_min": float(x_min),
        "x_max": float(x_max),
        "y_min": float(y_min),
        "y_max": float(y_max),
        "min_cells_reference": int(min_cells_reference),
        "min_cells_moving": int(min_cells_moving),
    }
    if x_min >= x_max or y_min >= y_max:
        metadata = {
            **base_metadata,
            "grid_origin_x_um": np.nan,
            "grid_origin_y_um": np.nan,
            "n_rows": 0,
            "n_cols": 0,
            "n_valid_grids": 0,
            "intersection_valid": False,
            "reference_mean_occupancy": 0.0,
            "moving_mean_occupancy": 0.0,
            "reference_cells_in_valid_grids": 0,
            "moving_cells_in_valid_grids": 0,
            "reference_mean_cells_per_valid_grid": np.nan,
            "moving_mean_cells_per_valid_grid": np.nan,
            "reference_median_cells_per_valid_grid": np.nan,
            "moving_median_cells_per_valid_grid": np.nan,
            "reference_min_cells_per_valid_grid": np.nan,
            "moving_min_cells_per_valid_grid": np.nan,
            "reference_max_cells_per_valid_grid": np.nan,
            "moving_max_cells_per_valid_grid": np.nan,
        }
        return (
            np.full(reference.shape[0], -1, dtype=np.int64),
            np.full(moving.shape[0], -1, dtype=np.int64),
            0,
            np.empty(0, dtype=np.int64),
            np.empty(0, dtype=np.int64),
            metadata,
        )

    origin_x = x_min - offset_x
    origin_y = y_min - offset_y
    n_cols = max(1, int(np.ceil((x_max - origin_x) / grid_size_um)))
    n_rows = max(1, int(np.ceil((y_max - origin_y) / grid_size_um)))
    n_grid_cells = n_rows * n_cols

    def flat_indices(coordinates: np.ndarray) -> np.ndarray:
        inside = (
            (coordinates[:, 0] >= x_min) & (coordinates[:, 0] <= x_max)
            & (coordinates[:, 1] >= y_min) & (coordinates[:, 1] <= y_max)
        )
        x_index = np.floor((coordinates[:, 0] - origin_x) / grid_size_um).astype(np.int64)
        y_index = np.floor((coordinates[:, 1] - origin_y) / grid_size_um).astype(np.int64)
        x_index = np.minimum(x_index, n_cols - 1)
        y_index = np.minimum(y_index, n_rows - 1)
        valid_index = (
            inside & (x_index >= 0) & (x_index < n_cols)
            & (y_index >= 0) & (y_index < n_rows)
        )
        return np.where(valid_index, y_index * n_cols + x_index, -1)

    reference_flat = flat_indices(reference)
    moving_flat = flat_indices(moving)
    reference_counts_all = np.bincount(
        reference_flat[reference_flat >= 0], minlength=n_grid_cells
    )
    moving_counts_all = np.bincount(
        moving_flat[moving_flat >= 0], minlength=n_grid_cells
    )
    valid_mask = (
        (reference_counts_all >= min_cells_reference)
        & (moving_counts_all >= min_cells_moving)
    )
    n_valid_grids = int(valid_mask.sum())
    flat_to_label = np.full(n_grid_cells, -1, dtype=np.int64)
    flat_to_label[valid_mask] = np.arange(n_valid_grids, dtype=np.int64)

    reference_labels = np.full(reference.shape[0], -1, dtype=np.int64)
    moving_labels = np.full(moving.shape[0], -1, dtype=np.int64)
    reference_inside = reference_flat >= 0
    moving_inside = moving_flat >= 0
    reference_labels[reference_inside] = flat_to_label[reference_flat[reference_inside]]
    moving_labels[moving_inside] = flat_to_label[moving_flat[moving_inside]]
    reference_valid_counts = reference_counts_all[valid_mask]
    moving_valid_counts = moving_counts_all[valid_mask]

    def occupancy_summary(counts: np.ndarray) -> tuple[int, float, float, float, float]:
        if counts.size == 0:
            return 0, np.nan, np.nan, np.nan, np.nan
        return (
            int(counts.sum()), float(counts.mean()), float(np.median(counts)),
            float(counts.min()), float(counts.max()),
        )

    reference_occupancy = occupancy_summary(reference_valid_counts)
    moving_occupancy = occupancy_summary(moving_valid_counts)
    metadata = {
        **base_metadata,
        "grid_origin_x_um": float(origin_x),
        "grid_origin_y_um": float(origin_y),
        "n_rows": n_rows,
        "n_cols": n_cols,
        "n_valid_grids": n_valid_grids,
        "intersection_valid": True,
        "reference_mean_occupancy": float(reference_counts_all.mean()),
        "moving_mean_occupancy": float(moving_counts_all.mean()),
        "reference_cells_in_valid_grids": reference_occupancy[0],
        "moving_cells_in_valid_grids": moving_occupancy[0],
        "reference_mean_cells_per_valid_grid": reference_occupancy[1],
        "moving_mean_cells_per_valid_grid": moving_occupancy[1],
        "reference_median_cells_per_valid_grid": reference_occupancy[2],
        "moving_median_cells_per_valid_grid": moving_occupancy[2],
        "reference_min_cells_per_valid_grid": reference_occupancy[3],
        "moving_min_cells_per_valid_grid": moving_occupancy[3],
        "reference_max_cells_per_valid_grid": reference_occupancy[4],
        "moving_max_cells_per_valid_grid": moving_occupancy[4],
    }
    return (
        reference_labels,
        moving_labels,
        n_valid_grids,
        reference_valid_counts,
        moving_valid_counts,
        metadata,
    )


# ---------------------------------------------------------------------------
# Xenium pairwise API
# ---------------------------------------------------------------------------

def compute_pair_gene_similarity(
    reference_adata: Any,
    moving_adata: Any,
    *,
    reference_coords: Any = None,
    moving_coords: Any = None,
    reference_coord_key: str = "spatial",
    moving_coord_key: str = "spatial",
    gene_list: Optional[Sequence[str]] = None,
    expression_layer: Optional[str] = None,
    grid_sizes_um: Sequence[float] = DEFAULT_GRID_SIZES_UM,
    grid_offsets_um: Optional[
        Sequence[Sequence[float]] | Mapping[float, Sequence[Sequence[float]]]
    ] = None,
    min_cells_reference: int = 5,
    min_cells_moving: int = 5,
    min_overlap_grids: int = 10,
    metrics: Sequence[str] = SUPPORTED_METRICS,
    random_state: int = 0,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Calculate multiscale per-gene similarity for one ordered slice pair.

    Returns ``(multiscale_gene_result, scale_gene_result, offset_result)``.
    Explicit coordinate arrays take precedence over ``obsm`` coordinates.
    """
    if min_cells_reference < 1 or min_cells_moving < 1:
        raise ValueError("Minimum cell thresholds must be positive integers")
    if min_overlap_grids < 1:
        raise ValueError("min_overlap_grids must be a positive integer")
    scales = _validate_grid_sizes(grid_sizes_um)
    metric_names = _validate_metrics(metrics)
    reference_coordinates, reference_coordinate_source = _resolve_coordinates(
        reference_adata, reference_coords, reference_coord_key, "reference_coords"
    )
    moving_coordinates, moving_coordinate_source = _resolve_coordinates(
        moving_adata, moving_coords, moving_coord_key, "moving_coords"
    )
    requested_genes, scored_genes, missing_genes = _resolve_genes(
        reference_adata, moving_adata, gene_list
    )
    if not scored_genes:
        raise ValueError(
            "No requested genes are present in both the reference and moving "
            "AnnData objects."
        )
    reference_matrix = _resolve_expression_matrix(reference_adata, expression_layer)
    moving_matrix = _resolve_expression_matrix(moving_adata, expression_layer)
    if reference_matrix.shape != (reference_adata.n_obs, reference_adata.n_vars):
        raise ValueError("Reference expression matrix has an unexpected shape")
    if moving_matrix.shape != (moving_adata.n_obs, moving_adata.n_vars):
        raise ValueError("Moving expression matrix has an unexpected shape")
    reference_index = {str(gene): index for index, gene in enumerate(reference_adata.var_names)}
    moving_index = {str(gene): index for index, gene in enumerate(moving_adata.var_names)}

    metric_columns = [
        f"{METRIC_PREFIX[metric]}_{variant}"
        for metric in metric_names for variant in ("all", "nz")
    ]
    grid_metadata_columns = [
        "x_min", "x_max", "y_min", "y_max", "grid_origin_x_um",
        "grid_origin_y_um", "n_rows", "n_cols", "n_valid_grids",
        "reference_cells_in_valid_grids", "moving_cells_in_valid_grids",
        "reference_mean_cells_per_valid_grid",
        "moving_mean_cells_per_valid_grid",
        "reference_median_cells_per_valid_grid",
        "moving_median_cells_per_valid_grid",
        "reference_min_cells_per_valid_grid",
        "moving_min_cells_per_valid_grid",
        "reference_max_cells_per_valid_grid",
        "moving_max_cells_per_valid_grid",
    ]
    grid_specs: list[dict[str, Any]] = []
    offsets_by_scale: dict[float, list[tuple[float, float]]] = {}
    for scale in scales:
        offsets = _offsets_for_scale(grid_offsets_um, scale)
        offsets_by_scale[scale] = offsets
        for offset_x, offset_y in offsets:
            (
                reference_labels,
                moving_labels,
                n_valid_grids,
                _,
                _,
                grid_metadata,
            ) = assign_grid_labels(
                reference_coordinates,
                moving_coordinates,
                grid_size_um=scale,
                offset_um=(offset_x, offset_y),
                min_cells_reference=min_cells_reference,
                min_cells_moving=min_cells_moving,
            )
            grid_specs.append(
                {
                    "grid_size_um": scale,
                    "offset_x_um": offset_x,
                    "offset_y_um": offset_y,
                    "reference_labels": reference_labels,
                    "moving_labels": moving_labels,
                    "n_valid_grids": n_valid_grids,
                    "support_sufficient": n_valid_grids >= min_overlap_grids,
                    "grid_metadata": grid_metadata,
                }
            )

    offset_rows: list[dict[str, Any]] = []
    for gene in scored_genes:
        reference_expression = _dense_column(
            reference_matrix, reference_index[gene]
        )
        moving_expression = _dense_column(
            moving_matrix, moving_index[gene]
        )
        for spec in grid_specs:
            n_valid_grids = spec["n_valid_grids"]
            grid_metadata = spec["grid_metadata"]
            row: dict[str, Any] = {
                "gene": gene,
                "grid_size_um": spec["grid_size_um"],
                "offset_x_um": spec["offset_x_um"],
                "offset_y_um": spec["offset_y_um"],
                "support_sufficient": spec["support_sufficient"],
                "n_grids_all": n_valid_grids,
                "n_grids_nz": 0,
                "intersection_valid": grid_metadata["intersection_valid"],
            }
            row.update(
                {column: grid_metadata[column] for column in grid_metadata_columns}
            )
            row.update({column: np.nan for column in metric_columns})
            if spec["support_sufficient"]:
                reference_means = _grid_means(
                    reference_expression, spec["reference_labels"], n_valid_grids
                )
                moving_means = _grid_means(
                    moving_expression, spec["moving_labels"], n_valid_grids
                )
                nonzero = (reference_means > 0) & (moving_means > 0)
                row["n_grids_nz"] = int(nonzero.sum())
                for metric in metric_names:
                    prefix = METRIC_PREFIX[metric]
                    row[f"{prefix}_all"] = _metric_value(
                        metric, reference_means, moving_means, random_state
                    )
                    row[f"{prefix}_nz"] = _metric_value(
                        metric,
                        reference_means[nonzero],
                        moving_means[nonzero],
                        random_state,
                    )
            offset_rows.append(row)

    offset_columns = [
        "gene", "grid_size_um", "offset_x_um", "offset_y_um",
        "support_sufficient", "intersection_valid", "n_grids_all", "n_grids_nz",
        *grid_metadata_columns,
        *metric_columns,
    ]
    offset_result = pd.DataFrame.from_records(offset_rows, columns=offset_columns)

    scale_rows: list[dict[str, Any]] = []
    for scale in scales:
        expected_offsets = len(offsets_by_scale[scale])
        for gene in scored_genes:
            group = offset_result[
                (offset_result["grid_size_um"] == scale)
                & (offset_result["gene"] == gene)
            ]
            valid_group = group[group["support_sufficient"]]
            all_median, all_min, all_max, _ = _finite_summary(
                valid_group["n_grids_all"].to_numpy(float)
            )
            nz_median, _, _, _ = _finite_summary(
                valid_group["n_grids_nz"].to_numpy(float)
            )
            row = {
                "gene": gene,
                "grid_size_um": scale,
                "n_offsets_requested": expected_offsets,
                "n_offsets_support_sufficient": int(valid_group.shape[0]),
                "n_grids_all_median": all_median,
                "n_grids_all_min": all_min,
                "n_grids_all_max": all_max,
                "n_grids_nz_median": nz_median,
            }
            for column in metric_columns:
                values = valid_group[column].to_numpy(float)
                median, minimum, maximum, standard_deviation = _finite_summary(
                    values
                )
                row[column] = median
                row[f"{column}_n_offsets_valid"] = int(np.isfinite(values).sum())
                row[f"{column}_offset_min"] = minimum
                row[f"{column}_offset_max"] = maximum
                row[f"{column}_offset_sd"] = standard_deviation
            scale_rows.append(row)
    scale_result = pd.DataFrame.from_records(scale_rows)

    multiscale_rows: list[dict[str, Any]] = []
    for gene in scored_genes:
        group = scale_result[scale_result["gene"] == gene]
        support_scales = group["n_offsets_support_sufficient"].to_numpy(int) > 0
        row = {
            "gene": gene,
            "n_scales_requested": len(scales),
            "n_scales_support_sufficient": int(support_scales.sum()),
        }
        for column in (
            "n_grids_all_median", "n_grids_nz_median", *metric_columns
        ):
            median, minimum, maximum, standard_deviation = _finite_summary(
                group.loc[support_scales, column].to_numpy(float)
            )
            row[column] = median
            if column in metric_columns:
                row[f"{column}_n_scales_valid"] = int(
                    np.isfinite(group.loc[support_scales, column].to_numpy(float)).sum()
                )
                row[f"{column}_scale_min"] = minimum
                row[f"{column}_scale_max"] = maximum
                row[f"{column}_scale_sd"] = standard_deviation
        multiscale_rows.append(row)
    multiscale_result = pd.DataFrame.from_records(multiscale_rows)

    common_attrs = {
        "grid_sizes_um": list(scales),
        "metrics": list(metric_names),
        "expression_layer": expression_layer,
        "reference_coordinate_source": reference_coordinate_source,
        "moving_coordinate_source": moving_coordinate_source,
        "n_genes_requested": len(requested_genes),
        "n_genes_scored": len(scored_genes),
        "requested_genes": requested_genes,
        "scored_genes": scored_genes,
        "missing_genes": missing_genes,
        "min_cells_reference": int(min_cells_reference),
        "min_cells_moving": int(min_cells_moving),
        "min_overlap_grids": int(min_overlap_grids),
        "random_state": int(random_state),
    }
    for result in (multiscale_result, scale_result, offset_result):
        result.attrs.update(common_attrs)
    return multiscale_result, scale_result, offset_result


def summarize_pair_gene_similarity(
    scale_gene_result: pd.DataFrame,
    *,
    primary_metric: str = "pcc_all",
    min_gene_nz_grids: int = 5,
    min_gene_valid_offsets: int = 1,
    min_gene_valid_scales: int = 1,
) -> dict[str, Any]:
    """Summarize reliable gene scores without removing detailed-result rows."""
    if min_gene_nz_grids < 0:
        raise ValueError("min_gene_nz_grids must be nonnegative")
    if min_gene_valid_offsets < 1:
        raise ValueError("min_gene_valid_offsets must be at least 1")
    if min_gene_valid_scales < 1:
        raise ValueError("min_gene_valid_scales must be at least 1")
    metric_offset_count = f"{primary_metric}_n_offsets_valid"
    required = {
        "gene", "grid_size_um", "n_offsets_support_sufficient",
        "n_grids_all_median", "n_grids_nz_median", primary_metric,
        metric_offset_count,
    }
    missing_columns = required - set(scale_gene_result.columns)
    if missing_columns:
        raise ValueError(f"scale_gene_result is missing columns: {sorted(missing_columns)}")
    scales = list(
        scale_gene_result.attrs.get(
            "grid_sizes_um", scale_gene_result["grid_size_um"].drop_duplicates()
        )
    )
    work = scale_gene_result.copy()
    work["primary_metric_scale_eligible"] = (
        np.isfinite(work[primary_metric].to_numpy(float))
        & (work["n_grids_nz_median"].to_numpy(float) >= min_gene_nz_grids)
        & (work[metric_offset_count].to_numpy(int) >= min_gene_valid_offsets)
    )
    eligible_scale_counts = work.groupby("gene", sort=False)[
        "primary_metric_scale_eligible"
    ].sum()
    eligible_genes = [
        gene for gene in work["gene"].drop_duplicates()
        if int(eligible_scale_counts.get(gene, 0)) >= min_gene_valid_scales
    ]
    eligible_gene_set = set(eligible_genes)
    finite_by_gene = work.groupby("gene", sort=False)[primary_metric].agg(
        lambda values: bool(np.isfinite(values.to_numpy(float)).any())
    )

    per_scale_score: dict[float, float] = {}
    per_scale_n_genes_eligible: dict[float, int] = {}
    per_scale_grids_all: dict[float, float] = {}
    per_scale_grids_nz: dict[float, float] = {}
    per_scale_offsets_support: dict[float, int] = {}
    for scale in scales:
        group = work[np.isclose(work["grid_size_um"], scale)]
        eligible = group[
            group["primary_metric_scale_eligible"]
            & group["gene"].isin(eligible_gene_set)
        ]
        per_scale_score[float(scale)] = _finite_summary(
            eligible[primary_metric].to_numpy(float)
        )[0]
        per_scale_n_genes_eligible[float(scale)] = int(eligible["gene"].nunique())
        per_scale_grids_all[float(scale)] = _finite_summary(
            group["n_grids_all_median"].to_numpy(float)
        )[0]
        per_scale_grids_nz[float(scale)] = _finite_summary(
            group["n_grids_nz_median"].to_numpy(float)
        )[0]
        support_offsets = group["n_offsets_support_sufficient"].to_numpy(int)
        per_scale_offsets_support[float(scale)] = (
            int(np.max(support_offsets)) if support_offsets.size else 0
        )

    n_genes_finite = int(finite_by_gene.sum())
    n_genes_eligible = len(eligible_genes)
    n_genes_requested = int(
        scale_gene_result.attrs.get(
            "n_genes_requested", scale_gene_result["gene"].nunique()
        )
    )
    multiscale_score = _finite_summary(list(per_scale_score.values()))[0]
    return {
        "primary_metric": primary_metric,
        "grid_sizes_um": [float(scale) for scale in scales],
        "min_gene_nz_grids": int(min_gene_nz_grids),
        "min_gene_valid_offsets": int(min_gene_valid_offsets),
        "min_gene_valid_scales": int(min_gene_valid_scales),
        "per_scale_median_primary_score": per_scale_score,
        "multiscale_median_primary_score": multiscale_score,
        "n_genes_requested": n_genes_requested,
        "n_genes_scored": int(scale_gene_result["gene"].nunique()),
        "n_genes_primary_metric_finite": n_genes_finite,
        "n_genes_primary_metric_eligible": n_genes_eligible,
        "finite_gene_fraction": (
            float(n_genes_finite / n_genes_requested) if n_genes_requested else np.nan
        ),
        "eligible_gene_fraction": (
            float(n_genes_eligible / n_genes_requested) if n_genes_requested else np.nan
        ),
        "eligible_genes": eligible_genes,
        "per_scale_n_genes_eligible": per_scale_n_genes_eligible,
        "per_scale_median_n_grids_all": per_scale_grids_all,
        "per_scale_median_n_grids_nz": per_scale_grids_nz,
        "n_offsets_support_sufficient": per_scale_offsets_support,
    }


# ---------------------------------------------------------------------------
# Legacy wrappers
# ---------------------------------------------------------------------------

def compute_grid_similarity(
    slice1: Any,
    slice2: Any,
    tech1: str,
    tech2: str,
    coord_key: str = "spatial",
    gene_list: Optional[Sequence[str]] = None,
    grid_size_table: Optional[dict[tuple[str, str], float]] = None,
    min_fraction: float = 0.5,
) -> pd.DataFrame:
    """Deprecated one-scale wrapper preserving the old result schema."""
    warnings.warn(
        "compute_grid_similarity() is deprecated; use "
        "compute_pair_gene_similarity() for Xenium evaluation.",
        DeprecationWarning,
        stacklevel=2,
    )
    if slice1.n_obs > slice2.n_obs:
        slice1, slice2 = slice2, slice1
        tech1, tech2 = tech2, tech1
    table = GRID_SIZE_PX if grid_size_table is None else grid_size_table
    key = _pair_key(tech1, tech2)
    if key not in table:
        raise KeyError(f"No grid size registered for technology pair {key}")
    grid_size = float(table[key])
    coordinates1 = _validate_coordinates(slice1.obsm[coord_key], slice1.n_obs, "slice1 coordinates")
    coordinates2 = _validate_coordinates(slice2.obsm[coord_key], slice2.n_obs, "slice2 coordinates")
    _, _, _, _, _, occupancy_metadata = assign_grid_labels(
        coordinates1,
        coordinates2,
        grid_size_um=grid_size,
        min_cells_reference=1,
        min_cells_moving=1,
    )
    threshold = max(
        1,
        int(
            np.ceil(
                min(
                    occupancy_metadata["reference_mean_occupancy"],
                    occupancy_metadata["moving_mean_occupancy"],
                )
                * min_fraction
            )
        ),
    )
    _, scale_result, _ = compute_pair_gene_similarity(
        slice1,
        slice2,
        reference_coord_key=coord_key,
        moving_coord_key=coord_key,
        gene_list=gene_list,
        grid_sizes_um=(grid_size,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=threshold,
        min_cells_moving=threshold,
        min_overlap_grids=2,
        metrics=SUPPORTED_METRICS,
    )
    legacy_columns = [
        "n_grids_all", "n_grids_nz", "pcc_all", "cos_sim_all", "ssim_all",
        "mi_all", "pcc_nz", "cos_sim_nz", "ssim_nz", "mi_nz",
    ]
    if (
        scale_result.empty
        or not scale_result["n_offsets_support_sufficient"].gt(0).any()
    ):
        result = pd.DataFrame(columns=legacy_columns)
        result.index.name = "gene"
    else:
        result = scale_result.rename(
            columns={
                "n_grids_all_median": "n_grids_all",
                "n_grids_nz_median": "n_grids_nz",
            }
        ).set_index("gene")[legacy_columns]
    result.attrs.update(
        {
            "cell_size_px": grid_size,
            "n_grids_total": (
                int(result["n_grids_all"].iloc[0]) if not result.empty else 0
            ),
            "tech_pair": key,
            "canonical_order": (tech1, tech2),
        }
    )
    return result


def compute_stack_similarity(
    adatas: list[Any],
    tech: str | Sequence[str],
    coord_key: str = "spatial",
    gene_list: Optional[Sequence[str]] = None,
    grid_size_table: Optional[dict[tuple[str, str], float]] = None,
    min_fraction: float = 0.5,
    min_grids_nz: int = 50,
) -> pd.DataFrame:
    """Legacy consecutive-stack summary retained for existing callers."""
    if len(adatas) < 2:
        raise ValueError("Need at least two slices")
    techs = [tech] * len(adatas) if isinstance(tech, str) else list(tech)
    if len(techs) != len(adatas):
        raise ValueError("The number of technologies must equal the number of slices")
    if gene_list is None:
        genes = [
            str(gene) for gene in adatas[0].var_names
            if all(str(gene) in set(map(str, value.var_names)) for value in adatas[1:])
        ]
    else:
        genes = list(map(str, gene_list))
    pair_results = [
        compute_grid_similarity(
            adatas[index],
            adatas[index + 1],
            tech1=techs[index],
            tech2=techs[index + 1],
            coord_key=coord_key,
            gene_list=genes,
            grid_size_table=grid_size_table,
            min_fraction=min_fraction,
        )
        for index in range(len(adatas) - 1)
    ]
    metric_all = ["pcc_all", "cos_sim_all", "ssim_all", "mi_all"]
    metric_nz = ["pcc_nz", "cos_sim_nz", "ssim_nz", "mi_nz"]
    rows = []
    for gene in genes:
        available = [result.loc[gene] for result in pair_results if gene in result.index]
        row: dict[str, Any] = {"gene": gene, "n_pairs_total": len(pair_results)}
        for column in metric_all:
            row[column] = _finite_summary([value[column] for value in available])[0]
        contributing = [value for value in available if value["n_grids_nz"] >= min_grids_nz]
        for column in metric_nz:
            row[column] = _finite_summary([value[column] for value in contributing])[0]
        row["n_grids_all"] = _finite_summary(
            [value["n_grids_all"] for value in available]
        )[0]
        row["n_grids_nz"] = _finite_summary(
            [value["n_grids_nz"] for value in contributing]
        )[0]
        row["n_contributing_pairs"] = len(contributing)
        rows.append(row)
    return pd.DataFrame.from_records(rows).set_index("gene")


# ---------------------------------------------------------------------------
# Synthetic validation
# ---------------------------------------------------------------------------

def run_self_tests() -> None:
    """Run deterministic synthetic validation of the reusable API."""
    import anndata as ad

    exact_boundary = np.array(
        [(x, y) for y in (0.0, 1000.0) for x in (0.0, 1000.0)]
    )
    _, _, _, _, _, exact_metadata = assign_grid_labels(
        exact_boundary, exact_boundary, grid_size_um=200.0,
        min_cells_reference=1, min_cells_moving=1,
    )
    assert exact_metadata["n_cols"] == 5
    assert exact_metadata["n_rows"] == 5
    nondivisible = np.array(
        [(x, y) for y in (0.0, 950.0) for x in (0.0, 950.0)]
    )
    _, _, _, _, _, nondivisible_metadata = assign_grid_labels(
        nondivisible, nondivisible, grid_size_um=200.0,
        min_cells_reference=1, min_cells_moving=1,
    )
    assert nondivisible_metadata["n_cols"] == 5
    assert nondivisible_metadata["n_rows"] == 5

    coordinates = np.array(
        [(x, y) for y in np.arange(0.0, 1200.0, 200.0)
         for x in np.arange(0.0, 1200.0, 200.0)],
        dtype=float,
    )
    gradient = coordinates[:, 0] + 2.0 * coordinates[:, 1] + 1.0
    expression = np.column_stack(
        [gradient, np.ones(len(coordinates)), np.zeros(len(coordinates))]
    )
    genes = ["gradient", "constant", "zero"]
    reference = ad.AnnData(expression.copy())
    moving = ad.AnnData(expression.copy())
    reference.var_names = genes
    moving.var_names = genes
    reference.obsm["spatial"] = coordinates.copy()
    moving.obsm["candidate"] = coordinates.copy()
    reference.layers["normcounts"] = expression.copy()
    moving.layers["normcounts"] = expression.copy()
    original_reference_coordinates = reference.obsm["spatial"].copy()

    perfect_multi, perfect_scale, perfect_offset = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coord_key="spatial",
        moving_coord_key="candidate",
        gene_list=genes,
        expression_layer="normcounts",
        grid_sizes_um=DEFAULT_GRID_SIZES_UM,
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
    )
    perfect_gradient = perfect_multi.set_index("gene").loc["gradient"]
    assert np.isclose(perfect_gradient["pcc_all"], 1.0)
    assert np.isclose(perfect_gradient["cos_sim_all"], 1.0)
    assert np.isclose(perfect_gradient["ssim_all"], 1.0)
    assert np.isnan(perfect_multi.set_index("gene").loc["zero", "pcc_all"])
    assert perfect_scale["grid_size_um"].drop_duplicates().tolist() == [200.0, 300.0, 400.0]
    assert np.array_equal(reference.obsm["spatial"], original_reference_coordinates)
    perfect_scale_indexed = perfect_scale.set_index(["gene", "grid_size_um"])
    for metric_column in (
        "pcc_all", "pcc_nz", "cos_sim_all", "cos_sim_nz",
        "ssim_all", "ssim_nz", "mi_all", "mi_nz",
    ):
        assert f"{metric_column}_n_offsets_valid" in perfect_scale.columns
        assert f"{metric_column}_n_scales_valid" in perfect_multi.columns
    for low_information_gene in ("constant", "zero"):
        row = perfect_scale_indexed.loc[(low_information_gene, 200.0)]
        assert row["n_offsets_support_sufficient"] > 0
        assert np.isnan(row["pcc_all"])
        assert row["pcc_all_n_offsets_valid"] == 0
    assert perfect_scale_indexed.loc[("gradient", 200.0), "pcc_all_n_offsets_valid"] == 1
    required_offset_diagnostics = {
        "x_min", "x_max", "y_min", "y_max", "grid_origin_x_um",
        "grid_origin_y_um", "n_rows", "n_cols", "n_valid_grids",
        "reference_cells_in_valid_grids", "moving_cells_in_valid_grids",
        "reference_mean_cells_per_valid_grid",
        "moving_median_cells_per_valid_grid",
    }
    assert required_offset_diagnostics <= set(perfect_offset.columns)

    translation = np.array([123.0, -77.0])
    translated_multi, _, _ = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates + translation,
        moving_coords=coordinates + translation,
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc", "ssim"),
    )
    base_value = perfect_scale.loc[
        (perfect_scale["gene"] == "gradient")
        & (perfect_scale["grid_size_um"] == 200.0), "pcc_all"
    ].iloc[0]
    assert np.isclose(translated_multi.loc[0, "pcc_all"], base_value)
    assert "cos_sim_all" not in translated_multi.columns

    moving.obsm["wrong"] = coordinates + np.array([10000.0, 10000.0])
    explicit_multi, _, _ = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates,
        moving_coords=coordinates,
        reference_coord_key="spatial",
        moving_coord_key="wrong",
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc",),
    )
    assert np.isclose(explicit_multi.loc[0, "pcc_all"], 1.0)
    assert explicit_multi.attrs["moving_coordinate_source"] == "explicit"

    rng = np.random.default_rng(7)
    random_moving = moving.copy()
    random_moving.X = expression[rng.permutation(len(expression))]
    random_multi, _, _ = compute_pair_gene_similarity(
        reference,
        random_moving,
        reference_coords=coordinates,
        moving_coords=coordinates,
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc",),
    )
    assert random_multi.loc[0, "pcc_all"] < 0.95

    partial_multi, _, partial_offset = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates,
        moving_coords=coordinates + np.array([400.0, 0.0]),
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=5,
        metrics=("pcc",),
    )
    assert partial_offset.loc[0, "support_sufficient"]
    assert np.isfinite(partial_multi.loc[0, "pcc_all"])

    insufficient_multi, _, insufficient_offset = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates,
        moving_coords=coordinates + np.array([950.0, 0.0]),
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=10,
        metrics=("pcc",),
    )
    assert not insufficient_offset.loc[0, "support_sufficient"]
    assert np.isnan(insufficient_multi.loc[0, "pcc_all"])

    no_overlap_multi, _, no_overlap_offset = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates,
        moving_coords=coordinates + np.array([2000.0, 0.0]),
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc",),
    )
    assert not no_overlap_offset.loc[0, "intersection_valid"]
    assert not no_overlap_offset.loc[0, "support_sufficient"]
    assert np.isnan(no_overlap_multi.loc[0, "pcc_all"])

    sparse_reference = reference.copy()
    sparse_moving = moving.copy()
    sparse_reference.X = sp.csr_matrix(reference.X)
    sparse_moving.X = sp.csc_matrix(moving.X)
    sparse_multi, _, _ = compute_pair_gene_similarity(
        sparse_reference,
        sparse_moving,
        reference_coords=coordinates,
        moving_coords=coordinates,
        gene_list=genes,
        grid_sizes_um=(200.0,),
        grid_offsets_um=((0.0, 0.0),),
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc", "cosine", "ssim"),
    )
    dense_single = perfect_scale[perfect_scale["grid_size_um"] == 200.0].set_index("gene")
    sparse_single = sparse_multi.set_index("gene")
    for column in ("pcc_all", "cos_sim_all", "ssim_all"):
        assert np.allclose(
            dense_single.loc["gradient", column],
            sparse_single.loc["gradient", column],
            equal_nan=True,
        )

    _, offset_scale, multiple_offsets = compute_pair_gene_similarity(
        reference,
        moving,
        reference_coords=coordinates,
        moving_coords=coordinates,
        gene_list=("gradient",),
        grid_sizes_um=(200.0,),
        grid_offsets_um=None,
        min_cells_reference=1,
        min_cells_moving=1,
        min_overlap_grids=2,
        metrics=("pcc",),
    )
    assert len(multiple_offsets) == 4
    assert offset_scale.loc[0, "n_offsets_requested"] == 4
    assert offset_scale.loc[0, "n_offsets_support_sufficient"] == 4
    assert offset_scale.loc[0, "pcc_all_n_offsets_valid"] == 4
    assert np.isclose(
        offset_scale.loc[0, "pcc_all"],
        np.nanmedian(multiple_offsets["pcc_all"]),
    )

    fake_rows = []
    for gene, scores, nz_grids, valid_offsets in (
        ("good", (0.1, 0.9, 0.4), (10, 10, 10), (2, 2, 2)),
        ("few_nz", (0.8, 0.8, 0.8), (2, 2, 2), (2, 2, 2)),
        ("few_offsets", (0.7, 0.7, 0.7), (10, 10, 10), (0, 0, 0)),
        ("one_scale", (0.6, np.nan, np.nan), (10, 10, 10), (2, 0, 0)),
    ):
        for scale, score, n_nz, n_valid_offsets in zip(
            (200.0, 300.0, 400.0), scores, nz_grids, valid_offsets
        ):
            fake_rows.append(
                {
                    "gene": gene, "grid_size_um": scale,
                    "n_offsets_support_sufficient": 2,
                    "n_grids_all_median": 20.0,
                    "n_grids_nz_median": n_nz,
                    "pcc_all": score,
                    "pcc_all_n_offsets_valid": n_valid_offsets,
                }
            )
    fake_scale = pd.DataFrame.from_records(fake_rows)
    fake_scale.attrs["grid_sizes_um"] = [200.0, 300.0, 400.0]
    summary = summarize_pair_gene_similarity(
        fake_scale, min_gene_nz_grids=5, min_gene_valid_offsets=1,
        min_gene_valid_scales=2,
    )
    assert np.isclose(summary["multiscale_median_primary_score"], 0.4)
    assert not np.isclose(summary["multiscale_median_primary_score"], 0.9)
    assert summary["eligible_genes"] == ["good"]
    assert summary["n_genes_primary_metric_finite"] == 4
    assert summary["n_genes_primary_metric_eligible"] == 1
    relaxed_scales = summarize_pair_gene_similarity(
        fake_scale, min_gene_nz_grids=5, min_gene_valid_offsets=1,
        min_gene_valid_scales=1,
    )
    assert relaxed_scales["eligible_genes"] == ["good", "one_scale"]
    relaxed_nz = summarize_pair_gene_similarity(
        fake_scale, min_gene_nz_grids=0, min_gene_valid_offsets=1,
        min_gene_valid_scales=2,
    )
    assert "few_nz" in relaxed_nz["eligible_genes"]
    assert "few_offsets" not in relaxed_nz["eligible_genes"]
    assert set(fake_scale["gene"]) == {"good", "few_nz", "few_offsets", "one_scale"}
    for invalid_kwargs in (
        {"min_gene_nz_grids": -1},
        {"min_gene_valid_offsets": 0},
        {"min_gene_valid_scales": 0},
    ):
        try:
            summarize_pair_gene_similarity(fake_scale, **invalid_kwargs)
        except ValueError:
            pass
        else:
            raise AssertionError(f"Invalid summary threshold accepted: {invalid_kwargs}")

    try:
        compute_pair_gene_similarity(
            reference, moving, reference_coords=coordinates,
            moving_coords=coordinates, gene_list=("gradient",),
            grid_sizes_um=(200.0,),
            grid_offsets_um=((0.0, 0.0), (200.0, 0.0)),
            min_cells_reference=1, min_cells_moving=1,
            min_overlap_grids=2, metrics=("pcc",),
        )
    except ValueError as exc:
        assert "Duplicate grid offsets" in str(exc)
    else:
        raise AssertionError("Duplicate normalized offsets were not rejected")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            compute_pair_gene_similarity(
                reference, moving, reference_coords=coordinates,
                moving_coords=coordinates, gene_list=("not_shared",),
                grid_sizes_um=(200.0,), grid_offsets_um=((0.0, 0.0),),
                min_cells_reference=1, min_cells_moving=1,
                min_overlap_grids=2, metrics=("pcc",),
            )
        except ValueError as exc:
            assert "No requested genes are present" in str(exc)
        else:
            raise AssertionError("An empty shared-gene set was not rejected")

    original_dense_column = _dense_column
    extraction_count = 0

    def counting_dense_column(matrix: Any, column_index: int) -> np.ndarray:
        nonlocal extraction_count
        extraction_count += 1
        return original_dense_column(matrix, column_index)

    globals()["_dense_column"] = counting_dense_column
    try:
        compute_pair_gene_similarity(
            reference, moving, reference_coords=coordinates,
            moving_coords=coordinates, gene_list=("gradient", "constant"),
            grid_sizes_um=(200.0, 300.0), grid_offsets_um=None,
            min_cells_reference=1, min_cells_moving=1,
            min_overlap_grids=2, metrics=("pcc",),
        )
    finally:
        globals()["_dense_column"] = original_dense_column
    assert extraction_count == 2 * 2

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        missing_multi, _, _ = compute_pair_gene_similarity(
            reference,
            moving,
            reference_coords=coordinates,
            moving_coords=coordinates,
            gene_list=("gradient", "absent"),
            grid_sizes_um=(200.0,),
            grid_offsets_um=((0.0, 0.0),),
            min_cells_reference=1,
            min_cells_moving=1,
            min_overlap_grids=2,
            metrics=("pcc",),
        )
    assert missing_multi["gene"].tolist() == ["gradient"]
    assert missing_multi.attrs["missing_genes"] == ["absent"]
    assert caught

    print("All feature_similarity self-tests passed")


if __name__ == "__main__":
    run_self_tests()
