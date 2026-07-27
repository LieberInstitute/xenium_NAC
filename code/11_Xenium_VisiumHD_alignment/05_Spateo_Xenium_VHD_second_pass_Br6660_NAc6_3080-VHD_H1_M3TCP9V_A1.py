#!/usr/bin/env python3
"""Pair-specific second-pass alignment of VHD to fixed Xenium.

This pair uses its first-pass coordinates to identify the cross-modality tissue
overlap. Only Xenium and VHD ``WM`` observations inside that overlap are used
to fit each residual rigid transformation. The first result is saved as
``intermediate``; overlap is recalculated from those coordinates before a
second refinement is fitted and saved as the final ``second_pass`` result.
Each fitted transformation is applied to the complete VHD slice.

The script also saves the standard Module 02 outputs and plots, plus separate
first-pass-input overviews for all Xenium CellType and VHD Spatial_Domain
labels.
"""

from __future__ import annotations

import json
import math
import os
import pickle
import warnings
from pathlib import Path

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")
warnings.filterwarnings(
    "ignore",
    message="The pynvml package is deprecated.*",
    category=FutureWarning,
)

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import spateo as st
import torch
from scipy.spatial import cKDTree


# =============================================================================
# PAIR-SPECIFIC CONFIGURATION
# =============================================================================
DONOR = "Br6660"
XENIUM_SAMPLE = "Br6660_NAc6_3080"
VHD_SAMPLE = "VHD_H1_M3TCP9V_A1"
VHD_Z_HEIGHT = 3040.0

REALIGN_LABEL = "WM"
INTERMEDIATE_OVERLAP_DISTANCE_UM = 100.0
FINAL_OVERLAP_DISTANCE_UM = 300.0

DEVICE = "cuda"
MAX_ITER = 200
PARTIAL_ROBUST_LEVEL = 100
# =============================================================================

ROOT = Path(__file__).resolve().parents[2]
PAIR_TAG = f"{XENIUM_SAMPLE}-{VHD_SAMPLE}"

XENIUM_H5AD = (
    ROOT
    / "processed-data"
    / "02_build_spe"
    / "h5ad"
    / "spe_NormCounts_nucleus_normcounts.h5ad"
)
CELLTYPE_CSV = (
    ROOT / "processed-data" / "05_Clustering" / "Banksy_cell_types.csv"
)
XENIUM_COORDINATES = (
    ROOT
    / "processed-data"
    / "05_Xenium_alignment"
    / f"module02_Spateo_second_pass_{DONOR}"
    / f"{DONOR}_second_pass_coordinates.csv"
)
VHD_H5AD = (
    ROOT
    / "processed-data"
    / "HD_Full_Analysis"
    / "h5ad"
    / f"{VHD_SAMPLE}_spaceranger_square008um_spatial_um.h5ad"
)
VHD_FIRST_PASS_COORDINATES = (
    ROOT
    / "processed-data"
    / "11_Xenium_VisiumHD_alignment"
    / "module_01_Spateo_first_pass"
    / PAIR_TAG
    / f"{PAIR_TAG}_first_pass_coordinates.csv"
)
PLOT_DIR = (
    ROOT
    / "plots"
    / "11_Xenium_VisiumHD_alignment"
    / "module_02_Spateo_second_pass"
    / PAIR_TAG
)
OUTPUT_DIR = (
    ROOT
    / "processed-data"
    / "11_Xenium_VisiumHD_alignment"
    / "module_02_Spateo_second_pass"
    / PAIR_TAG
)

TRANSFORMATION_PATH = OUTPUT_DIR / f"{PAIR_TAG}_second_pass.pkl"
COORDINATE_PATH = OUTPUT_DIR / f"{PAIR_TAG}_second_pass_coordinates.csv"
METADATA_PATH = OUTPUT_DIR / f"{PAIR_TAG}_second_pass_metadata.json"
CELLTYPE_QC_PATH = OUTPUT_DIR / f"{PAIR_TAG}_realign_celltype_qc.csv"

INTERMEDIATE_TRANSFORMATION_PATH = (
    OUTPUT_DIR / f"{PAIR_TAG}_intermediate.pkl"
)
INTERMEDIATE_COORDINATE_PATH = (
    OUTPUT_DIR / f"{PAIR_TAG}_intermediate_coordinates.csv"
)
INTERMEDIATE_METADATA_PATH = (
    OUTPUT_DIR / f"{PAIR_TAG}_intermediate_metadata.json"
)
INTERMEDIATE_CELLTYPE_QC_PATH = (
    OUTPUT_DIR / f"{PAIR_TAG}_intermediate_realign_celltype_qc.csv"
)

FIRST_PASS_OVERLAY_PATH = PLOT_DIR / f"{PAIR_TAG}_first_pass_input_um.png"
INTERMEDIATE_OVERLAY_PATH = (
    PLOT_DIR / f"{PAIR_TAG}_intermediate_spateo_um.png"
)
SECOND_PASS_OVERLAY_PATH = PLOT_DIR / f"{PAIR_TAG}_second_spateo_um.png"
XENIUM_SANITY_PATH = (
    PLOT_DIR / f"{DONOR}_spateo_Xenium_sanity_check.png"
)

XENIUM_COLOR = "#4C78A8"
VHD_COLOR = "#F58518"
POINT_SIZE = 0.6
POINT_ALPHA = 0.72
N_COLUMNS = 4
SHARED_LABEL_PLOTS = (
    "D1_Island_A",
    "D1_Island_B",
    "Ependymal",
    "Excitatory",
    "WM",
)


def validate_input_paths() -> None:
    required_paths = (
        XENIUM_H5AD,
        CELLTYPE_CSV,
        XENIUM_COORDINATES,
        VHD_H5AD,
        VHD_FIRST_PASS_COORDINATES,
    )
    missing = [path for path in required_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(
            "Required plotting input(s) do not exist: "
            + ", ".join(map(str, missing))
        )


def load_xenium() -> tuple[ad.AnnData, list[ad.AnnData], float]:
    """Load aligned donor Xenium data, target slice, and physical depth."""
    xenium_all = ad.read_h5ad(XENIUM_H5AD)
    required_obs = {"Donor", "Sample", "cell_id"}
    missing_obs = required_obs - set(xenium_all.obs.columns)
    if missing_obs:
        raise ValueError(
            f"Xenium H5AD is missing obs columns: {sorted(missing_obs)}"
        )

    donor_mask = xenium_all.obs["Donor"].astype(str).eq(DONOR)
    donor_xenium = xenium_all[donor_mask].copy()
    if donor_xenium.n_obs == 0:
        raise ValueError(f"No Xenium cells found for {DONOR}")

    celltypes = pd.read_csv(
        CELLTYPE_CSV,
        usecols=["cell_id", "CellType"],
        dtype={"cell_id": "string"},
    )
    if celltypes["cell_id"].isna().any() or celltypes["cell_id"].duplicated().any():
        raise ValueError("Banksy cell-type CSV contains invalid cell IDs")
    celltype_map = celltypes.set_index("cell_id")["CellType"]
    xenium_ids = donor_xenium.obs["cell_id"].astype("string")
    donor_xenium.obs["CellType"] = xenium_ids.map(celltype_map)

    coordinates = pd.read_csv(
        XENIUM_COORDINATES,
        usecols=[
            "donor",
            "sample",
            "cell_id",
            "x_second_pass",
            "y_second_pass",
            "z_height",
        ],
        dtype={"donor": "string", "sample": "string", "cell_id": "string"},
    )
    coordinates = coordinates[coordinates["donor"].eq(DONOR)].copy()
    if coordinates["cell_id"].isna().any() or coordinates["cell_id"].duplicated().any():
        raise ValueError("Final Xenium coordinate CSV contains invalid cell IDs")
    coordinate_map = coordinates.set_index("cell_id")
    missing_ids = xenium_ids[~xenium_ids.isin(coordinate_map.index)]
    if len(missing_ids):
        raise ValueError(
            f"{len(missing_ids)} Xenium cells are missing final coordinates"
        )
    donor_xenium.obsm["plot_spatial"] = coordinate_map.loc[
        xenium_ids, ["x_second_pass", "y_second_pass"]
    ].to_numpy(dtype=float)
    donor_xenium.obsm["align_spatial"] = np.asarray(
        donor_xenium.obsm["plot_spatial"], dtype=float
    ).copy()

    sample_depths = (
        coordinates.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
    )
    ordered_slices = []
    for sample in sample_depths.index.astype(str):
        mask = donor_xenium.obs["Sample"].astype(str).eq(sample)
        ordered_slices.append(donor_xenium[mask].copy())

    target_mask = donor_xenium.obs["Sample"].astype(str).eq(XENIUM_SAMPLE)
    target = donor_xenium[target_mask].copy()
    if target.n_obs == 0:
        raise ValueError(f"No Xenium cells found for {XENIUM_SAMPLE}")
    target_depth = sample_depths.loc[XENIUM_SAMPLE]
    target.obsm["spatial"] = np.asarray(
        target.obsm["plot_spatial"], dtype=float
    ).copy()
    return target, ordered_slices, float(target_depth)


def load_vhd() -> ad.AnnData:
    """Load one VHD slice with its Module 01 first-pass coordinates."""
    vhd = ad.read_h5ad(VHD_H5AD)
    if "Spatial_Domain" not in vhd.obs:
        raise ValueError("VHD H5AD is missing obs['Spatial_Domain']")

    coordinates = pd.read_csv(
        VHD_FIRST_PASS_COORDINATES,
        usecols=[
            "donor",
            "xenium_sample",
            "vhd_sample",
            "vhd_obs_id",
            "x_first_pass",
            "y_first_pass",
        ],
        dtype={
            "donor": "string",
            "xenium_sample": "string",
            "vhd_sample": "string",
            "vhd_obs_id": "string",
        },
    )
    expected = {
        "donor": DONOR,
        "xenium_sample": XENIUM_SAMPLE,
        "vhd_sample": VHD_SAMPLE,
    }
    for column, value in expected.items():
        observed = set(coordinates[column].dropna().astype(str))
        if observed != {value}:
            raise ValueError(
                f"First-pass CSV {column} mismatch: {sorted(observed)}"
            )
    if (
        coordinates["vhd_obs_id"].isna().any()
        or coordinates["vhd_obs_id"].duplicated().any()
    ):
        raise ValueError("First-pass VHD coordinate CSV contains invalid IDs")

    coordinate_map = coordinates.set_index("vhd_obs_id")
    vhd_ids = pd.Index(vhd.obs_names.astype(str))
    missing_ids = vhd_ids[~vhd_ids.isin(coordinate_map.index)]
    unexpected_ids = coordinate_map.index[
        ~coordinate_map.index.isin(vhd_ids)
    ]
    if len(missing_ids) or len(unexpected_ids):
        raise ValueError(
            "VHD ID mismatch between H5AD and first-pass CSV: "
            f"missing={len(missing_ids)}, unexpected={len(unexpected_ids)}"
        )
    vhd.obsm["plot_spatial"] = coordinate_map.loc[
        vhd_ids, ["x_first_pass", "y_first_pass"]
    ].to_numpy(dtype=float)
    vhd.obsm["spatial"] = np.asarray(
        vhd.obsm["plot_spatial"], dtype=float
    ).copy()
    return vhd


def ordered_labels(values: pd.Series) -> list[str]:
    """Return nonmissing labels in deterministic, case-insensitive order."""
    labels = values.astype("string").dropna().astype(str)
    labels = labels[labels.str.strip().ne("")]
    return sorted(labels.unique().tolist(), key=str.casefold)


def padded_limits(values: np.ndarray) -> tuple[float, float]:
    lower = float(np.nanmin(values))
    upper = float(np.nanmax(values))
    padding = max((upper - lower) * 0.025, 1.0)
    return lower - padding, upper + padding


def plot_label_subplots(
    adata,
    *,
    spatial_key: str,
    label_key: str,
    modality_name: str,
    color: str,
    output_path: Path,
) -> None:
    """Save one figure containing one spatial subplot per annotation label."""
    labels = ordered_labels(adata.obs[label_key])
    if not labels:
        raise ValueError(f"No usable labels found in obs[{label_key!r}]")

    xy = np.asarray(adata.obsm[spatial_key], dtype=float)[:, :2]
    if not np.isfinite(xy).all():
        raise ValueError(f"{modality_name} coordinates contain nonfinite values")
    x_limits = padded_limits(xy[:, 0])
    y_limits = padded_limits(xy[:, 1])

    n_columns = min(N_COLUMNS, len(labels))
    n_rows = math.ceil(len(labels) / n_columns)
    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(4.0 * n_columns, 3.7 * n_rows),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    label_values = adata.obs[label_key].astype("string")

    for axis, label in zip(axes.flat, labels):
        mask = (
            label_values.eq(label)
            .fillna(False)
            .to_numpy(dtype=bool)
        )
        label_xy = xy[mask]
        axis.scatter(
            label_xy[:, 0],
            label_xy[:, 1],
            s=POINT_SIZE,
            alpha=POINT_ALPHA,
            color=color,
            linewidths=0,
            label=modality_name,
        )
        axis.set_xlim(x_limits)
        axis.set_ylim(y_limits)
        axis.set_aspect("equal", adjustable="box")
        axis.set_title(f"{label} (n={len(label_xy):,})", fontsize=9)
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.legend(
            loc="upper right",
            fontsize=6,
            markerscale=7,
            frameon=False,
        )

    for axis in axes.flat[len(labels):]:
        axis.axis("off")

    figure.suptitle(
        f"{modality_name} first-pass alignment input by {label_key}",
        fontsize=13,
        y=1.002,
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(f"Saved {modality_name} {label_key} overview: {output_path}")


def plot_pair_overlay(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    spatial_key: str,
    output_path: Path,
    stage_label: str,
) -> None:
    """Plot all Xenium and VHD observations in one pairwise overlay."""
    xenium_xy = np.asarray(xenium.obsm[spatial_key], dtype=float)[:, :2]
    vhd_xy = np.asarray(vhd.obsm[spatial_key], dtype=float)[:, :2]
    figure, axis = plt.subplots(figsize=(5.0, 4.5))
    axis.scatter(
        xenium_xy[:, 0],
        xenium_xy[:, 1],
        s=0.2,
        alpha=0.45,
        color=XENIUM_COLOR,
        linewidths=0,
        label=XENIUM_SAMPLE,
    )
    axis.scatter(
        vhd_xy[:, 0],
        vhd_xy[:, 1],
        s=0.2,
        alpha=0.45,
        color=VHD_COLOR,
        linewidths=0,
        label=VHD_SAMPLE,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    axis.set_title(
        f"{stage_label}: {XENIUM_SAMPLE} → {VHD_SAMPLE}",
        fontsize=9,
    )
    axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_xenium_sanity(
    slices: list[ad.AnnData],
    output_path: Path,
) -> None:
    """Plot the standard consecutive final-Xenium slice overlays."""
    n_panels = len(slices) - 1
    n_columns = min(3, n_panels)
    n_rows = math.ceil(n_panels / n_columns)
    figure, axes = plt.subplots(
        n_rows,
        n_columns,
        figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )
    names = [
        str(model.obs["Sample"].astype(str).iloc[0]) for model in slices
    ]
    for index, axis in enumerate(axes.flat[:n_panels]):
        reference_xy = np.asarray(
            slices[index].obsm["align_spatial"], dtype=float
        )[:, :2]
        moving_xy = np.asarray(
            slices[index + 1].obsm["align_spatial"], dtype=float
        )[:, :2]
        axis.scatter(
            reference_xy[:, 0],
            reference_xy[:, 1],
            s=0.2,
            alpha=0.45,
            color=XENIUM_COLOR,
            linewidths=0,
            label=names[index],
        )
        axis.scatter(
            moving_xy[:, 0],
            moving_xy[:, 1],
            s=0.2,
            alpha=0.45,
            color=VHD_COLOR,
            linewidths=0,
            label=names[index + 1],
        )
        axis.set_aspect("equal")
        axis.set_xlabel("x (um)")
        axis.set_ylabel("y (um)")
        axis.set_title(
            f"Final Spateo aligned: {names[index]} → {names[index + 1]}",
            fontsize=9,
        )
        axis.legend(
            loc="upper right",
            fontsize=6,
            markerscale=8,
            frameon=False,
        )
    for axis in axes.flat[n_panels:]:
        axis.axis("off")
    figure.suptitle(
        "Final second-pass Xenium consecutive-slice overlays",
        y=1.002,
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_by_label_overlay(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    *,
    display_name: str,
    xenium_labels: tuple[str, ...],
    vhd_labels: tuple[str, ...],
    output_path: Path,
    stage_label: str = "Second-pass aligned",
) -> None:
    """Plot a selected cross-modality label correspondence."""
    xenium_mask = (
        xenium.obs["CellType"]
        .astype("string")
        .isin(xenium_labels)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    vhd_mask = (
        vhd.obs["Spatial_Domain"]
        .astype("string")
        .isin(vhd_labels)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    xenium_xy = np.asarray(
        xenium.obsm["align_spatial"], dtype=float
    )[xenium_mask, :2]
    vhd_xy = np.asarray(
        vhd.obsm["align_spatial"], dtype=float
    )[vhd_mask, :2]

    figure, axis = plt.subplots(figsize=(5.0, 4.5))
    axis.scatter(
        vhd_xy[:, 0],
        vhd_xy[:, 1],
        s=0.6,
        alpha=0.35,
        color=VHD_COLOR,
        linewidths=0,
        label=VHD_SAMPLE,
        zorder=1,
    )
    axis.scatter(
        xenium_xy[:, 0],
        xenium_xy[:, 1],
        s=0.6,
        alpha=0.60,
        color=XENIUM_COLOR,
        linewidths=0,
        label=XENIUM_SAMPLE,
        zorder=2,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    axis.set_title(
        f"{display_name}: {stage_label}: "
        f"{XENIUM_SAMPLE} → {VHD_SAMPLE}",
        fontsize=9,
    )
    axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)
    print(
        f"{display_name}: Xenium={len(xenium_xy):,}; "
        f"VHD={len(vhd_xy):,}"
    )


def plot_overlap_wm_overlay(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    *,
    spatial_key: str,
    stage_label: str,
    output_path: Path,
) -> None:
    """Plot exactly the current-overlap WM observations used for fitting."""
    xenium_mask = (
        xenium.obs["CellType"]
        .astype("string")
        .eq(REALIGN_LABEL)
        .fillna(False)
        .to_numpy(dtype=bool)
        & xenium.obs["_in_alignment_overlap"].to_numpy(dtype=bool)
    )
    vhd_mask = (
        vhd.obs["Spatial_Domain"]
        .astype("string")
        .eq(REALIGN_LABEL)
        .fillna(False)
        .to_numpy(dtype=bool)
        & vhd.obs["_in_alignment_overlap"].to_numpy(dtype=bool)
    )
    xenium_xy = np.asarray(
        xenium.obsm[spatial_key], dtype=float
    )[xenium_mask, :2]
    vhd_xy = np.asarray(
        vhd.obsm[spatial_key], dtype=float
    )[vhd_mask, :2]

    figure, axis = plt.subplots(figsize=(5.0, 4.5))
    axis.scatter(
        vhd_xy[:, 0],
        vhd_xy[:, 1],
        s=0.6,
        alpha=0.35,
        color=VHD_COLOR,
        linewidths=0,
        label=VHD_SAMPLE,
        zorder=1,
    )
    axis.scatter(
        xenium_xy[:, 0],
        xenium_xy[:, 1],
        s=0.6,
        alpha=0.60,
        color=XENIUM_COLOR,
        linewidths=0,
        label=XENIUM_SAMPLE,
        zorder=2,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    axis.set_title(
        f"WM in current overlap: {stage_label}",
        fontsize=9,
    )
    axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def prepare_alignment_subsets(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    overlap_mode: str,
    overlap_distance_um: float,
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    """Select WM using tissue-wide or WM-specific cross-modal overlap."""
    if overlap_mode not in {"tissue_then_WM", "WM_to_WM"}:
        raise ValueError(f"Unsupported overlap mode: {overlap_mode}")

    xenium_xy = np.asarray(xenium.obsm["plot_spatial"], dtype=float)[:, :2]
    vhd_xy = np.asarray(vhd.obsm["plot_spatial"], dtype=float)[:, :2]
    xenium_is_wm = (
        xenium.obs["CellType"]
        .astype("string")
        .eq(REALIGN_LABEL)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    vhd_is_wm = (
        vhd.obs["Spatial_Domain"]
        .astype("string")
        .eq(REALIGN_LABEL)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    if not xenium_is_wm.any() or not vhd_is_wm.any():
        raise ValueError("WM must occur in both modalities")

    if overlap_mode == "tissue_then_WM":
        xenium_overlap = (
            cKDTree(vhd_xy)
            .query(xenium_xy, k=1, workers=-1)[0]
            <= overlap_distance_um
        )
        vhd_overlap = (
            cKDTree(xenium_xy)
            .query(vhd_xy, k=1, workers=-1)[0]
            <= overlap_distance_um
        )
    else:
        xenium_overlap = np.zeros(xenium.n_obs, dtype=bool)
        vhd_overlap = np.zeros(vhd.n_obs, dtype=bool)
        xenium_overlap[xenium_is_wm] = (
            cKDTree(vhd_xy[vhd_is_wm])
            .query(xenium_xy[xenium_is_wm], k=1, workers=-1)[0]
            <= overlap_distance_um
        )
        vhd_overlap[vhd_is_wm] = (
            cKDTree(xenium_xy[xenium_is_wm])
            .query(vhd_xy[vhd_is_wm], k=1, workers=-1)[0]
            <= overlap_distance_um
        )

    xenium.obs["_in_alignment_overlap"] = xenium_overlap
    vhd.obs["_in_alignment_overlap"] = vhd_overlap

    vhd_filtered = vhd.copy()
    sc.pp.filter_cells(vhd_filtered, min_genes=10)

    xenium_mask = xenium_is_wm & xenium_overlap
    filtered_vhd_is_wm = (
        vhd_filtered.obs["Spatial_Domain"]
        .astype("string")
        .eq(REALIGN_LABEL)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    vhd_mask = (
        filtered_vhd_is_wm
        & vhd_filtered.obs["_in_alignment_overlap"].to_numpy(dtype=bool)
    )
    qc = pd.DataFrame(
        [
            {
                "alignment_label": REALIGN_LABEL,
                "overlap_mode": overlap_mode,
                "overlap_distance_um": overlap_distance_um,
                "n_xenium_total": xenium.n_obs,
                "n_vhd_total": vhd.n_obs,
                "n_xenium_in_overlap": int(xenium_overlap.sum()),
                "n_vhd_in_overlap": int(vhd_overlap.sum()),
                "n_xenium_WM_total": int(xenium_is_wm.sum()),
                "n_vhd_WM_total": int(vhd_is_wm.sum()),
                "n_xenium_WM_in_overlap": int(xenium_mask.sum()),
                "n_vhd_WM_in_overlap_after_filter": int(vhd_mask.sum()),
            }
        ]
    )
    if not xenium_mask.any() or not vhd_mask.any():
        raise ValueError(
            "WM must occur in the current overlap of both modalities:\n"
            + qc.to_string(index=False)
        )

    xenium_sub = xenium[xenium_mask].copy()
    vhd_sub = vhd_filtered[vhd_mask].copy()
    xenium_sub.obsm["spatial"] = np.asarray(
        xenium_sub.obsm["plot_spatial"], dtype=float
    ).copy()
    vhd_sub.obsm["spatial"] = np.asarray(
        vhd_sub.obsm["plot_spatial"], dtype=float
    ).copy()
    xenium_sub.obs["CellType"] = pd.Categorical(
        [REALIGN_LABEL] * xenium_sub.n_obs,
        categories=[REALIGN_LABEL],
    )
    vhd_sub.obs["CellType"] = pd.Categorical(
        [REALIGN_LABEL] * vhd_sub.n_obs,
        categories=[REALIGN_LABEL],
    )
    return xenium_sub, vhd_sub, qc


def apply_transformation(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    transformation,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Apply the residual transform while verifying Xenium stays fixed."""
    fixed = xenium.copy()
    moving = vhd.copy()
    for model in (fixed, moving):
        if "align_spatial" in model.obsm:
            del model.obsm["align_spatial"]
    aligned = st.align.morpho_align_apply_transformation(
        models=[fixed, moving],
        spatial_key="spatial",
        key_added="align_spatial",
        transformation=transformation,
    )
    if not np.allclose(
        np.asarray(aligned[0].obsm["align_spatial"]),
        np.asarray(aligned[0].obsm["spatial"]),
    ):
        raise RuntimeError("Fixed Xenium coordinates changed unexpectedly")
    return aligned[0], aligned[1]


def fit_residual_transformation(
    xenium_sub: ad.AnnData,
    vhd_sub: ad.AnnData,
    stage_name: str,
    overlap_distance_um: float,
):
    """Fit one overlap-restricted WM residual transformation."""
    print(
        f"Fitting {stage_name} overlap-restricted "
        f"{REALIGN_LABEL} ↔ {REALIGN_LABEL} residual transformation "
        f"({overlap_distance_um:g} um support)."
    )
    return st.align.morpho_align_transformation(
        models=[xenium_sub, vhd_sub],
        spatial_key="spatial",
        key_added="align_spatial",
        device=DEVICE,
        verbose=True,
        max_iter=MAX_ITER,
        partial_robust_level=PARTIAL_ROBUST_LEVEL,
        SVI_mode=False,
        rep_layer=["CellType"],
        rep_field=["obs"],
        dissimilarity=["label"],
    )


def prepare_next_round_inputs(
    xenium_aligned: ad.AnnData,
    vhd_aligned: ad.AnnData,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Use the preceding round's aligned coordinates as the next input."""
    xenium_next = xenium_aligned.copy()
    vhd_next = vhd_aligned.copy()
    for model in (xenium_next, vhd_next):
        current_xy = np.asarray(
            model.obsm["align_spatial"], dtype=float
        ).copy()
        model.obsm["plot_spatial"] = current_xy.copy()
        model.obsm["spatial"] = current_xy.copy()
        if "_in_alignment_overlap" in model.obs:
            del model.obs["_in_alignment_overlap"]
    return xenium_next, vhd_next


def compose_transformations(first, second):
    """Compose first-pass→intermediate and intermediate→final rigid maps."""
    if len(first) != 1 or len(second) != 1:
        raise ValueError("Expected one transformation for each pairwise round")
    rotation_1 = np.asarray(first[0]["Rotation"], dtype=float)
    translation_1 = np.asarray(first[0]["Translation"], dtype=float)
    rotation_2 = np.asarray(second[0]["Rotation"], dtype=float)
    translation_2 = np.asarray(second[0]["Translation"], dtype=float)
    rotation = rotation_2 @ rotation_1
    translation = translation_1 @ rotation_2.T + translation_2
    return [{"Rotation": rotation, "Translation": translation}]


def validate_composed_transformation(
    original_vhd: ad.AnnData,
    final_vhd: ad.AnnData,
    transformation,
) -> None:
    """Verify the saved composite PKL reproduces the final coordinates."""
    original_xy = np.asarray(
        original_vhd.obsm["spatial"], dtype=float
    )[:, :2]
    expected = (
        original_xy @ transformation[0]["Rotation"].T
        + transformation[0]["Translation"]
    )
    observed = np.asarray(
        final_vhd.obsm["align_spatial"], dtype=float
    )[:, :2]
    if not np.allclose(expected, observed, rtol=1e-5, atol=1e-4):
        maximum_error = float(np.max(np.abs(expected - observed)))
        raise RuntimeError(
            "Composite transformation does not reproduce final coordinates; "
            f"maximum absolute error={maximum_error}"
        )


def save_coordinates(
    vhd_aligned: ad.AnnData,
    xenium_z_height: float,
    output_path: Path,
    stage: str,
) -> None:
    """Save intermediate or final VHD coordinates in Module 02."""
    if stage not in {"intermediate", "second_pass"}:
        raise ValueError(f"Unsupported coordinate stage: {stage}")
    first_pass = pd.read_csv(VHD_FIRST_PASS_COORDINATES)
    coordinate_map = first_pass.set_index("vhd_obs_id")
    vhd_ids = pd.Index(vhd_aligned.obs_names.astype(str))
    matched = coordinate_map.loc[vhd_ids]
    aligned_xy = np.asarray(
        vhd_aligned.obsm["align_spatial"], dtype=float
    )[:, :2]
    output = {
        "donor": DONOR,
        "xenium_sample": XENIUM_SAMPLE,
        "vhd_sample": VHD_SAMPLE,
        "vhd_obs_id": vhd_ids,
        "x_original_um": matched["x_original_um"].to_numpy(),
        "y_original_um": matched["y_original_um"].to_numpy(),
        "x_alignment_input_um": matched[
            "x_alignment_input_um"
        ].to_numpy(),
        "y_alignment_input_flipped_um": matched[
            "y_alignment_input_flipped_um"
        ].to_numpy(),
        f"x_{stage}": aligned_xy[:, 0],
        f"y_{stage}": aligned_xy[:, 1],
        "xenium_z_height": xenium_z_height,
        "vhd_z_height": VHD_Z_HEIGHT,
        "Spatial_Domain": (
            vhd_aligned.obs["Spatial_Domain"]
            .astype("string")
            .fillna("")
            .astype(str)
            .to_numpy()
        ),
    }
    pd.DataFrame(output).to_csv(output_path, index=False)


def save_pre_alignment_label_overviews(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
) -> None:
    """Save the two separate first-pass-input inspection figures."""

    plot_label_subplots(
        xenium,
        spatial_key="plot_spatial",
        label_key="CellType",
        modality_name="Xenium",
        color=XENIUM_COLOR,
        output_path=(
            PLOT_DIR
            / f"{PAIR_TAG}_Xenium_by_CellType_first_pass_input_um.png"
        ),
    )
    plot_label_subplots(
        vhd,
        spatial_key="plot_spatial",
        label_key="Spatial_Domain",
        modality_name="VHD",
        color=VHD_COLOR,
        output_path=(
            PLOT_DIR
            / f"{PAIR_TAG}_VHD_by_Spatial_Domain_first_pass_input_um.png"
        ),
    )


def plot_aligned_stage_celltypes(
    xenium: ad.AnnData,
    vhd: ad.AnnData,
    by_celltype_dir: Path,
    *,
    filename_stage: str,
    display_stage: str,
) -> None:
    """Save the standard by-cell-type plots for one aligned stage."""
    for label in SHARED_LABEL_PLOTS:
        plot_by_label_overlay(
            xenium,
            vhd,
            display_name=label,
            xenium_labels=(label,),
            vhd_labels=(label,),
            output_path=(
                by_celltype_dir
                / f"{PAIR_TAG}_{label}_{filename_stage}_spateo_um.png"
            ),
            stage_label=display_stage,
        )

    plot_by_label_overlay(
        xenium,
        vhd,
        display_name="CHAT_PVALB_Inh",
        xenium_labels=("CHAT", "Inh_PVALB"),
        vhd_labels=("CHAT_PVALB",),
        output_path=(
            by_celltype_dir
            / f"{PAIR_TAG}_CHAT_PVALB_Inh_{filename_stage}_spateo_um.png"
        ),
        stage_label=display_stage,
    )
    plot_overlap_wm_overlay(
        xenium,
        vhd,
        spatial_key="align_spatial",
        stage_label=display_stage,
        output_path=(
            by_celltype_dir
            / f"{PAIR_TAG}_WM_overlap_{filename_stage}_spateo_um.png"
        ),
    )


def main() -> None:
    validate_input_paths()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    by_celltype_dir = PLOT_DIR / "by_celltype"
    by_celltype_dir.mkdir(parents=True, exist_ok=True)

    if DEVICE == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("CUDA requested, but CUDA is unavailable")

    print("Loading final aligned Xenium and first-pass VHD coordinates.")
    xenium, xenium_slices, xenium_z_height = load_xenium()
    vhd = load_vhd()
    print("Xenium cells:", f"{xenium.n_obs:,}")
    print("VHD observations:", f"{vhd.n_obs:,}")
    print("Xenium z height (um):", xenium_z_height)
    print("VHD z height (um):", VHD_Z_HEIGHT)

    print("Saving separate Xenium and VHD label overviews.")
    save_pre_alignment_label_overviews(xenium, vhd)
    plot_xenium_sanity(xenium_slices, XENIUM_SANITY_PATH)
    plot_pair_overlay(
        xenium,
        vhd,
        "plot_spatial",
        FIRST_PASS_OVERLAY_PATH,
        "First-pass input",
    )

    # Round 1: first-pass input -> intermediate.
    xenium_sub, vhd_sub, intermediate_qc = prepare_alignment_subsets(
        xenium,
        vhd,
        overlap_mode="tissue_then_WM",
        overlap_distance_um=INTERMEDIATE_OVERLAP_DISTANCE_UM,
    )
    intermediate_qc.to_csv(INTERMEDIATE_CELLTYPE_QC_PATH, index=False)
    print("Intermediate-round QC:")
    print(intermediate_qc.to_string(index=False))
    plot_overlap_wm_overlay(
        xenium,
        vhd,
        spatial_key="plot_spatial",
        stage_label="First-pass input",
        output_path=(
            by_celltype_dir
            / f"{PAIR_TAG}_WM_overlap_first_pass_input_um.png"
        ),
    )
    intermediate_transformation = fit_residual_transformation(
        xenium_sub,
        vhd_sub,
        stage_name="intermediate",
        overlap_distance_um=INTERMEDIATE_OVERLAP_DISTANCE_UM,
    )
    with INTERMEDIATE_TRANSFORMATION_PATH.open("wb") as handle:
        pickle.dump(
            intermediate_transformation,
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    xenium_intermediate, vhd_intermediate = apply_transformation(
        xenium,
        vhd,
        intermediate_transformation,
    )
    save_coordinates(
        vhd_intermediate,
        xenium_z_height,
        INTERMEDIATE_COORDINATE_PATH,
        stage="intermediate",
    )
    plot_pair_overlay(
        xenium_intermediate,
        vhd_intermediate,
        "align_spatial",
        INTERMEDIATE_OVERLAY_PATH,
        "Intermediate aligned",
    )
    plot_aligned_stage_celltypes(
        xenium_intermediate,
        vhd_intermediate,
        by_celltype_dir,
        filename_stage="intermediate",
        display_stage="Intermediate aligned",
    )

    intermediate_metadata = {
        "donor": DONOR,
        "stage": "intermediate",
        "input_coordinate_space": "module-01 first-pass coordinates",
        "overlap_definition": {
            "mode": "tissue_then_WM",
            "method": "nearest observation in the other modality",
            "maximum_distance_um": INTERMEDIATE_OVERLAP_DISTANCE_UM,
        },
        "realignment_label": REALIGN_LABEL,
        "realign_counts": intermediate_qc.to_dict(orient="records"),
        "model_order": [XENIUM_SAMPLE, VHD_SAMPLE],
        "fixed_model_unchanged": True,
        "device": DEVICE,
        "max_iter": MAX_ITER,
        "partial_robust_level": PARTIAL_ROBUST_LEVEL,
        "transformation_output": str(
            INTERMEDIATE_TRANSFORMATION_PATH.resolve()
        ),
        "coordinate_output": str(
            INTERMEDIATE_COORDINATE_PATH.resolve()
        ),
        "plot_directory": str(PLOT_DIR.resolve()),
    }
    INTERMEDIATE_METADATA_PATH.write_text(
        json.dumps(intermediate_metadata, indent=2) + "\n"
    )

    # Round 2: recompute overlap from intermediate coordinates, then refine.
    xenium_round_2, vhd_round_2 = prepare_next_round_inputs(
        xenium_intermediate,
        vhd_intermediate,
    )
    xenium_sub_2, vhd_sub_2, final_qc = prepare_alignment_subsets(
        xenium_round_2,
        vhd_round_2,
        overlap_mode="WM_to_WM",
        overlap_distance_um=FINAL_OVERLAP_DISTANCE_UM,
    )
    final_qc.to_csv(CELLTYPE_QC_PATH, index=False)
    print("Final-round QC:")
    print(final_qc.to_string(index=False))
    plot_overlap_wm_overlay(
        xenium_round_2,
        vhd_round_2,
        spatial_key="plot_spatial",
        stage_label="Intermediate input to final round",
        output_path=(
            by_celltype_dir
            / f"{PAIR_TAG}_WM_overlap_intermediate_input_um.png"
        ),
    )
    final_residual_transformation = fit_residual_transformation(
        xenium_sub_2,
        vhd_sub_2,
        stage_name="final second-pass",
        overlap_distance_um=FINAL_OVERLAP_DISTANCE_UM,
    )
    xenium_final, vhd_final = apply_transformation(
        xenium_round_2,
        vhd_round_2,
        final_residual_transformation,
    )

    final_composite_transformation = compose_transformations(
        intermediate_transformation,
        final_residual_transformation,
    )
    validate_composed_transformation(
        vhd,
        vhd_final,
        final_composite_transformation,
    )
    with TRANSFORMATION_PATH.open("wb") as handle:
        pickle.dump(
            final_composite_transformation,
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )
    save_coordinates(
        vhd_final,
        xenium_z_height,
        COORDINATE_PATH,
        stage="second_pass",
    )
    plot_pair_overlay(
        xenium_final,
        vhd_final,
        "align_spatial",
        SECOND_PASS_OVERLAY_PATH,
        "Final second-pass aligned",
    )
    plot_aligned_stage_celltypes(
        xenium_final,
        vhd_final,
        by_celltype_dir,
        filename_stage="second",
        display_stage="Final second-pass aligned",
    )

    metadata = {
        "donor": DONOR,
        "fixed_reference": {
            "technology": "Xenium",
            "sample": XENIUM_SAMPLE,
            "xenium_z_height": xenium_z_height,
            "coordinate_source": str(XENIUM_COORDINATES.resolve()),
        },
        "moving_model": {
            "technology": "VisiumHD",
            "sample": VHD_SAMPLE,
            "vhd_z_height": VHD_Z_HEIGHT,
            "coordinate_source": str(
                VHD_FIRST_PASS_COORDINATES.resolve()
            ),
        },
        "number_of_internal_residual_rounds": 2,
        "overlap_definitions": {
            "intermediate_round": {
                "coordinate_space": "module-01 first-pass coordinates",
                "mode": "tissue_then_WM",
                "method": (
                    "nearest observation in the other modality, then WM"
                ),
                "maximum_distance_um": INTERMEDIATE_OVERLAP_DISTANCE_UM,
            },
            "final_round": {
                "coordinate_space": "intermediate coordinates",
                "mode": "WM_to_WM",
                "method": (
                    "filter WM first, then nearest WM in the other modality"
                ),
                "maximum_distance_um": FINAL_OVERLAP_DISTANCE_UM,
            },
        },
        "realignment_labels": {
            "xenium_CellType": REALIGN_LABEL,
            "vhd_Spatial_Domain": REALIGN_LABEL,
        },
        "intermediate_round_counts": intermediate_qc.to_dict(
            orient="records"
        ),
        "final_round_counts": final_qc.to_dict(orient="records"),
        "model_order": [XENIUM_SAMPLE, VHD_SAMPLE],
        "fixed_model_unchanged": True,
        "intermediate_coordinate_input": str(
            INTERMEDIATE_COORDINATE_PATH.resolve()
        ),
        "intermediate_transformation_output": str(
            INTERMEDIATE_TRANSFORMATION_PATH.resolve()
        ),
        "final_residual_composed_with_intermediate": True,
        "saved_second_pass_transformation_scope": (
            "module-01 first-pass coordinates to final second-pass coordinates"
        ),
        "device": DEVICE,
        "max_iter": MAX_ITER,
        "partial_robust_level": PARTIAL_ROBUST_LEVEL,
        "transformation_output": str(TRANSFORMATION_PATH.resolve()),
        "coordinate_output": str(COORDINATE_PATH.resolve()),
        "plot_directory": str(PLOT_DIR.resolve()),
    }
    METADATA_PATH.write_text(json.dumps(metadata, indent=2) + "\n")

    print("Saved intermediate transformation:", INTERMEDIATE_TRANSFORMATION_PATH)
    print("Saved intermediate coordinates:", INTERMEDIATE_COORDINATE_PATH)
    print("Saved final composite transformation:", TRANSFORMATION_PATH)
    print("Saved final second-pass coordinates:", COORDINATE_PATH)
    print("Saved plots:", PLOT_DIR)
    print("Done.")


if __name__ == "__main__":
    main()


# Run directly:
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
# python 05_Spateo_Xenium_VHD_second_pass_Br6660_NAc6_3080-VHD_H1_M3TCP9V_A1.py
