#!/usr/bin/env python3
"""First-pass Spateo alignment of one VHD slice to one aligned Xenium slice.

Xenium is always the fixed reference (models[0]); VHD is always the moving
slice (models[1]). The fitted rigid transformation and all transformed VHD
coordinates are saved for downstream joint Xenium/VHD 3D reconstruction.
"""

from __future__ import annotations

import argparse
import importlib.metadata as md
import json
import os
import pickle
import re
import subprocess
import warnings
from pathlib import Path
from typing import Sequence

os.environ["PYVISTA_OFF_SCREEN"] = "true"
os.environ.setdefault("PYVISTA_EGL", "true")
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"
os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib")

warnings.filterwarnings(
    "ignore",
    message="The pynvml package is deprecated.*",
    category=FutureWarning,
)
warnings.filterwarnings("ignore", message="pkg_resources is deprecated as an API")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import anndata as ad
import numpy as np
import pandas as pd
import pyvista as pv
import scanpy as sc
import scipy.sparse as sp
import spateo as st
import torch
from numba.core.errors import NumbaWarning

warnings.filterwarnings("ignore", category=NumbaWarning)

try:
    pv.start_xvfb()
except Exception:
    pass

try:
    st.__version__ = md.version("spateo-release")
except md.PackageNotFoundError:
    pass


D1_CELL_TYPES = ("D1_Island_A", "D1_Island_B")
REFERENCE_COLOR = "#4C78A8"
MOVING_COLOR = "#F58518"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Align one moving VisiumHD slice to one fixed, final-aligned "
            "Xenium slice using Spateo."
        )
    )
    parser.add_argument("--donor", required=True, help="Xenium donor, e.g. Br6660")
    parser.add_argument(
        "--xenium-sample",
        required=True,
        help="Exact Xenium Sample value, e.g. Br6660_Nac10_4080",
    )
    parser.add_argument(
        "--vhd-sample",
        required=True,
        help=(
            "VHD sample label, e.g. VHD_H1_XKYDCP3_D1 or H1-XKYDCP3_D1"
        ),
    )
    parser.add_argument(
        "--vhd-z-height",
        required=True,
        type=float,
        help="Physical depth of the VHD section in microns, e.g. 4120.",
    )
    parser.add_argument(
        "--xenium-h5ad",
        type=Path,
        help="Override the normalized Xenium H5AD input.",
    )
    parser.add_argument(
        "--celltype-csv",
        type=Path,
        help="Override Banksy_cell_types.csv.",
    )
    parser.add_argument(
        "--aligned-coordinates",
        type=Path,
        help="Override the donor final second-pass coordinate CSV.",
    )
    parser.add_argument(
        "--vhd-h5ad",
        type=Path,
        help="Override the per-sample micron-coordinate VHD H5AD.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Override the module processed-data base directory; a pair-named "
            "subfolder is added automatically."
        ),
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help=(
            "Override the module plot base directory; a pair-named subfolder "
            "is added automatically."
        ),
    )
    parser.add_argument(
        "--device",
        choices=("auto", "cpu", "cuda"),
        default="auto",
        help="Spateo compute device (default: auto).",
    )
    parser.add_argument("--max-iter", type=int, default=500)
    parser.add_argument("--partial-robust-level", type=int, default=100)
    parser.add_argument("--chunk-capacity", type=int, default=2)
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Allow existing pair-specific PKL/CSV outputs to be replaced.",
    )
    return parser.parse_args()


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def safe_filename(value: str) -> str:
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", str(value)).strip("._")
    if not safe:
        raise ValueError(f"Cannot create a safe filename from {value!r}")
    return safe


def canonical_vhd_label(value: str) -> str:
    label = str(value).strip().replace("-", "_")
    if not label.startswith("VHD_"):
        label = f"VHD_{label}"
    return label


def resolve_paths(args: argparse.Namespace, root: Path) -> None:
    vhd_label = canonical_vhd_label(args.vhd_sample)
    args.vhd_label = vhd_label
    args.pair_tag = (
        f"{safe_filename(args.xenium_sample)}-{safe_filename(vhd_label)}"
    )
    args.xenium_h5ad = args.xenium_h5ad or (
        root
        / "processed-data"
        / "02_build_spe"
        / "h5ad"
        / "spe_NormCounts_nucleus_normcounts.h5ad"
    )
    args.celltype_csv = args.celltype_csv or (
        root
        / "processed-data"
        / "05_Clustering"
        / "Banksy_cell_types.csv"
    )
    args.aligned_coordinates = args.aligned_coordinates or (
        root
        / "processed-data"
        / "05_Xenium_alignment"
        / f"module02_Spateo_second_pass_{args.donor}"
        / f"{args.donor}_second_pass_coordinates.csv"
    )
    args.vhd_h5ad = args.vhd_h5ad or (
        root
        / "processed-data"
        / "HD_Full_Analysis"
        / "h5ad"
        / f"{vhd_label}_spaceranger_square008um_spatial_um.h5ad"
    )
    module_name = "module_01_Spateo_first_pass"
    output_base_dir = args.output_dir or (
        root / "processed-data" / "11_Xenium_VisiumHD_alignment" / module_name
    )
    plot_base_dir = args.plot_dir or (
        root / "plots" / "11_Xenium_VisiumHD_alignment" / module_name
    )
    args.output_dir = output_base_dir / args.pair_tag
    args.plot_dir = plot_base_dir / args.pair_tag
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)

    for path_name in (
        "xenium_h5ad",
        "celltype_csv",
        "aligned_coordinates",
        "vhd_h5ad",
    ):
        path = Path(getattr(args, path_name))
        if not path.is_file():
            raise FileNotFoundError(f"{path_name} does not exist: {path}")

    if not args.xenium_sample.startswith(f"{args.donor}_"):
        raise ValueError(
            f"--xenium-sample {args.xenium_sample!r} does not start with "
            f"the donor prefix {args.donor + '_'!r}"
        )
    if not np.isfinite(args.vhd_z_height) or args.vhd_z_height < 0:
        raise ValueError("--vhd-z-height must be a finite, non-negative value")


def prepare_output_paths(args: argparse.Namespace) -> dict[str, Path]:
    pair_tag = args.pair_tag
    paths = {
        "transformation": args.output_dir / f"{pair_tag}_first_pass.pkl",
        "coordinates": (
            args.output_dir / f"{pair_tag}_first_pass_coordinates.csv"
        ),
        "metadata": args.output_dir / f"{pair_tag}_first_pass_metadata.json",
        "gene_qc": args.output_dir / f"{pair_tag}_gene_overlap_qc.csv",
        "xenium_sanity": (
            args.plot_dir / f"{safe_filename(args.donor)}_spateo_Xenium_sanity_check.png"
        ),
        "before": args.plot_dir / f"{pair_tag}_before_spateo_um.png",
        "pca": args.plot_dir / f"{pair_tag}_pca.png",
        "initial": args.plot_dir / f"{pair_tag}_initial_spateo_um.png",
    }
    protected = [paths["transformation"], paths["coordinates"]]
    existing = [path for path in protected if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Pair-specific output(s) already exist; rerun with --overwrite: "
            + ", ".join(map(str, existing))
        )
    return paths


def load_xenium(
    args: argparse.Namespace,
) -> tuple[ad.AnnData, list[ad.AnnData], pd.DataFrame]:
    adata_all = ad.read_h5ad(args.xenium_h5ad)
    required_obs = {"Donor", "Sample", "cell_id"}
    missing_obs = required_obs - set(adata_all.obs.columns)
    if missing_obs:
        raise ValueError(f"Xenium H5AD is missing obs columns: {sorted(missing_obs)}")

    annotations = pd.read_csv(
        args.celltype_csv, dtype={"cell_id": "string"}
    )
    if not {"cell_id", "CellType"}.issubset(annotations.columns):
        raise ValueError(
            "Cell-type CSV must contain columns 'cell_id' and 'CellType'"
        )
    annotation_ids = annotations["cell_id"].astype("string")
    if annotation_ids.isna().any() or annotation_ids.duplicated().any():
        raise ValueError("Banksy cell-type CSV cell_id values must be unique")
    ct_map = annotations.set_index("cell_id")["CellType"]

    adata_all.obs["cell_id"] = adata_all.obs["cell_id"].astype("string")
    adata_all.obs["CellType"] = adata_all.obs["cell_id"].map(ct_map)
    donor_mask = adata_all.obs["Donor"].astype(str).eq(args.donor)
    donor_adata = adata_all[donor_mask].copy()
    if donor_adata.n_obs == 0:
        raise ValueError(f"No Xenium cells found for donor {args.donor!r}")
    donor_adata.obs = donor_adata.obs.reset_index(drop=True)
    donor_adata.obsm["spatial"] = np.asarray(
        getattr(donor_adata.obsm["spatial"], "values", donor_adata.obsm["spatial"])
    )

    coordinate_columns = [
        "donor",
        "sample",
        "cell_id",
        "x_second_pass",
        "y_second_pass",
        "z_height",
    ]
    coordinates = pd.read_csv(
        args.aligned_coordinates,
        usecols=coordinate_columns,
        dtype={"cell_id": "string", "sample": "string", "donor": "string"},
    )
    if coordinates["cell_id"].isna().any():
        raise ValueError("Final coordinate CSV contains missing cell IDs")
    if coordinates["cell_id"].duplicated().any():
        raise ValueError("Final coordinate CSV contains duplicate cell IDs")
    observed_donors = set(coordinates["donor"].dropna().astype(str))
    if observed_donors != {args.donor}:
        raise ValueError(
            f"Coordinate CSV donor mismatch: expected {args.donor}, "
            f"observed {sorted(observed_donors)}"
        )
    if coordinates[["x_second_pass", "y_second_pass", "z_height"]].isna().any().any():
        raise ValueError("Final coordinate CSV contains missing coordinates")

    coordinate_map = coordinates.set_index("cell_id")
    donor_ids = donor_adata.obs["cell_id"].astype("string")
    missing_ids = donor_ids[~donor_ids.isin(coordinate_map.index)]
    if len(missing_ids):
        raise ValueError(
            f"{len(missing_ids)} donor Xenium cells are missing from final coordinates"
        )
    matched = coordinate_map.loc[donor_ids]
    donor_adata.obsm["align_spatial"] = matched[
        ["x_second_pass", "y_second_pass"]
    ].to_numpy(dtype=float)
    donor_adata.obs["z_height"] = matched["z_height"].to_numpy(dtype=float)

    sample_order = (
        coordinates.groupby("sample", observed=True)["z_height"]
        .median()
        .sort_values()
        .index.astype(str)
        .tolist()
    )
    observed_samples = set(donor_adata.obs["Sample"].astype(str))
    if set(sample_order) != observed_samples:
        raise ValueError(
            "Sample mismatch between Xenium H5AD and final coordinate CSV"
        )
    ordered_slices = [
        donor_adata[
            donor_adata.obs["Sample"].astype(str).eq(sample)
        ].copy()
        for sample in sample_order
    ]

    selected = donor_adata[
        donor_adata.obs["Sample"].astype(str).eq(args.xenium_sample)
    ].copy()
    if selected.n_obs == 0:
        raise ValueError(
            f"Xenium sample {args.xenium_sample!r} is absent; "
            f"available samples: {sample_order}"
        )
    selected_z = coordinates.loc[
        coordinates["sample"].astype(str).eq(args.xenium_sample), "z_height"
    ].unique()
    if len(selected_z) != 1:
        raise ValueError(
            f"Expected one z_height for {args.xenium_sample}, "
            f"observed {selected_z.tolist()}"
        )
    args.xenium_z_height = float(selected_z[0])
    selected.obsm["spatial"] = np.asarray(
        selected.obsm["align_spatial"], dtype=float
    ).copy()
    selected.obs["sample_id"] = args.xenium_sample
    return selected, ordered_slices, coordinates


def load_vhd(args: argparse.Namespace) -> ad.AnnData:
    vhd = ad.read_h5ad(args.vhd_h5ad)
    if "spatial" not in vhd.obsm:
        raise ValueError("VHD H5AD does not contain obsm['spatial']")
    if "Spatial_Domain" not in vhd.obs:
        raise ValueError("VHD H5AD does not contain obs['Spatial_Domain']")

    xy_original = np.asarray(vhd.obsm["spatial"], dtype=float)[:, :2].copy()
    vhd.obsm["spatial_original_um"] = xy_original.copy()
    xy = xy_original.copy()
    y_min, y_max = float(xy[:, 1].min()), float(xy[:, 1].max())
    xy[:, 1] = y_min + y_max - xy[:, 1]
    vhd.obsm["spatial"] = xy
    vhd.obs["CellType"] = vhd.obs["Spatial_Domain"].astype("string")
    vhd.obs["sample_id"] = args.vhd_label
    return vhd


def display_name(sample: str, donor: str) -> str:
    prefix = f"{donor}_"
    return sample[len(prefix):] if sample.startswith(prefix) else sample


def subplot_grid(n_panels: int) -> tuple[plt.Figure, np.ndarray]:
    if n_panels < 1:
        raise ValueError("At least one panel is required")
    n_columns = min(3, n_panels)
    n_rows = int(np.ceil(n_panels / n_columns))
    return plt.subplots(
        n_rows,
        n_columns,
        figsize=(5.0 * n_columns, 4.5 * n_rows),
        squeeze=False,
    )


def plot_consecutive_xenium_overlays(
    slices: Sequence[ad.AnnData],
    donor: str,
    output_path: Path,
) -> None:
    figure, axes = subplot_grid(len(slices) - 1)
    axes_flat = axes.ravel()
    names = [
        display_name(str(model.obs["Sample"].astype(str).iloc[0]), donor)
        for model in slices
    ]
    for index, axis in enumerate(axes_flat[: len(slices) - 1]):
        reference_xy = np.asarray(slices[index].obsm["align_spatial"])[:, :2]
        moving_xy = np.asarray(slices[index + 1].obsm["align_spatial"])[:, :2]
        axis.scatter(
            reference_xy[:, 0],
            reference_xy[:, 1],
            s=0.2,
            alpha=0.45,
            color=REFERENCE_COLOR,
            linewidths=0,
            label=names[index],
        )
        axis.scatter(
            moving_xy[:, 0],
            moving_xy[:, 1],
            s=0.2,
            alpha=0.45,
            color=MOVING_COLOR,
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
            loc="upper right", fontsize=6, markerscale=8, frameon=False
        )
    for axis in axes_flat[len(slices) - 1 :]:
        axis.axis("off")
    figure.suptitle(
        "Final second-pass Xenium consecutive-slice overlays", y=1.002
    )
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def plot_pair_overlay(
    xen: ad.AnnData,
    vhd: ad.AnnData,
    spatial_key: str,
    xenium_name: str,
    vhd_name: str,
    output_path: Path,
    stage_label: str,
    cell_type: str | None = None,
) -> None:
    xen_xy = np.asarray(xen.obsm[spatial_key], dtype=float)[:, :2]
    vhd_xy = np.asarray(vhd.obsm[spatial_key], dtype=float)[:, :2]
    if cell_type is not None:
        xen_mask = (
            xen.obs["CellType"]
            .astype("string")
            .eq(cell_type)
            .fillna(False)
            .to_numpy(dtype=bool)
        )
        vhd_mask = (
            vhd.obs["Spatial_Domain"]
            .astype("string")
            .eq(cell_type)
            .fillna(False)
            .to_numpy(dtype=bool)
        )
        xen_xy = xen_xy[xen_mask]
        vhd_xy = vhd_xy[vhd_mask]
        print(
            f"{cell_type}: Xenium cells={len(xen_xy):,}; "
            f"VHD observations={len(vhd_xy):,}"
        )

    point_size = 0.6 if cell_type is not None else 0.2
    point_alpha = 0.8 if cell_type is not None else 0.45
    figure, axis = plt.subplots(figsize=(5.0, 4.5))
    axis.scatter(
        xen_xy[:, 0],
        xen_xy[:, 1],
        s=point_size,
        alpha=point_alpha,
        color=REFERENCE_COLOR,
        linewidths=0,
        label=xenium_name,
    )
    axis.scatter(
        vhd_xy[:, 0],
        vhd_xy[:, 1],
        s=point_size,
        alpha=point_alpha,
        color=MOVING_COLOR,
        linewidths=0,
        label=vhd_name,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    title = f"{stage_label}: {xenium_name} → {vhd_name}"
    if cell_type is not None:
        title = f"{cell_type}: {title}"
    axis.set_title(title, fontsize=9)
    axis.legend(loc="upper right", fontsize=6, markerscale=8, frameon=False)
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def prepare_alignment_inputs(
    xen_full: ad.AnnData, vhd_full: ad.AnnData
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    xen = xen_full.copy()
    vhd = vhd_full.copy()
    xen.var_names_make_unique()
    vhd.var_names_make_unique()

    vhd_genes_before = pd.Index(vhd.var_names)
    xen_genes = pd.Index(xen.var_names)
    common_before = xen_genes.intersection(vhd_genes_before)

    sc.pp.filter_cells(vhd, min_genes=10)
    vhd.layers["counts_full_before_intersection"] = vhd.X.copy()
    if sp.issparse(vhd.X):
        library_size = np.asarray(vhd.X.sum(axis=1)).ravel()
    else:
        library_size = np.asarray(vhd.X.sum(axis=1)).ravel()
    vhd.obs["lib_size_before_intersection"] = library_size
    sc.pp.filter_genes(vhd, min_cells=3)

    vhd_genes_after = pd.Index(vhd.var_names)
    common_genes = xen_genes.intersection(vhd_genes_after)
    if len(common_genes) < 2:
        raise ValueError(
            f"Only {len(common_genes)} common genes remain after VHD filtering"
        )

    qc = pd.DataFrame(
        [
            {
                "xenium_genes": len(xen_genes),
                "vhd_genes_before_filter": len(vhd_genes_before),
                "common_genes_before_filter": len(common_before),
                "vhd_genes_after_filter": len(vhd_genes_after),
                "common_genes_used": len(common_genes),
                "xenium_genes_missing_in_vhd_after_filter": len(
                    xen_genes.difference(vhd_genes_after)
                ),
            }
        ]
    )

    xen = xen[:, common_genes].copy()
    vhd = vhd[:, common_genes].copy()
    xen.layers["lognorm"] = xen.X.copy()
    vhd.layers["counts"] = vhd.X.copy()

    target_sum = 1e4
    library_size = vhd.obs["lib_size_before_intersection"].to_numpy(dtype=float)
    library_size[library_size == 0] = np.nan
    if sp.issparse(vhd.layers["counts"]):
        normalized = vhd.layers["counts"].multiply(
            target_sum / library_size[:, None]
        ).tocsr()
        normalized.data = np.log1p(normalized.data)
    else:
        normalized = np.log1p(
            vhd.layers["counts"] * (target_sum / library_size[:, None])
        )
    vhd.X = normalized
    vhd.layers["lognorm"] = vhd.X.copy()

    # PCA order matches the alignment order: fixed Xenium, then moving VHD.
    st.align.group_pca([xen, vhd], pca_key="X_pca")
    for model in (xen, vhd):
        if "align_spatial" in model.obsm:
            del model.obsm["align_spatial"]
    return xen, vhd, qc


def plot_joint_pca(
    xen: ad.AnnData,
    vhd: ad.AnnData,
    xenium_name: str,
    vhd_name: str,
    output_path: Path,
) -> None:
    joint_pca = np.vstack([xen.obsm["X_pca"], vhd.obsm["X_pca"]])
    variance = np.var(joint_pca, axis=0, ddof=1)
    variance_ratio = variance / variance.sum()
    figure, axis = plt.subplots(figsize=(6, 5))
    for model, name, color in (
        (xen, xenium_name, REFERENCE_COLOR),
        (vhd, vhd_name, MOVING_COLOR),
    ):
        pca = np.asarray(model.obsm["X_pca"])
        axis.scatter(
            pca[:, 0],
            pca[:, 1],
            s=1,
            alpha=0.4,
            label=name,
            color=color,
            linewidths=0,
        )
    axis.set_xlabel(f"PC1 ({variance_ratio[0] * 100:.1f}% var.)")
    axis.set_ylabel(f"PC2 ({variance_ratio[1] * 100:.1f}% var.)")
    axis.set_title("Joint PCA: fixed Xenium and moving VHD")
    axis.legend(loc="upper right", markerscale=5, frameon=False)
    figure.tight_layout()
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def save_vhd_coordinates(
    vhd: ad.AnnData,
    args: argparse.Namespace,
    output_path: Path,
) -> None:
    original = np.asarray(vhd.obsm["spatial_original_um"], dtype=float)[:, :2]
    alignment_input = np.asarray(vhd.obsm["spatial"], dtype=float)[:, :2]
    aligned = np.asarray(vhd.obsm["align_spatial"], dtype=float)[:, :2]
    output = pd.DataFrame(
        {
            "donor": args.donor,
            "xenium_sample": args.xenium_sample,
            "vhd_sample": args.vhd_label,
            "vhd_obs_id": vhd.obs_names.astype(str),
            "x_original_um": original[:, 0],
            "y_original_um": original[:, 1],
            "x_alignment_input_um": alignment_input[:, 0],
            "y_alignment_input_flipped_um": alignment_input[:, 1],
            "x_first_pass": aligned[:, 0],
            "y_first_pass": aligned[:, 1],
            "xenium_z_height": args.xenium_z_height,
            "vhd_z_height": args.vhd_z_height,
            "Spatial_Domain": (
                vhd.obs["Spatial_Domain"]
                .astype("string")
                .fillna("")
                .astype(str)
                .to_numpy()
            ),
        }
    )
    if output["vhd_obs_id"].duplicated().any():
        raise ValueError("VHD observation IDs are not unique")
    if output[["x_first_pass", "y_first_pass"]].isna().any().any():
        raise ValueError("Aligned VHD coordinates contain missing values")
    output.to_csv(output_path, index=False)


def main() -> None:
    args = parse_args()
    root = git_root()
    resolve_paths(args, root)
    paths = prepare_output_paths(args)
    device = (
        ("cuda" if torch.cuda.is_available() else "cpu")
        if args.device == "auto"
        else args.device
    )
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("--device cuda was requested, but CUDA is unavailable")

    print("Spateo version:", st.__version__)
    print("Device:", device)
    print("Donor:", args.donor)
    print("Fixed Xenium reference:", args.xenium_sample)
    print("Moving VHD slice:", args.vhd_label)
    print("VHD z height (um):", args.vhd_z_height)
    print("Final Xenium coordinates:", args.aligned_coordinates)
    print("VHD H5AD:", args.vhd_h5ad)

    xen_full, xenium_slices, coordinate_table = load_xenium(args)
    vhd_full = load_vhd(args)
    print("Xenium z height (um):", args.xenium_z_height)

    plt.ioff()
    plot_consecutive_xenium_overlays(
        xenium_slices, args.donor, paths["xenium_sanity"]
    )
    plot_pair_overlay(
        xen_full,
        vhd_full,
        "spatial",
        args.xenium_sample,
        args.vhd_label,
        paths["before"],
        stage_label="Before Spateo",
    )

    xen_align, vhd_align, gene_qc = prepare_alignment_inputs(
        xen_full, vhd_full
    )
    gene_qc.to_csv(paths["gene_qc"], index=False)
    plot_joint_pca(
        xen_align,
        vhd_align,
        args.xenium_sample,
        args.vhd_label,
        paths["pca"],
    )

    spatial_key = "spatial"
    key_added = "align_spatial"
    # The order is intentional and must be identical during fit and apply:
    # Xenium is fixed (models[0]); VHD is transformed (models[1]).
    transformation = st.align.morpho_align_transformation(
        models=[xen_align, vhd_align],
        spatial_key=spatial_key,
        key_added=key_added,
        device=device,
        verbose=True,
        rep_layer="X_pca",
        rep_field="obsm",
        dissimilarity="cos",
        max_iter=args.max_iter,
        partial_robust_level=args.partial_robust_level,
        sparse_calculation_mode=True,
        use_chunk=True,
        chunk_capacity=args.chunk_capacity,
    )

    with paths["transformation"].open("wb") as handle:
        pickle.dump(transformation, handle, protocol=pickle.HIGHEST_PROTOCOL)

    # Apply the fitted transformation to the complete, unfiltered VHD object so
    # downstream reconstruction receives coordinates for every VHD observation.
    for model in (xen_full, vhd_full):
        if key_added in model.obsm:
            del model.obsm[key_added]
    aligned_slices = st.align.morpho_align_apply_transformation(
        models=[xen_full, vhd_full],
        spatial_key=spatial_key,
        key_added=key_added,
        transformation=transformation,
    )
    aligned_xen, aligned_vhd = aligned_slices

    if not np.allclose(
        np.asarray(aligned_xen.obsm[key_added]),
        np.asarray(aligned_xen.obsm[spatial_key]),
    ):
        raise RuntimeError("Fixed Xenium coordinates changed unexpectedly")

    save_vhd_coordinates(aligned_vhd, args, paths["coordinates"])
    plot_pair_overlay(
        aligned_xen,
        aligned_vhd,
        key_added,
        args.xenium_sample,
        args.vhd_label,
        paths["initial"],
        stage_label="Initial Spateo aligned",
    )

    by_celltype_dir = args.plot_dir / "by_celltype"
    for cell_type in D1_CELL_TYPES:
        plot_pair_overlay(
            aligned_xen,
            aligned_vhd,
            key_added,
            args.xenium_sample,
            args.vhd_label,
            by_celltype_dir
            / (
                f"{safe_filename(args.xenium_sample)}-"
                f"{safe_filename(args.vhd_label)}_"
                f"{safe_filename(cell_type)}_initial_spateo_um.png"
            ),
            stage_label="Initial Spateo aligned",
            cell_type=cell_type,
        )

    metadata = {
        "donor": args.donor,
        "fixed_reference": {
            "technology": "Xenium",
            "sample": args.xenium_sample,
            "coordinate_source": str(args.aligned_coordinates.resolve()),
            "coordinate_columns": ["x_second_pass", "y_second_pass"],
            "xenium_z_height": args.xenium_z_height,
        },
        "moving_model": {
            "technology": "VisiumHD",
            "sample": args.vhd_label,
            "vhd_z_height": args.vhd_z_height,
            "input": str(args.vhd_h5ad.resolve()),
            "aligned_coordinate_output": str(paths["coordinates"].resolve()),
        },
        "model_order_during_fit": [args.xenium_sample, args.vhd_label],
        "model_order_during_apply": [args.xenium_sample, args.vhd_label],
        "transformation_formula": (
            "vhd_align_spatial = vhd_spatial @ Rotation.T + Translation"
        ),
        "xenium_cells": int(aligned_xen.n_obs),
        "vhd_observations": int(aligned_vhd.n_obs),
        "xenium_donor_slices": int(len(xenium_slices)),
        "xenium_sample_order": [
            str(model.obs["Sample"].astype(str).iloc[0])
            for model in xenium_slices
        ],
        "common_genes_used": int(gene_qc.loc[0, "common_genes_used"]),
        "spateo_version": str(st.__version__),
        "device": device,
        "max_iter": args.max_iter,
        "partial_robust_level": args.partial_robust_level,
        "chunk_capacity": args.chunk_capacity,
        "transformation_output": str(paths["transformation"].resolve()),
        "plot_directory": str(args.plot_dir.resolve()),
    }
    paths["metadata"].write_text(json.dumps(metadata, indent=2) + "\n")

    print("Saved transformation:", paths["transformation"])
    print("Saved aligned VHD coordinates:", paths["coordinates"])
    print("Saved plots:", args.plot_dir)
    print("Done.")


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/11_Xenium_VisiumHD_alignment
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# python 04_Spateo_Xenium_VHD_first_pass.py \
#   --donor Br6660 \
#   --xenium-sample Br6660_Nac10_4080 \
#   --vhd-sample VHD_H1_XKYDCP3_D1 \
#   --vhd-z-height 4120 \
#   --overwrite
#
# Optional path/device overrides:
# python 04_Spateo_Xenium_VHD_first_pass.py \
#   --donor Br6660 \
#   --xenium-sample Br6660_Nac10_4080 \
#   --vhd-sample VHD_H1_XKYDCP3_D1 \
#   --vhd-z-height 4120 \
#   --aligned-coordinates /path/to/final_coordinates.csv \
#   --vhd-h5ad /path/to/vhd_spatial_um.h5ad \
#   --device cuda \
#   --overwrite
