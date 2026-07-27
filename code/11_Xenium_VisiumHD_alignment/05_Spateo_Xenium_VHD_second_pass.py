#!/usr/bin/env python3
"""Label-guided second-pass refinement of VHD-to-Xenium alignment.

Module 01 ``x_first_pass``/``y_first_pass`` VHD coordinates are the direct
input. A rigid transformation is fitted on selected labels in that coordinate
space and applied to the complete VHD slice; Xenium remains fixed throughout.
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import os
import pickle
import subprocess
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


def load_first_pass_helpers():
    """Load plotting/data helpers from the paired Module 01 CLI script."""
    helper_path = Path(__file__).with_name(
        "04_Spateo_Xenium_VHD_first_pass.py"
    )
    spec = importlib.util.spec_from_file_location(
        "_xenium_vhd_first_pass_helpers", helper_path
    )
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load Module 01 helpers from {helper_path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FIRST_PASS = load_first_pass_helpers()
SHARED_LABEL_PLOTS = (
    "D1_Island_A",
    "D1_Island_B",
    "Ependymal",
    "Excitatory",
    "WM",
)
CHAT_PVALB_INH_NAME = "CHAT_PVALB_Inh"
CHAT_PVALB_INH_XENIUM_LABELS = ("CHAT", "Inh_PVALB")
CHAT_PVALB_INH_VHD_LABELS = ("CHAT_PVALB",)
BY_CELLTYPE_XENIUM_POINT_SIZE = 0.6
BY_CELLTYPE_XENIUM_ALPHA = 0.60
BY_CELLTYPE_VHD_POINT_SIZE = 0.6
BY_CELLTYPE_VHD_ALPHA = 0.35


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Refine first-pass VHD coordinates in the fixed Xenium space using "
            "selected cell-type/spatial-domain labels."
        )
    )
    parser.add_argument("--donor", required=True)
    parser.add_argument("--xenium-sample", required=True)
    parser.add_argument(
        "--vhd-sample",
        required=True,
        help="VHD label, e.g. VHD_H1_XKYDCP3_D1.",
    )
    parser.add_argument(
        "--vhd-z-height",
        required=True,
        type=float,
        help="Physical depth of the VHD section in microns.",
    )
    parser.add_argument(
        "--realign-celltypes",
        nargs="+",
        required=True,
        help=(
            "One or more shared labels used for refinement. Xenium labels come "
            "from obs['CellType']; VHD labels come from obs['Spatial_Domain']. "
            "The special label CHAT_PVALB_Inh maps Xenium CHAT + Inh_PVALB "
            "to VHD CHAT_PVALB."
        ),
    )
    parser.add_argument("--xenium-h5ad", type=Path)
    parser.add_argument("--celltype-csv", type=Path)
    parser.add_argument("--aligned-coordinates", type=Path)
    parser.add_argument("--vhd-h5ad", type=Path)
    parser.add_argument(
        "--first-pass-coordinates",
        type=Path,
        help="Override the pair-specific Module 01 VHD coordinate CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        help=(
            "Override the Module 02 processed-data base directory; a "
            "pair-named subfolder is added automatically."
        ),
    )
    parser.add_argument(
        "--plot-dir",
        type=Path,
        help=(
            "Override the Module 02 plot base directory; a pair-named "
            "subfolder is added automatically."
        ),
    )
    parser.add_argument(
        "--device", choices=("auto", "cpu", "cuda"), default="auto"
    )
    parser.add_argument("--max-iter", type=int, default=200)
    parser.add_argument("--partial-robust-level", type=int, default=100)
    parser.add_argument("--overwrite", action="store_true")
    return parser.parse_args()


def git_root() -> Path:
    return Path(
        subprocess.check_output(
            ["git", "rev-parse", "--show-toplevel"], text=True
        ).strip()
    )


def resolve_paths(args: argparse.Namespace, root: Path) -> None:
    args.vhd_label = FIRST_PASS.canonical_vhd_label(args.vhd_sample)
    args.pair_tag = (
        f"{FIRST_PASS.safe_filename(args.xenium_sample)}-"
        f"{FIRST_PASS.safe_filename(args.vhd_label)}"
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
        / (
            f"{args.vhd_label}_"
            "spaceranger_square008um_spatial_um.h5ad"
        )
    )
    first_pass_pair_dir = (
        root
        / "processed-data"
        / "11_Xenium_VisiumHD_alignment"
        / "module_01_Spateo_first_pass"
        / args.pair_tag
    )
    args.first_pass_coordinates = (
        args.first_pass_coordinates
        or first_pass_pair_dir
        / f"{args.pair_tag}_first_pass_coordinates.csv"
    )

    module_name = "module_02_Spateo_second_pass"
    output_base = args.output_dir or (
        root
        / "processed-data"
        / "11_Xenium_VisiumHD_alignment"
        / module_name
    )
    plot_base = args.plot_dir or (
        root / "plots" / "11_Xenium_VisiumHD_alignment" / module_name
    )
    args.output_dir = output_base / args.pair_tag
    args.plot_dir = plot_base / args.pair_tag
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.plot_dir.mkdir(parents=True, exist_ok=True)

    for name in (
        "xenium_h5ad",
        "celltype_csv",
        "aligned_coordinates",
        "vhd_h5ad",
        "first_pass_coordinates",
    ):
        path = Path(getattr(args, name))
        if not path.is_file():
            raise FileNotFoundError(f"{name} does not exist: {path}")
    if not args.xenium_sample.startswith(f"{args.donor}_"):
        raise ValueError(
            f"--xenium-sample does not match donor {args.donor}: "
            f"{args.xenium_sample}"
        )
    if not np.isfinite(args.vhd_z_height) or args.vhd_z_height < 0:
        raise ValueError("--vhd-z-height must be finite and non-negative")

    # Remove duplicates while preserving the user's requested order.
    args.realign_celltypes = list(
        dict.fromkeys(map(str, args.realign_celltypes))
    )
    if any(not label.strip() for label in args.realign_celltypes):
        raise ValueError("--realign-celltypes cannot contain empty labels")


def output_paths(args: argparse.Namespace) -> dict[str, Path]:
    tag = args.pair_tag
    paths = {
        "transformation": args.output_dir / f"{tag}_second_pass.pkl",
        "coordinates": (
            args.output_dir / f"{tag}_second_pass_coordinates.csv"
        ),
        "metadata": args.output_dir / f"{tag}_second_pass_metadata.json",
        "celltype_qc": args.output_dir / f"{tag}_realign_celltype_qc.csv",
        "first_pass_plot": args.plot_dir / f"{tag}_first_pass_input_um.png",
        "second_pass_plot": args.plot_dir / f"{tag}_second_spateo_um.png",
        "xenium_sanity": (
            args.plot_dir
            / f"{FIRST_PASS.safe_filename(args.donor)}_spateo_Xenium_sanity_check.png"
        ),
    }
    protected = [
        paths["transformation"],
        paths["coordinates"],
    ]
    existing = [path for path in protected if path.exists()]
    if existing and not args.overwrite:
        raise FileExistsError(
            "Second-pass output(s) already exist; rerun with --overwrite: "
            + ", ".join(map(str, existing))
        )
    return paths


def load_first_pass_coordinates(
    xen: ad.AnnData,
    vhd: ad.AnnData,
    args: argparse.Namespace,
) -> tuple[ad.AnnData, ad.AnnData]:
    """Load Module 01 VHD first-pass x/y coordinates directly from CSV."""
    required = [
        "donor",
        "xenium_sample",
        "vhd_sample",
        "vhd_obs_id",
        "x_first_pass",
        "y_first_pass",
        "xenium_z_height",
        "vhd_z_height",
    ]
    coordinates = pd.read_csv(
        args.first_pass_coordinates,
        usecols=required,
        dtype={
            "donor": "string",
            "xenium_sample": "string",
            "vhd_sample": "string",
            "vhd_obs_id": "string",
        },
    )
    if coordinates["vhd_obs_id"].isna().any():
        raise ValueError("First-pass coordinate CSV contains missing VHD IDs")
    if coordinates["vhd_obs_id"].duplicated().any():
        raise ValueError("First-pass coordinate CSV contains duplicate VHD IDs")
    expected_values = {
        "donor": args.donor,
        "xenium_sample": args.xenium_sample,
        "vhd_sample": args.vhd_label,
    }
    for column, expected in expected_values.items():
        observed = set(coordinates[column].dropna().astype(str))
        if observed != {expected}:
            raise ValueError(
                f"First-pass CSV {column} mismatch: expected {expected!r}, "
                f"observed {sorted(observed)}"
            )
    if not np.allclose(
        coordinates["xenium_z_height"].to_numpy(dtype=float),
        args.xenium_z_height,
    ):
        raise ValueError("First-pass CSV xenium_z_height does not match Xenium")
    if not np.allclose(
        coordinates["vhd_z_height"].to_numpy(dtype=float),
        args.vhd_z_height,
    ):
        raise ValueError("First-pass CSV vhd_z_height does not match CLI")
    if coordinates[["x_first_pass", "y_first_pass"]].isna().any().any():
        raise ValueError("First-pass CSV contains missing aligned coordinates")

    coordinate_map = coordinates.set_index("vhd_obs_id")
    vhd_ids = pd.Index(vhd.obs_names.astype(str))
    missing = vhd_ids[~vhd_ids.isin(coordinate_map.index)]
    unexpected = coordinate_map.index[~coordinate_map.index.isin(vhd_ids)]
    if len(missing) or len(unexpected):
        raise ValueError(
            "VHD ID mismatch between H5AD and first-pass CSV: "
            f"missing={len(missing)}, unexpected={len(unexpected)}"
        )
    first_xy = coordinate_map.loc[
        vhd_ids, ["x_first_pass", "y_first_pass"]
    ].to_numpy(dtype=float)

    xen_first = xen.copy()
    vhd_first = vhd.copy()
    xen_first.obsm["align_spatial"] = np.asarray(
        xen_first.obsm["spatial"], dtype=float
    ).copy()
    vhd_first.obsm["align_spatial"] = first_xy.copy()
    # Second-pass fitting/application uses first-pass coordinates as its input
    # spatial key; original VHD spatial remains available in the separate
    # vhd object passed to the output writer.
    xen_first.obsm["spatial"] = np.asarray(
        xen_first.obsm["align_spatial"], dtype=float
    ).copy()
    vhd_first.obsm["spatial"] = first_xy.copy()
    return xen_first, vhd_first


def apply_transformation(
    xen: ad.AnnData,
    vhd: ad.AnnData,
    transformation: list[dict[str, np.ndarray]],
) -> tuple[ad.AnnData, ad.AnnData]:
    fixed = xen.copy()
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


def prepare_label_subsets(
    xen_first: ad.AnnData,
    vhd_first: ad.AnnData,
    celltypes: list[str],
) -> tuple[ad.AnnData, ad.AnnData, pd.DataFrame]:
    # Preserve the only expression-based filtering step that affected the old
    # label refinement. Gene intersection/PCA are unnecessary for label mode.
    vhd_fit = vhd_first.copy()
    sc.pp.filter_cells(vhd_fit, min_genes=10)

    xen_labels = xen_first.obs["CellType"].astype("string")
    vhd_labels = vhd_fit.obs["Spatial_Domain"].astype("string")
    xen_alignment_labels = xen_labels.mask(
        xen_labels.isin(CHAT_PVALB_INH_XENIUM_LABELS),
        CHAT_PVALB_INH_NAME,
    )
    vhd_alignment_labels = vhd_labels.mask(
        vhd_labels.isin(CHAT_PVALB_INH_VHD_LABELS),
        CHAT_PVALB_INH_NAME,
    )
    qc_rows = []
    for cell_type in celltypes:
        qc_rows.append(
            {
                "cell_type": cell_type,
                "n_xenium": int(
                    xen_alignment_labels.eq(cell_type).fillna(False).sum()
                ),
                "n_vhd": int(
                    vhd_alignment_labels.eq(cell_type).fillna(False).sum()
                ),
            }
        )
    qc = pd.DataFrame(qc_rows)
    missing = qc[(qc["n_xenium"] == 0) | (qc["n_vhd"] == 0)]
    if not missing.empty:
        raise ValueError(
            "Every realignment label must occur in both modalities:\n"
            + missing.to_string(index=False)
        )

    xen_mask = (
        xen_alignment_labels.isin(celltypes)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    vhd_mask = (
        vhd_alignment_labels.isin(celltypes)
        .fillna(False)
        .to_numpy(dtype=bool)
    )
    xen_sub = xen_first[xen_mask].copy()
    vhd_sub = vhd_fit[vhd_mask].copy()

    # Fit the second-pass transform directly in the first-pass coordinate space.
    xen_sub.obsm["spatial"] = np.asarray(
        xen_sub.obsm["align_spatial"], dtype=float
    ).copy()
    vhd_sub.obsm["spatial"] = np.asarray(
        vhd_sub.obsm["align_spatial"], dtype=float
    ).copy()
    xen_sub.obs["CellType"] = pd.Categorical(
        xen_alignment_labels[xen_mask].astype(str).to_numpy(),
        categories=celltypes,
    )
    vhd_sub.obs["CellType"] = pd.Categorical(
        vhd_alignment_labels[vhd_mask].astype(str).to_numpy(),
        categories=celltypes,
    )
    return xen_sub, vhd_sub, qc


def save_second_pass_coordinates(
    vhd_original: ad.AnnData,
    vhd_second: ad.AnnData,
    args: argparse.Namespace,
    path: Path,
) -> None:
    original = np.asarray(
        vhd_original.obsm["spatial_original_um"], dtype=float
    )[:, :2]
    alignment_input = np.asarray(
        vhd_original.obsm["spatial"], dtype=float
    )[:, :2]
    second_xy = np.asarray(
        vhd_second.obsm["align_spatial"], dtype=float
    )[:, :2]
    if not (len(original) == len(second_xy) == vhd_original.n_obs):
        raise ValueError("VHD observation count changed across alignment passes")

    output = pd.DataFrame(
        {
            "donor": args.donor,
            "xenium_sample": args.xenium_sample,
            "vhd_sample": args.vhd_label,
            "vhd_obs_id": vhd_original.obs_names.astype(str),
            "x_original_um": original[:, 0],
            "y_original_um": original[:, 1],
            "x_alignment_input_um": alignment_input[:, 0],
            "y_alignment_input_flipped_um": alignment_input[:, 1],
            "x_second_pass": second_xy[:, 0],
            "y_second_pass": second_xy[:, 1],
            "xenium_z_height": args.xenium_z_height,
            "vhd_z_height": args.vhd_z_height,
            "Spatial_Domain": (
                vhd_original.obs["Spatial_Domain"]
                .astype("string")
                .fillna("")
                .astype(str)
                .to_numpy()
            ),
        }
    )
    if output["vhd_obs_id"].duplicated().any():
        raise ValueError("VHD observation IDs are not unique")
    if output[["x_second_pass", "y_second_pass"]].isna().any().any():
        raise ValueError("Saved VHD coordinates contain missing values")
    output.to_csv(path, index=False)


def plot_by_celltype_overlay(
    xen: ad.AnnData,
    vhd: ad.AnnData,
    spatial_key: str,
    xenium_name: str,
    vhd_name: str,
    output_path: Path,
    display_name: str,
    xenium_labels: tuple[str, ...],
    vhd_labels: tuple[str, ...],
) -> None:
    """Plot selected labels with VHD beneath Xenium to show overlap."""
    xen_mask = (
        xen.obs["CellType"]
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
    xen_xy = np.asarray(xen.obsm[spatial_key], dtype=float)[xen_mask, :2]
    vhd_xy = np.asarray(vhd.obsm[spatial_key], dtype=float)[vhd_mask, :2]
    print(
        f"{display_name}: Xenium {list(xenium_labels)} "
        f"cells={len(xen_xy):,}; VHD {list(vhd_labels)} "
        f"observations={len(vhd_xy):,}"
    )

    figure, axis = plt.subplots(figsize=(5.0, 4.5))
    # Use equal point sizes and distinguish the layers only by opacity.
    axis.scatter(
        vhd_xy[:, 0],
        vhd_xy[:, 1],
        s=BY_CELLTYPE_VHD_POINT_SIZE,
        alpha=BY_CELLTYPE_VHD_ALPHA,
        color=FIRST_PASS.MOVING_COLOR,
        linewidths=0,
        label=vhd_name,
        zorder=1,
    )
    axis.scatter(
        xen_xy[:, 0],
        xen_xy[:, 1],
        s=BY_CELLTYPE_XENIUM_POINT_SIZE,
        alpha=BY_CELLTYPE_XENIUM_ALPHA,
        color=FIRST_PASS.REFERENCE_COLOR,
        linewidths=0,
        label=xenium_name,
        zorder=2,
    )
    axis.set_aspect("equal")
    axis.set_xlabel("x (um)")
    axis.set_ylabel("y (um)")
    axis.set_title(
        f"{display_name}: Second-pass aligned: "
        f"{xenium_name} → {vhd_name}",
        fontsize=9,
    )
    axis.legend(
        loc="upper right", fontsize=6, markerscale=8, frameon=False
    )
    figure.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    args = parse_args()
    root = git_root()
    resolve_paths(args, root)
    paths = output_paths(args)
    device = (
        ("cuda" if torch.cuda.is_available() else "cpu")
        if args.device == "auto"
        else args.device
    )
    if device == "cuda" and not torch.cuda.is_available():
        raise RuntimeError("--device cuda requested, but CUDA is unavailable")

    print("Device:", device)
    print("Donor:", args.donor)
    print("Fixed Xenium reference:", args.xenium_sample)
    print("Moving VHD slice:", args.vhd_label)
    print("VHD z height (um):", args.vhd_z_height)
    print("Realignment cell types:", args.realign_celltypes)
    print("First-pass VHD coordinates:", args.first_pass_coordinates)

    xen_full, xenium_slices, _ = FIRST_PASS.load_xenium(args)
    vhd_full = FIRST_PASS.load_vhd(args)
    print("Xenium z height (um):", args.xenium_z_height)

    FIRST_PASS.plot_consecutive_xenium_overlays(
        xenium_slices, args.donor, paths["xenium_sanity"]
    )

    xen_first, vhd_first = load_first_pass_coordinates(
        xen_full, vhd_full, args
    )
    FIRST_PASS.plot_pair_overlay(
        xen_first,
        vhd_first,
        "align_spatial",
        args.xenium_sample,
        args.vhd_label,
        paths["first_pass_plot"],
        stage_label="First-pass input",
    )

    xen_sub, vhd_sub, celltype_qc = prepare_label_subsets(
        xen_first, vhd_first, args.realign_celltypes
    )
    celltype_qc.to_csv(paths["celltype_qc"], index=False)
    print(celltype_qc.to_string(index=False))

    second_pass_transformation = st.align.morpho_align_transformation(
        models=[xen_sub, vhd_sub],
        spatial_key="spatial",
        key_added="align_spatial",
        device=device,
        verbose=True,
        max_iter=args.max_iter,
        partial_robust_level=args.partial_robust_level,
        SVI_mode=False,
        rep_layer=["CellType"],
        rep_field=["obs"],
        dissimilarity=["label"],
    )

    with paths["transformation"].open("wb") as handle:
        pickle.dump(
            second_pass_transformation,
            handle,
            protocol=pickle.HIGHEST_PROTOCOL,
        )

    # Apply the second-pass transformation directly to the Module 01
    # first-pass coordinates.
    xen_second, vhd_second = apply_transformation(
        xen_first, vhd_first, second_pass_transformation
    )
    save_second_pass_coordinates(
        vhd_full,
        vhd_second,
        args,
        paths["coordinates"],
    )
    FIRST_PASS.plot_pair_overlay(
        xen_second,
        vhd_second,
        "align_spatial",
        args.xenium_sample,
        args.vhd_label,
        paths["second_pass_plot"],
        stage_label="Second-pass aligned",
    )

    by_celltype_dir = args.plot_dir / "by_celltype"
    plot_celltypes = list(
        dict.fromkeys([*SHARED_LABEL_PLOTS, *args.realign_celltypes])
    )
    for cell_type in plot_celltypes:
        plot_by_celltype_overlay(
            xen_second,
            vhd_second,
            "align_spatial",
            args.xenium_sample,
            args.vhd_label,
            (
                by_celltype_dir
                / (
                    f"{args.pair_tag}_"
                    f"{FIRST_PASS.safe_filename(cell_type)}_"
                    "second_spateo_um.png"
                )
            ),
            display_name=cell_type,
            xenium_labels=(cell_type,),
            vhd_labels=(cell_type,),
        )

    plot_by_celltype_overlay(
        xen_second,
        vhd_second,
        "align_spatial",
        args.xenium_sample,
        args.vhd_label,
        (
            by_celltype_dir
            / (
                f"{args.pair_tag}_{CHAT_PVALB_INH_NAME}_"
                "second_spateo_um.png"
            )
        ),
        display_name=CHAT_PVALB_INH_NAME,
        xenium_labels=CHAT_PVALB_INH_XENIUM_LABELS,
        vhd_labels=CHAT_PVALB_INH_VHD_LABELS,
    )

    metadata = {
        "donor": args.donor,
        "fixed_reference": {
            "technology": "Xenium",
            "sample": args.xenium_sample,
            "xenium_z_height": args.xenium_z_height,
            "coordinate_source": str(args.aligned_coordinates.resolve()),
        },
        "moving_model": {
            "technology": "VisiumHD",
            "sample": args.vhd_label,
            "vhd_z_height": args.vhd_z_height,
            "input": str(args.vhd_h5ad.resolve()),
        },
        "realign_celltypes": args.realign_celltypes,
        "realign_counts": celltype_qc.to_dict(orient="records"),
        "by_celltype_plots": [
            *plot_celltypes,
            CHAT_PVALB_INH_NAME,
        ],
        "chat_pvalb_inh_mapping": {
            "xenium_CellType": list(CHAT_PVALB_INH_XENIUM_LABELS),
            "vhd_Spatial_Domain": list(CHAT_PVALB_INH_VHD_LABELS),
        },
        "first_pass_coordinate_input": str(
            args.first_pass_coordinates.resolve()
        ),
        "second_pass_transformation_output": str(
            paths["transformation"].resolve()
        ),
        "second_pass_coordinate_output": str(
            paths["coordinates"].resolve()
        ),
        "model_order": [args.xenium_sample, args.vhd_label],
        "fixed_model_unchanged": True,
        "refinement_space": "first-pass aligned Xenium coordinate space",
        "second_pass_formula": (
            "vhd_second_pass = "
            "vhd_first_pass @ Rotation.T + Translation"
        ),
        "device": device,
        "max_iter": args.max_iter,
        "partial_robust_level": args.partial_robust_level,
        "plot_directory": str(args.plot_dir.resolve()),
    }
    paths["metadata"].write_text(json.dumps(metadata, indent=2) + "\n")

    print("Saved second-pass transformation:", paths["transformation"])
    print("Saved second-pass VHD coordinates:", paths["coordinates"])
    print("Saved plots:", args.plot_dir)
    print("Done.")


if __name__ == "__main__":
    main()


# Example usage
# -------------
# cd /dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC/code/11_Xenium_VisiumHD_alignment
# conda activate /dcs04/hicks/data/multi-sample-alignment-benchmark/envs/Spateo/
#
# python 05_Spateo_Xenium_VHD_second_pass.py \
#   --donor Br6660 \
#   --xenium-sample Br6660_Nac10_4080 \
#   --vhd-sample VHD_H1_XKYDCP3_D1 \
#   --vhd-z-height 4120 \
#   --realign-celltypes D1_Island_B \
#   --device cuda \
#   --partial-robust-level 100 \
#   --overwrite
#
# Multiple labels can be supplied after one flag:
#   --realign-celltypes D1_Island_A D1_Island_B
