"""Score every ChEMBL drug's target gene set on Visium-HD cells with drug2cell.

drug2cell (Kanemaru et al., Nature 2023) computes, per cell, the mean
log-normalized expression of each drug's target genes using the packaged
ChEMBL-derived drug-target dictionary (pChEMBL-active human targets).

The bundled dictionary (drug2cell v0.1.2) only includes Phase-4 (approved)
drugs, so psychedelics like psilocybin (Phase 3), LSD (Phase 2), and
mescaline (Phase 1) are absent.  We manually inject them below using
target-gene annotations from ChEMBL (queried 2026-09-01).

Input : VHD_sfe_logcounts.h5ad (X = logcounts; obs has Spatial_Domain, labels)
Output: processed-data/17_drug_analysis/drug2cell/VHD_drug2cell_scores.h5ad
        (obs = cells with groupings, var = drugs)
"""

from collections import ChainMap
from pathlib import Path

import drug2cell as d2c
import scanpy as sc
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
H5AD = BASE / "processed-data/HD_Full_Analysis/h5ad/VHD_sfe_logcounts.h5ad"
OUT = BASE / "processed-data/17_drug_analysis/drug2cell/VHD_drug2cell_scores.h5ad"

print(f"Loading {H5AD}")
adata = sc.read_h5ad(H5AD)
print(adata)

# Remove cells outside the NAc (Hypo and Excitatory spatial domains)
# and Excitatory cell-type labels from any remaining spatial domains
sd_exclude = ["Hypo", "Excitatory"]
label_exclude = ["Excitatory"]
keep = ~adata.obs["Spatial_Domain"].isin(sd_exclude) & ~adata.obs["labels"].isin(label_exclude)
print(f"Excluding {(~keep).sum()} cells (Spatial_Domain in {sd_exclude} or labels in {label_exclude})")
adata = adata[keep].copy()

# ---------------------------------------------------------------------------
# Build augmented drug-target dictionary
# ---------------------------------------------------------------------------
# Load the default nested dict (ATC categories -> drugs -> gene lists)
default_targets = d2c.data.chembl()
# Flatten the nested dict into {drug: [genes]}
flat_targets = dict(ChainMap(*[default_targets[cat] for cat in default_targets]))

# Custom psychedelic entries — ChEMBL IDs and human target gene symbols from
# ChEMBL binding assay data (pChEMBL ≥ 5, human targets).
# Psilocin (CHEMBL65547) is the active metabolite of psilocybin; it has far
# more target data than psilocybin itself, so we include both.
CUSTOM_PSYCHEDELICS = {
    "CHEMBL194378|PSILOCYBIN": [
        "HTR2A", "HTR2B", "HTR1A",
    ],
    "CHEMBL65547|PSILOCIN": [
        "HTR2A", "HTR2B", "HTR2C", "HTR1A", "HTR1B", "HTR1D", "HTR1E",
        "HTR5A", "HTR6", "HTR7", "DRD1", "DRD3", "ADORA2A", "ADORA2B",
    ],
    "CHEMBL263881|LYSERGIDE": [
        "HTR1A", "HTR2A", "HTR2B", "HTR2C", "HTR1D", "HTR1E",
        "HTR5A", "HTR6", "HTR7",
        "DRD2", "DRD3", "DRD4", "DRD5",
        "ADRB1", "ADRB2", "HRH1",
    ],
    "CHEMBL26687|MESCALINE": [
        "HTR2A", "HTR2C", "HTR1A",
    ],
    "CHEMBL12420|N,N-DIMETHYLTRYPTAMINE": [
        "HTR2A", "HTR2B", "HTR2C", "HTR1A",
    ],
}
flat_targets.update(CUSTOM_PSYCHEDELICS)
print(f"Added {len(CUSTOM_PSYCHEDELICS)} custom psychedelic entries "
      f"({len(flat_targets)} total drugs)")

# X already holds logcounts in this export
d2c.score(adata, targets=flat_targets, use_raw=False)

drug_adata = adata.uns["drug2cell"]
# Carry groupings over for the enrichment step
for col in ["Sample", "Spatial_Domain", "labels", "X_coords", "Y_coords"]:
    if col in adata.obs.columns:
        drug_adata.obs[col] = adata.obs[col].values

print(f"Scored {drug_adata.shape[1]} drugs x {drug_adata.shape[0]} cells")
drug_adata.write_h5ad(OUT)
print(f"Wrote {OUT}")

session_info.show()
