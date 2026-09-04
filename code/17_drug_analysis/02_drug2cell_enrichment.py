"""Identify drugs whose target-gene scores are enriched in each spatial domain
and cell type (Wilcoxon rank-sum, one group vs rest, on the drug2cell score
matrix), and flag a curated panel of drugs of abuse / psychiatric drugs.

Outputs (processed-data/17_drug_analysis/drug2cell/):
  - drug_enrichment_{Spatial_Domain,labels}.tsv.gz : full Wilcoxon results
  - drug_enrichment_{...}_highlight.tsv            : curated-drug subset
"""

from pathlib import Path

import pandas as pd
import scanpy as sc
import session_info

BASE = Path("/dcs05/lieber/marmaypag/xenium_NAC_LIBD4125/xenium_NAC")
D2C_DIR = BASE / "processed-data/17_drug_analysis/drug2cell"
GROUP_COLS = ["Spatial_Domain", "labels"]

# Curated highlight panel: substring match against ChEMBL preferred names
# (drug2cell var names look like "CHEMBL521|IBUPROFEN").
HIGHLIGHT = {
    "opioid": [
        "MORPHINE", "FENTANYL", "OXYCODONE", "HYDROCODONE", "HEROIN",
        "METHADONE", "BUPRENORPHINE", "NALTREXONE", "NALOXONE", "TRAMADOL",
    ],
    "stimulant": [
        "COCAINE", "AMPHETAMINE", "METHAMPHETAMINE", "METHYLPHENIDATE",
        "MDMA", "MODAFINIL",
    ],
    "nicotine": ["NICOTINE", "VARENICLINE", "CYTISINE"],
    "alcohol_related": ["DISULFIRAM", "ACAMPROSATE", "NALMEFENE"],
    "cannabinoid": ["DRONABINOL", "NABILONE", "CANNABIDIOL", "RIMONABANT"],
    "antipsychotic": [
        "HALOPERIDOL", "RISPERIDONE", "CLOZAPINE", "OLANZAPINE",
        "QUETIAPINE", "ARIPIPRAZOLE", "CHLORPROMAZINE","CARIPRAZINE",
	"BREXPIPRAZOLE",
    ],
    "antidepressant": [
        "FLUOXETINE", "SERTRALINE", "CITALOPRAM", "ESCITALOPRAM",
        "PAROXETINE", "VENLAFAXINE", "DULOXETINE", "BUPROPION",
        "AMITRIPTYLINE", "MIRTAZAPINE",
    ],
    "mood_stabilizer": ["LITHIUM", "VALPROIC", "LAMOTRIGINE", "CARBAMAZEPINE"],
    "sedative": ["DIAZEPAM", "ALPRAZOLAM", "LORAZEPAM", "ZOLPIDEM", "MIDAZOLAM"],
    "dissociative_psychedelic": [
        "KETAMINE", "ESKETAMINE",
        "PSILOCYBIN", "PSILOCIN", "LYSERGIDE", "MESCALINE",
        "DIMETHYLTRYPTAMINE",
    ],
}


def classify(drug_name: str) -> str:
    upper = drug_name.upper()
    for cls, names in HIGHLIGHT.items():
        if any(n in upper for n in names):
            return cls
    return ""


drug_adata = sc.read_h5ad(D2C_DIR / "VHD_drug2cell_scores.h5ad")
print(drug_adata)

# Remove cells outside the NAc (Hypo and Excitatory spatial domains)
# and Excitatory cell-type labels from any remaining spatial domains
sd_exclude = ["Hypo", "Excitatory"]
label_exclude = ["Excitatory"]
keep = ~drug_adata.obs["Spatial_Domain"].isin(sd_exclude) & ~drug_adata.obs["labels"].isin(label_exclude)
print(f"Excluding {(~keep).sum()} cells (Spatial_Domain in {sd_exclude} or labels in {label_exclude})")
drug_adata = drug_adata[keep].copy()

for group in GROUP_COLS:
    print(f"--- Enrichment by {group}")
    drug_adata.obs[group] = drug_adata.obs[group].astype("category")
    sc.tl.rank_genes_groups(drug_adata, groupby=group, method="wilcoxon")
    res = sc.get.rank_genes_groups_df(drug_adata, group=None)
    res = res.rename(columns={"names": "drug", "group": group})
    res["highlight_class"] = [classify(d) for d in res["drug"]]

    res.to_csv(
        D2C_DIR / f"drug_enrichment_{group}.tsv.gz",
        sep="\t", index=False, compression="gzip",
    )
    res[res["highlight_class"] != ""].to_csv(
        D2C_DIR / f"drug_enrichment_{group}_highlight.tsv",
        sep="\t", index=False,
    )
    print(f"Wrote {group}: {res.shape[0]} rows, "
          f"{(res['highlight_class'] != '').sum()} highlighted")

print("Done.")

session_info.show()
