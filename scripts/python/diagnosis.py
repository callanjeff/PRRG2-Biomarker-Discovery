#!/usr/bin/env python3
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT = SCRIPT_DIR.parent.parent
SURV = PROJECT / "data" / "metadata" / "SURVIVAL_PANCANCER-9_Cohorts"
EXPR = PROJECT / "data" / "processed" / "rna_seq-PANCANCER.tsv"

surv = pd.read_csv(SURV, dtype=str)
expr = pd.read_csv(EXPR, sep="\t", index_col=0)

print("Survival patients:", surv["_PATIENT"].nunique())
print("Expression matrix (genes × samples):", expr.shape)

# Keep the full sample IDs for parsing
sample_ids = pd.Series(expr.columns, name="sample_id")


# Parse TCGA pieces
def tcga_patient(s):
    return s[:12]


def tcga_sample_type(s):
    return s[13:15]  # '01', '11', etc.


def tcga_vial(s):
    return s[15:16]  # 'A','B',...


meta = pd.DataFrame(
    {
        "sample_id": sample_ids,
        "patient": sample_ids.map(tcga_patient),
        "type": sample_ids.map(tcga_sample_type),
        "vial": sample_ids.map(tcga_vial),
    }
)

# Focus on Primary Tumor (01). If you want tumor+normal pairs, change this.
meta_tumor = meta[meta["type"] == "01"].copy()

# If multiple tumor aliquots per patient, prefer vial A; otherwise first by lexical
meta_tumor["vial_rank"] = meta_tumor["vial"].rank(method="dense")  # A<B<C...
choice = meta_tumor.sort_values(["patient", "vial"]).drop_duplicates(  # A before B
    "patient", keep="first"
)

# Subset expression to chosen tumor samples
expr_tumor = expr.loc[:, choice["sample_id"].tolist()].copy()

# Rename columns to patient IDs (now 1 col per patient)
expr_tumor.columns = choice["patient"].values

# Align to survival patients (intersection)
surv_patients = set(surv["_PATIENT"].astype(str))
expr_patients = set(expr_tumor.columns.astype(str))
overlap = sorted(expr_patients & surv_patients)

expr_aligned = expr_tumor.loc[:, overlap].copy()

# Diagnostics
print("Unique tumor samples chosen:", expr_tumor.shape[1])
print("Patients in survival only:", len(surv_patients - expr_patients))
print("Patients in RNA-seq only:", len(expr_patients - surv_patients))
print("Overlap patients (post-dedup):", len(overlap))
print("Final aligned matrix (genes × patients):", expr_aligned.shape)

# Optional: save for downstream survival/GSVA
out = PROJECT / "results" / "tables" / "expr_tumor_aligned_to_survival.tsv"
out.parent.mkdir(parents=True, exist_ok=True)
expr_aligned.to_csv(out, sep="\t")
print(f"Saved: {out}")
