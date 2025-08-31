#!/usr/bin/env python3
from pathlib import Path
import pandas as pd
import sys

PROJECT = Path("/Volumes/OWC_EnvoyPRO/PRRG2-Biomarker-Discovery")

# 1) Load the deduped, tumor-only, patient-aligned expression matrix (genes x patients=2991)
expr = pd.read_csv(PROJECT/"results/tables/expr_tumor_aligned_to_survival.tsv",
                   sep="\t", index_col=0)

# 2) Load survival (no extension needed)
surv_path = PROJECT/"data/metadata/SURVIVAL_PANCANCER-9_Cohorts"
surv = pd.read_csv(surv_path, sep=None, engine="python", dtype=str)

# 3) Safety checks
print("Expr shape (genes x patients):", expr.shape)  # expect (20531, 2991)
print("Survival raw rows:", surv.shape[0])

# Make sure survival is unique per patient
dup_count = surv.duplicated(subset=["_PATIENT"]).sum()
if dup_count:
    print(f"⚠️ Survival has {dup_count} duplicate _PATIENT rows; keeping first per patient.")
    surv = surv.drop_duplicates(subset=["_PATIENT"], keep="first")

# 4) Transpose expression -> patients as rows
expr_T = expr.T
expr_T.index.name = "_PATIENT"

# 5) Inner merge by _PATIENT (guarantees 1 row per patient)
merged = surv.merge(expr_T, on="_PATIENT", how="inner")

print("Merged shape (patients x [meta+genes]):", merged.shape)  # expect ~ (2991, 20531 + meta)
print("Unique patients in merged:", merged["_PATIENT"].nunique())

# 6) OPTIONAL: keep only PRRG2 if you want a slim table
# prrg2_only = merged[["_PATIENT", "cancer type abbreviation", "TIME", "STATUS", "PRRG2"]]
# prrg2_only.to_csv(PROJECT/"results/tables/survival_with_PRRG2.tsv", sep="\t", index=False)

# 7) Save full merged matrix
out = PROJECT/"results/tables/survival_with_expr_FULL.tsv"
merged.to_csv(out, sep="\t", index=False)
print(f"✅ Saved: {out}")

