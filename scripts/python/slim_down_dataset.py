#!/usr/bin/env python3
import sys
from pathlib import Path

import pandas as pd

PROJECT = Path("/Volumes/OWC_EnvoyPRO/PRRG2-Biomarker-Discovery")
INP = PROJECT / "results" / "tables" / "survival_with_expr_FULL.tsv"
OUT = PROJECT / "results" / "tables" / "survival_with_PRRG2.tsv"

# Read (be liberal about types; this also removes the DtypeWarning)
df = pd.read_csv(INP, sep="\t", dtype=str, low_memory=False)


# --- helpers to find columns no matter their exact names ---
def find_first(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None


# Common survival column name variants
time_candidates = [
    "TIME",
    "time",
    "OS_TIME",
    "OS_MONTHS",
    "OS.time",
    "OS_days",
    "days_to_event",
    "days_to_last_followup",
    "days_to_death",
]
status_candidates = [
    "STATUS",
    "status",
    "OS_STATUS",
    "OS",
    "OS.event",
    "event",
    "vital_status",
]
cohort_candidates = [
    "cancer type abbreviation",
    "COHORT",
    "cohort",
    "project",
    "disease_type",
]

pat_col = (
    "_PATIENT"
    if "_PATIENT" in df.columns
    else find_first(df, ["PATIENT", "patient", "barcode"])
)
if pat_col is None:
    sys.exit("❌ Could not find patient ID column (expected '_PATIENT').")

time_col = find_first(df, time_candidates)
status_col = find_first(df, status_candidates)
cohort_col = find_first(df, cohort_candidates)

# If no explicit TIME/STATUS, derive them from typical TCGA fields
if time_col is None:
    # derive time (in days) from days_to_death or days_to_last_followup
    dtd = "days_to_death" if "days_to_death" in df.columns else None
    dtlf = "days_to_last_followup" if "days_to_last_followup" in df.columns else None
    if dtd or dtlf:
        # choose available; prefer death if present, else last follow-up
        df["_TIME_DAYS"] = (
            df[dtd].fillna(df[dtlf]) if (dtd and dtlf) else df[dtd or dtlf]
        )
        time_col = "_TIME_DAYS"
    else:
        print("⚠️ No explicit or derivable TIME column found; proceeding without it.")

if status_col is None and "vital_status" in df.columns:
    # Map vital_status (Alive/Dead) to 0/1
    df["_STATUS"] = (
        df["vital_status"].str.strip().str.upper().map({"DEAD": "1", "ALIVE": "0"})
    )
    status_col = "_STATUS"
elif status_col is None:
    print("⚠️ No explicit or derivable STATUS column found; proceeding without it.")

# Ensure PRRG2 exists
if "PRRG2" not in df.columns:
    # sometimes gene symbols carry spaces or weird casing—report available matches
    matches = [c for c in df.columns if c.upper() == "PRRG2"]
    if matches:
        df.rename(columns={matches[0]: "PRRG2"}, inplace=True)
    else:
        # show a few tail columns to help debug
        tail = df.columns[-10:].tolist()
        sys.exit(f"❌ Column 'PRRG2' not found. Last columns: {tail}")

# Build the output column list
keep = [pat_col]
if cohort_col:
    keep.append(cohort_col)
if time_col:
    keep.append(time_col)
if status_col:
    keep.append(status_col)
keep.append("PRRG2")

out_df = df[keep].copy()

# Normalize header names for convenience downstream
rename_map = {pat_col: "_PATIENT", "PRRG2": "PRRG2"}
if cohort_col:
    rename_map[cohort_col] = "COHORT"
if time_col:
    rename_map[time_col] = "TIME"
if status_col:
    rename_map[status_col] = "STATUS"
out_df.rename(columns=rename_map, inplace=True)

# Save
OUT.parent.mkdir(parents=True, exist_ok=True)
out_df.to_csv(OUT, sep="\t", index=False)
print(f"✅ Saved slim file: {OUT}")
print("   Columns:", list(out_df.columns))
print("   Shape:", out_df.shape)
