#!/usr/bin/env python3
"""
Survival QC + Kaplan–Meier for TCGA-style table.

Supported input: CSV or TSV (auto-detected).
Expected columns include (case/spacing/punct. will be normalized):
row ID, sample, _PATIENT, cancer type abbreviation, age_at_initial_pathologic_diagnosis,
gender, race, ajcc_pathologic_tumor_stage, clinical_stage, histological_type,
histological_grade, initial_pathologic_dx_year, menopause_status, birth_days_to,
vital_status, tumor_status, last_contact_days_to, death_days_to, cause_of_death,
new_tumor_event_type, new_tumor_event_site, new_tumor_event_site_other,
new_tumor_event_dx_days_to, treatment_outcome_first_course, margin_status,
residual_tumor, OS, OS.time, DSS, DSS.time, DFI, DFI.time, PFI, PFI.time, Redaction

Outputs:
results/tables/survival_summary_overall.tsv
results/tables/survival_summary_by_cancer.tsv
results/tables/survival_long.tsv
results/figures/km_<EP>_overall.png
results/figures/km_<EP>_by_cancer.png

Usage:
python scripts/python/survival_qc_and_km.py \
    --input data/metadata/SURVIVAL_PANCANCER-9_Cohorts.csv \
    --outdir results \
    --endpoints OS DSS PFI \
    --min-n 50
"""

from __future__ import annotations

import argparse
import os
import re
from typing import List, Tuple

import pandas as pd

# Optional plotting (lifelines/matplotlib). Script still runs without them.
try:
    import matplotlib.pyplot as plt
    from lifelines import KaplanMeierFitter

    HAVE_LIFELINES = True
except Exception:
    HAVE_LIFELINES = False


def pick_col(df: pd.DataFrame, candidates: list[str], required_name: str) -> str:
    """
    Return the first column from `candidates` that exists in df.columns.
    Raise a helpful error if none found.
    """
    for c in candidates:
        if c in df.columns:
            return c
    raise KeyError(
        f"None of the candidate columns for '{required_name}' "
        f"found. Tried: {candidates}. Available: {list(df.columns)[:50]} ..."
    )


def norm_cols(df: pd.DataFrame) -> pd.DataFrame:
    """Normalize column names: spaces/periods/odd chars -> underscores, lowercase."""
    df = df.copy()
    df.columns = [
        re.sub(r"_+", "_", re.sub(r"[^0-9a-zA-Z]+", "_", c)).strip("_").lower()
        for c in df.columns
    ]
    return df


def coerce_numeric(df: pd.DataFrame, cols: List[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = pd.to_numeric(out[c], errors="coerce")
    return out


def read_table_auto(path: str) -> pd.DataFrame:
    """
    Auto-detect delimiter (comma or tab) and read as strings.
    Uses pandas engine='python' to better handle quoted fields.
    """
    # sep=None triggers the python engine's built-in sniffer
    return pd.read_csv(path, sep=None, engine="python", dtype=str)


def build_long(df: pd.DataFrame, endpoints: List[str]) -> pd.DataFrame:
    """
    Create long-format survival table with columns:
    ['patient', 'sample', 'cancer', 'endpoint', 'time', 'event']
    """

    # Handle common naming variants after normalization
    patient_col_src = pick_col(
        df,
        candidates=["patient", "_patient", "bcr_patient_barcode", "case_submitter_id"],
        required_name="patient",
    )
    sample_col_src = pick_col(
        df,
        candidates=["sample", "sample_id", "aliquot_barcode"],
        required_name="sample",
    )
    cancer_col_src = pick_col(
        df,
        candidates=["cancer_type_abbreviation", "cancer", "project_id", "disease"],
        required_name="cancer type abbreviation",
    )

    keep = []
    for ep in endpoints:
        ep = ep.upper()
        time_col = f"{ep.lower()}_time"
        ep_col = ep.lower()
        if ep_col not in df.columns or time_col not in df.columns:
            # silently skip missing endpoints
            continue

        sub = (
            df[[patient_col_src, sample_col_src, cancer_col_src, ep_col, time_col]]
            .rename(
                columns={
                    patient_col_src: "patient",
                    sample_col_src: "sample",
                    cancer_col_src: "cancer",
                    ep_col: "event",
                    time_col: "time",
                }
            )
            .copy()
        )

        sub["time"] = pd.to_numeric(sub["time"], errors="coerce")
        sub["event"] = pd.to_numeric(sub["event"], errors="coerce")
        sub.loc[sub["time"].isna() | (sub["time"] <= 0), "time"] = pd.NA
        sub["endpoint"] = ep
        keep.append(sub)

    if not keep:
        raise ValueError("None of the requested endpoints were found in the table.")
    out = pd.concat(keep, axis=0, ignore_index=True)

    out = out.dropna(subset=["time"])
    out["event"] = out["event"].fillna(0).clip(lower=0, upper=1).astype(int)
    out["cancer"] = out["cancer"].astype(str).str.upper()
    return out


def summarize_survival(surv_long: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Basic counts + censoring summaries overall and by cancer."""

    def agg(df):
        return pd.Series(
            {
                "n": int(df.shape[0]),
                "events": int(df["event"].sum()),
                "censored": int((1 - df["event"]).sum()),
                "event_rate": float(df["event"].mean()),
                "median_time": float(df["time"].median(skipna=True)),
            }
        )

    overall = surv_long.groupby(["endpoint"], dropna=False).apply(agg).reset_index()
    by_cancer = (
        surv_long.groupby(["endpoint", "cancer"], dropna=False).apply(agg).reset_index()
    )
    return overall, by_cancer


def km_plot_overall(surv_long: pd.DataFrame, endpoint: str, outpath: str):
    if not HAVE_LIFELINES:
        return
    df = surv_long.query("endpoint == @endpoint").copy()
    if df.empty:
        return
    kmf = KaplanMeierFitter()
    kmf.fit(
        durations=df["time"], event_observed=df["event"], label=f"{endpoint} (overall)"
    )
    plt.figure()
    kmf.plot(ci_show=True)
    plt.title(f"Kaplan–Meier: {endpoint} (Overall)")
    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def km_plot_by_cancer(
    surv_long: pd.DataFrame,
    endpoint: str,
    outpath: str,
    min_n: int = 50,
    top_k: int = 9,
):
    if not HAVE_LIFELINES:
        return
    df = surv_long.query("endpoint == @endpoint").copy()
    counts = df.groupby("cancer")["patient"].count().sort_values(ascending=False)
    cancers = [c for c, n in counts.items() if n >= min_n][:top_k]
    if not cancers:
        return

    plt.figure()
    for c in cancers:
        sub = df[df["cancer"] == c]
        if sub.empty:
            continue
        kmf = KaplanMeierFitter()
        kmf.fit(sub["time"], event_observed=sub["event"], label=c)
        kmf.plot(ci_show=False)

    plt.title(f"Kaplan–Meier: {endpoint} by cancer (N≥{min_n})")
    plt.xlabel("Time (days)")
    plt.ylabel("Survival probability")
    plt.legend(title="Cancer", fontsize=8)
    plt.tight_layout()
    plt.savefig(outpath, dpi=180)
    plt.close()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--input", required=True, help="Path to SURVIVAL_PANCANCER-9_Cohorts.(csv|tsv)"
    )
    ap.add_argument(
        "--outdir", default="results", help="Base output directory (default: results)"
    )
    ap.add_argument(
        "--endpoints",
        nargs="+",
        default=["OS"],
        help="Endpoints (e.g., OS DSS DFI PFI)",
    )
    ap.add_argument(
        "--min-n", type=int, default=50, help="Min N per cancer for KM by-cancer plot"
    )
    args = ap.parse_args()

    # I/O setup
    tables_dir = os.path.join(args.outdir, "tables")
    figs_dir = os.path.join(args.outdir, "figures")
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(figs_dir, exist_ok=True)

    # Load & normalize
    df = read_table_auto(args.input)
    df = norm_cols(df)

    # Coerce numeric for known fields if present
    numeric_cols = [
        "age_at_initial_pathologic_diagnosis",
        "birth_days_to",
        "last_contact_days_to",
        "death_days_to",
        "os",
        "os_time",
        "dss",
        "dss_time",
        "dfi",
        "dfi_time",
        "pfi",
        "pfi_time",
    ]
    df = coerce_numeric(df, [c for c in numeric_cols if c in df.columns])

    # Build tidy survival table
    endpoints = [e.upper() for e in args.endpoints]
    surv_long = build_long(df, endpoints=endpoints)

    # Summaries
    overall, by_cancer = summarize_survival(surv_long)

    # Write tables
    overall_path = os.path.join(tables_dir, "survival_summary_overall.tsv")
    by_cancer_path = os.path.join(tables_dir, "survival_summary_by_cancer.tsv")
    long_path = os.path.join(tables_dir, "survival_long.tsv")
    overall.to_csv(overall_path, sep="\t", index=False)
    by_cancer.to_csv(by_cancer_path, sep="\t", index=False)
    surv_long.to_csv(long_path, sep="\t", index=False)

    # Plots
    for ep in endpoints:
        km_plot_overall(surv_long, ep, os.path.join(figs_dir, f"km_{ep}_overall.png"))
        km_plot_by_cancer(
            surv_long,
            ep,
            os.path.join(figs_dir, f"km_{ep}_by_cancer.png"),
            min_n=args.min_n,
            top_k=9,
        )

    print(f"[OK] Wrote:\n  {overall_path}\n  {by_cancer_path}\n  {long_path}")
    if HAVE_LIFELINES:
        print(f"[OK] Plots in {figs_dir}")
    else:
        print("[WARN] lifelines/matplotlib not available; skipped plots.")


if __name__ == "__main__":
    main()
