import pandas as pd

# load survival table
surv = pd.read_csv("data/metadata/SURVIVAL_PANCANCER-9_Cohorts.csv", dtype=str)

# load expression matrix
expr = pd.read_csv("data/metadata/rna_seq-PANCANCER.tsv", sep="\t", index_col=0)

print("Survival patients:", surv["_PATIENT"].nunique())
print("Expression samples:", expr.shape)

# normalize sample IDs
expr.columns = expr.columns.str.slice(0,12)   # shorten TCGA-XX-XXXX-XX... → TCGA-XX-XXXX

# check overlap
overlap = set(expr.columns) & set(surv["_PATIENT"])
print("Overlap patients:", len(overlap))

# filter expression to overlap only
expr_sub = expr.loc[:, expr.columns.isin(overlap)]

print("Subset expression shape:", expr_sub.shape)

