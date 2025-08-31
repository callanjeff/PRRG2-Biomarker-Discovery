# PRRG2-Biomarker_Discovery

Independent graduate-level research project exploring **PRRG2**
(Proline Rich and Gla Domain 2) as a potential cancer biomarker
using multi-omics datasets.


## Goals

- Compare PRRG2 expression in tumor vs normal (TCGA, GTEx).
- Assess immune infiltration (TIMER, CIBERSORTx).
- Evaluate survival, stage, and clinical correlations.
- Extend to proteomics (CPTAC) or validation (GEO) when possible.


## Workflow

1. **KNIME** → preprocessing
2. **Python** → wrangling, QC, joins
3. **R/Bioconductor** → DE, survival, GSVA, immune plots


## Outputs

- `results/figures/`
- `results/tables/`
