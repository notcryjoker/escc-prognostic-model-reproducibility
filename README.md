# Apparent versus generalizable performance of an ESCC transcriptomic prognostic model

Analysis code and intermediate results accompanying the manuscript *"Apparent versus generalizable performance of an ESCC transcriptomic prognostic model: quantifying selection leakage and overfitting, with a TRIPOD-informed audit and PROBAST assessment of 205 studies"* (submitted to *Journal of Translational Medicine*).

**Repository:** <https://github.com/notcryjoker/escc-prognostic-model-reproducibility>

This repository contains everything needed to reproduce every analysis reported in the manuscript. Raw expression data are **not** included (they are public and can be re-downloaded in minutes); all derived intermediate results (65 MB) **are** included so that results can be verified without re-running the full pipeline.

## Repository structure

```
scripts/                 Analysis scripts, numbered in execution order
  01–09                  TCGA + GEO preprocessing, DEG analysis, prognostic screening, LASSO
  10–13                  Main figures, flowchart, figure organisation, supplementary analysis
  14–17                  Split validation, fixed-model split, nested cross-validation
  18–27                  Literature audit (205 studies), TRIPOD/PROBAST assessment, figures
  28–35                  PDF builders, pipeline bootstrap, control split, gate-artifact regeneration
results/                 All intermediate result files (TSV/RDS/PNG/logs) produced by scripts 01–33
data/msigdb/             MSigDB gene-set subsets used for enrichment (GMT)
manuscript/supplementary/  16 supplementary files (BMC "Additional file N" convention)
```

**Scope note.** This repository is the reproducibility package only: analysis code, derived intermediate results, gene-set inputs and the supplementary data files. The manuscript text, cover letter, submission checklist and internal audit/decision records are deliberately **not** published here — the article is unpublished, and those materials contain pre-decisional response strategy. The complete package is also archived at Zenodo (<https://doi.org/10.5281/zenodo.23643162>); the raw TCGA/GEO expression data are public and excluded (see below).

## Data availability

All analysed data are publicly available:

- **TCGA-ESCA** — GDC Data Portal, <https://portal.gdc.cancer.gov>
- **GSE20347, GSE45670, GSE53624** — Gene Expression Omnibus, <https://www.ncbi.nlm.nih.gov/geo>
- **MSigDB** gene sets (KEGG, GO BP, Hallmark subsets) — `data/msigdb/` (included)

Raw data download and preprocessing follow `scripts/01_TCGA_preprocess.R`, `scripts/02_GEO_preprocess.R`, `scripts/01b_TCGA_to_symbol.R`, `scripts/02b_GSE53624_preprocess.R` and `scripts/download_gse53624.sh`. The GSE53624 external-validation attempt failed for documented technical reasons (platform annotation without gene symbols; see `scripts/extract_platform_table.py` output and the manuscript Limitations).

## Reproducing the pipeline

Environment: R (>= 4.3) with `survival`, `glmnet`, `limma`, `clusterProfiler`, `GSVA`, `pheatmap`, `survminer`, `rms`, `ggplot2`; Python (>= 3.11) with `pandas`. Gene-symbol annotation via BioMart.

Run scripts in numeric order (each script reads only outputs of lower-numbered scripts):

```bash
# 1. Download raw data (see Data availability above) into data/tcga/ and data/geo/
# 2. Run the R pipeline in numeric order (01 … 17), then the Python audit scripts (18, 20, 21, 23, 24, 26, 27)
for s in scripts/01_*.R scripts/01b_*.R scripts/02_*.R scripts/02b_*.R scripts/0[3-9]_*.R scripts/1*.R; do Rscript "$s"; done
for p in 18_literature_audit 20_probst_fulltext 21_tripod_audit 23_get_fulltext 24_probst_assessment 26_rule_agreement 27_make_table_s8; do python3 scripts/$p.py; done
# 3. Figures and review PDFs (28–31); statistics (32–33); gate artifacts (35)
```

Random seeds are fixed inside each script (nested CV seed 20231020; pipeline bootstrap B = 500, 489/500 completed). Full expected outputs, including the bootstrap and control-split summaries, are in `results/`.

## Literature-audit dataset

The complete audit of 205 published ESCC prognostic-model studies (PubMed, 2015–2026) — extracted reporting classifications, TRIPOD item coding and PROBAST assessments — is included as `results/literature_audit/` and as supplementary **Additional file 1** with the manuscript. A second, independent rule set was used for agreement estimation (Additional file 3; mean Cohen's κ = 0.802).

## License

Code: MIT License. Data derived from TCGA/GEO remain subject to their original terms. Audit dataset: CC0.

## Contact

Yang Xu — Department of Thoracic Surgery, The First People's Hospital of Yinchuan (Second Clinical Medical School, Ningxia Medical University), Yinchuan, Ningxia, China. Corresponding author: 2581917560@qq.com.
