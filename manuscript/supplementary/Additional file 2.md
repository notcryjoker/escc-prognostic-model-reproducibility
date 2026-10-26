# Additional file 2 — Literature-audit dataset: notes and column dictionary

**Source dataset for the TRIPOD-informed audit of 205 published ESCC prognostic-model studies (PubMed, 2015–2026).**

## Column dictionary

All T-item columns are binary indicators of what the study's **abstract** explicitly reported (1 = yes, 0 = no). Classification used the predefined keyword rule set in `scripts/21_tripod_audit.py` (primary rule set); agreement with an independent secondary rule set is reported in Additional file 3 (mean Cohen's κ = 0.802).

| Column | Meaning |
|---|---|
| PMID | PubMed identifier |
| Year | Publication year |
| Journal | Journal of publication |
| Title | Article title (as indexed in PubMed) |
| T1a_External_validation_reported | Abstract states independent external validation (TRIPOD item 16, external type) |
| T1b_Internal_split_or_CV_reported | Abstract states internal split-sample or cross-validation (TRIPOD item 16, internal type) |
| T1c_No_validation_mentioned | Abstract mentions neither form of validation |
| T2_Pipeline_rerun_stated | Abstract states the modelling pipeline, including gene/variable selection, was re-run during validation (beyond the TRIPOD minimum) |
| T3_Discrimination_AUC_Cindex_reported | Abstract reports a discrimination measure (AUC or C-index; TRIPOD item 13b) |
| T4_Calibration_reported | Abstract reports calibration (TRIPOD item 13c) |
| T5_Full_model_reported | Abstract reports the full model (coefficients, formula or nomogram; TRIPOD item 15) |
| T6_Sample_size_reported | Abstract reports the sample size (TRIPOD item 17) |
| T7_Multicentre_or_multisource | Abstract reports a multicentre or multi-source design (TRIPOD item 4) |
| Max_reported_AUC_or_Cindex | Highest AUC/C-index value stated in the abstract, where extractable (blank = not extractable) |
| Max_sample_size | Largest sample size stated in the abstract, where extractable |

## Key counts (recomputed from this table)

| Item | Studies | Proportion |
|---|---|---|
| No validation mentioned | 142/205 | 69.3% |
| External validation | 55/205 | 26.8% |
| Internal split/CV | 22/205 | 10.7% |
| Pipeline re-run stated | 0/205 | 0.0% |
| Discrimination reported | 81/205 | 39.5% |
| Calibration reported | 15/205 | 7.3% |
| Full model reported | 85/205 | 41.5% |
| Sample size reported | 35/205 | 17.1% |

## Caveats

- The audit assessed **reporting adequacy in abstracts**, not study conduct: open-access full texts were retrievable for only 3/205 records, so items that cannot be judged from an abstract were not scored.
- Records are sorted by year then PMID.
