# Supplementary Material S9: Experimental validation plan

## Functional-experiment and protein-level validation design based on the findings of this study

This plan translates the bioinformatic findings of this study into an experimental programme that can be tested directly in tissue samples, ordered by priority.

---

## 1. Validation objectives and testable hypotheses

| No. | Finding from this study | Hypothesis to be tested | Validation method |
|---|---|---|---|
| H1 | Infiltration of 19 immune cell types was elevated in the high-risk group, with pDC showing the strongest correlation (ρ = 0.471) | Stromal infiltration by pDC, monocytes/MDSC and immature dendritic cells is significantly higher in high-risk than in low-risk ESCC tumours | Multiplex immunofluorescence (mIF) / multicolour IHC |
| H2 | TIGIT, HAVCR2 (TIM-3), IDO1, LGALS9 and PDCD1 were positively correlated with the risk score | High-risk tumours display a multi-checkpoint co-inhibition phenotype (protein-level corroboration of the mRNA-level finding) | IHC (TIGIT, TIM-3, IDO1, Gal-9, PD-1, PD-L1) + RNA (qRT-PCR / RNAscope) |
| H3 | P3H4, RHOU, TNFSF10, BEX4 and CH25H were selected most frequently across repeated modelling (P3H4 in 100%) | Protein expression of these genes is associated with ESCC prognosis, and their localisation and the functions implied by enrichment analysis (ECM/stroma, immune, metabolic) are concordant | IHC + survival analysis |
| H4 | Enrichment analysis pointed to cell-cycle dysregulation (E2F/G2M/MYC) | High-risk tumours show a higher Ki-67 index and higher expression of cell-cycle proteins (Cyclin B1, p-Rb) | IHC (Ki-67, Cyclin B1, pHH3) |

---

## 2. Sample design (tissue microarray, TMA)

- **Discovery–validation two-arm design**:
  - Cohort A (internal validation): ESCC surgical specimens matched to the population of the discovery cohort, **n ≥ 150 cases** (requirements: follow-up ≥ 5 years, ≥ 60 events, so that the events-per-variable ratio and statistical power are acceptable; both the 17-gene risk score of this study and the IHC markers will be tested for association with survival in this cohort).
  - Cohort B (independent external validation): a second centre and a second population, **n ≥ 100 cases**, in which only the pre-specified markers will be tested (to avoid secondary data mining).
- **TMA construction**: three tumour cores (1.0–1.5 mm diameter) plus one paired adjacent-normal core per case; tumour content (> 60% tumour cells) verified by HE review of the first-section slides.
- **Ethics**: use of archival FFPE samples requires ethics approval and (where applicable) an informed-consent waiver; neoadjuvant treatment history will be recorded as a covariate.

## 3. mIF panel design (testing H1, H2)

**Panel 1 (immune-cell composition) — 6 colours + DAPI:**

| Marker | Target cells | Rationale |
|---|---|---|
| CD303 (BDCA-2) + CD123 | Plasmacytoid dendritic cells (pDC) | The immune cell type with the most significant difference in the high-risk group |
| CD68 + CD14 | Monocytes / macrophages / M-MDSC | Monocytes were elevated in the high-risk group |
| CD11c + CD1a | Conventional DC / immature DC | Immature DC were elevated in the high-risk group |
| CD8 + GZMB | Cytotoxic T cells | Effector-memory CD8⁺ T cells were elevated, but their functional state is uncertain |
| FOXP3 | Regulatory T cells | Mechanistic validation of pDC-induced Treg differentiation |
| Pan-CK | Tumour compartment delineation | Distinguishing tumour epithelium from stroma |

**Panel 2 (checkpoint co-inhibition status) — 4–6 colours:** PD-1, TIM-3, TIGIT, IDO1, PD-L1 (clone 22C3 or SP142), Gal-9.

**Quantification and interpretation:**
- Whole-slide scanning (40×) with quantification in QuPath or HALO;
- Reported metrics: density of each cell population per mm² of stroma, tumour–stroma interface distance grading (> 100 μm / 20–100 μm / < 20 μm), and the proportion of checkpoint-molecule-positive cells within each population;
- pDC defined by dual markers (CD303⁺CD123⁺) to avoid the specificity limitations of single-marker definitions;
- Blinded reading: two pathologists interpret independently, with an intraclass correlation coefficient (ICC) ≥ 0.75 required for acceptability.

## 4. Target-protein IHC (testing H3)

Priority order (by stability across repeated modelling in this study):

1. **P3H4** (stability 100%) — collagen prolyl 3-hydroxylase, localised to the ER/secretory pathway; expected stromal and tumour-cell cytoplasmic positivity.
2. **RHOU** (15 selections) — atypical Rho GTPase, expected cytoplasmic/membranous positivity, associated with the invasion front.
3. **TNFSF10/TRAIL** (9 selections), **BEX4** (8 selections), **CH25H** (7 selections).
4. Included in the same batch: Ki-67 (proliferation control), Cyclin B1 (H4 hypothesis).

**Scoring**: H-score (0–300; intensity 0–3 × percentage of positive cells), dichotomised by X-tile or by the median; Kaplan–Meier and multivariable Cox analysis (adjusted for age, stage and treatment) against overall survival.

## 5. mRNA-level validation

- RNAscope (ACD) for transcript localisation of P3H4, RHOU and CH25H (FFPE-compatible, single-cell resolution);
- or qRT-PCR on fresh-frozen paired samples (17-gene panel plus five checkpoint genes), with GUSB/HPRT1 as reference genes;
- if clinically feasible, the 17-gene score will be recomputed in this cohort and tested for correlation (Spearman) with the IHC/mIF markers and for association with survival — **the most direct independent test of the risk score of this study**.

## 6. Functional experiments (optional, mechanistic)

- Cell lines: ESCC lines KYSE30, KYSE150 and TE-1 (screened for a P3H4/RHOU expression gradient);
- Interventions: siRNA/shRNA or CRISPRi knockdown of P3H4 and RHOU; overexpression rescue;
- Phenotypes: proliferation (EdU, CCK-8), migration and invasion (Transwell, scratch assay), apoptosis (Annexin V/PI; sensitivity to TNFSF10 treatment);
- Mechanism: ECM density (collagen mapping) and secretome (quantification of the CH25H-related product 25-HC by LC–MS).

## 7. Statistics and sample size

- Primary endpoint: the stromal pDC density difference between the high- and low-risk groups (mIF); secondary endpoint: P3H4 H-score versus overall survival.
- Sample size (cohort A, n = 150, 60 events): with two-sided α = 0.05 and power = 0.8, a log-rank hazard ratio of approximately 0.55 is detectable; with 40 events, a hazard ratio of approximately 0.45 is detectable.
- All analyses will follow the pre-specified plan (the analysis plan frozen in writing before unblinding) to avoid multiple data mining.

## 8. Links to the results of this study

| Output of this study | Role in this plan |
|---|---|
| 719 consensus DEGs | Candidate target pool (source of IHC/RNAscope targets) |
| High-frequency stable genes (P3H4, RHOU, XIST, TNFSF10, BEX4, CH25H) | Priority validation targets |
| 17-gene risk score | Recomputed in the validation cohort and subjected to independent testing |
| Immune-infiltration and checkpoint co-inhibition map | Basis for the mIF panel design |
| Cell-cycle theme | Ki-67 / Cyclin B1 validation module |
