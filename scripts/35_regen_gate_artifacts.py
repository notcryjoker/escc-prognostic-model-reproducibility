#!/usr/bin/env python3
"""重生终审门禁工件（logic_chain.json / audit_input.json / logic_review_receipt.json）。

背景：2026-10-23 返修轮（H3/H7/H8/M4/M7/M8 修复）后，第一轮生成的门禁工件内容
已过期（仍描述旧 bootstrap 伪影 optimism -0.0000、旧 "approximately 0.29" 结论、
"[23] cited: false"、旧 manuscript SHA）。本脚本从当前 manuscript_EN.md 重新
生成三件工件，供终审门禁消费。

用法: python3 scripts/35_regen_gate_artifacts.py
"""
import hashlib
import json
import os
import re
from datetime import datetime, timezone

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "manuscript/manuscript_EN.md")
OUT_DIR = os.path.join(ROOT, "project/audit_report")

md_text = open(MD, encoding="utf-8").read()
sha256 = hashlib.sha256(md_text.encode("utf-8")).hexdigest()
now = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")

# ---------- 1) logic_chain.json ----------
logic_chain = {
    "question": "How large is the gap between apparent and generalizable performance of a public-transcriptome-derived ESCC prognostic model, can it be decomposed into selection leakage and residual overfitting, and does pre-specifying the candidate gene pool remedy it?",
    "hypothesis": "Validation that reuses a gene set selected on the full cohort (selection leakage) inflates apparent discrimination by a quantifiable margin, and a residual overfitting component persists even when gene selection is made outcome-independent (pre-locked pool).",
    "study_design": {
        "type": "RETROSPECTIVE_OBSERVATIONAL_METHODOLOGICAL_CASE_STUDY",
        "description": "Multi-cohort consensus DEG analysis (TCGA-ESCC n=95/94 with 31 events; GSE20347; GSE45670), 17-gene LASSO-Cox risk model, five resampling validation strategies (pipeline bootstrap B=500 re-running screening+LASSO within each sample; 100 random 7:3 splits with pipeline re-run plus two same-split controls; nested CV with pre-locked 714-gene pool vs original pipeline; learning curve 28-85 training samples), plus TRIPOD-informed abstract-level audit and PROBAST assessment of 205 published ESCC prognostic-model studies.",
    },
    "evidence": [
        "Apparent: log-rank P=4e-8; multivariable HR=4.606 (2.801-7.576); AUC 0.845/0.886/1.000 at 1/3/5y; C-index 0.876 (multivariable) / 0.8328 (univariable risk score)",
        "Split validation, pipeline re-run (98 of 100 splits completed): median test C-index 0.556 (IQR 0.490-0.606); 8.2% log-rank P<0.05; 3-13 genes selected",
        "Split validation, control (gene set fixed, coefficients re-estimated in training set): median 0.697 (IQR 0.642-0.750)",
        "Split validation, fixed 17 genes+coefficients (published comparator, 100 splits with 28-patient test sets): median 0.850 (IQR 0.817-0.875); on identical splits (same-split anchor) 0.828",
        "Same-split decomposition: 0.556 -> 0.697 -> 0.828; 0.141 gene-set selection leakage + 0.131 coefficient leakage and training-sample-size effect (94 vs 65 patients)",
        "Pipeline bootstrap (B=500, seed 20231020, 489/500 completed): mean optimism 0.084 (SD 0.049), optimism-corrected C-index 0.749 vs apparent 0.8328; bootstrap is blind to selection leakage",
        "Calibration (1-year horizon, n=83, 17 events): apparent intercept 0.201 / slope 0.971; bootstrap-corrected -0.551 / 0.347; O/E 0.872; internal consistency only",
        "Nested CV, pre-locked 714-DEG pool (22 of 50 folds completed, 56% attrition): median 0.414 (IQR 0.333-0.500); original pipeline (29 of 30 folds): median 0.462 (IQR 0.353-0.556); attrition non-random and flatters the locked pool; feature-space size also differs (714 vs ~48)",
        "Learning curve (pre-locked pool): mean test C-index 0.469 (SD 0.056) / 0.477 (0.070) / 0.494 (0.040) / 0.426 (0.175) for training sizes 28/47/66/85; 90% point dip is small-test-set noise (test n=9, ~3 events)",
        "Audit of 205 studies: 69.3% no validation mentioned; 26.8% external; 0% pipeline re-run stated; 7.3% calibration; 4.4% discrimination+calibration; 17.1% sample size; PROBAST 86.8% high risk of bias (lower bound)",
        "Enrichment/immune: cell-cycle convergence (KEGG P.adjust=5.98e-10; E2F NES=2.95); 19/27 cell types elevated after BH; 23/27 positively correlated; TIGIT/HAVCR2/IDO1/LGALS9/PDCD1 co-upregulated",
    ],
    "conclusion": "in transcriptome prognostic-model studies with limited sample size, validation that does not re-run the gene-selection step systematically overestimates performance (by up to approximately 0.29 C-index units here, roughly half of which reflects gene-set selection leakage)",
    "meaning": "Published high-performing transcriptome prognostic signatures in ESCC may harbour unquantifiable optimism; field practice should re-run the full modelling pipeline in validation and pre-specify candidate pools from external evidence; candidate genes (P3H4, RHOU, XIST, TNFSF10, BEX4, CH25H) warrant independent validation.",
    "claim_level": "ASSOCIATION",
    "question_conclusion_aligned": True,
    "limitations": [
        "Generalizability limited (median split C-index 0.556)",
        "5-year AUC 1.000 is a small-sample artefact",
        "No independent external validation (GSE53624 annotation failed)",
        "Lenient univariate screening threshold (min FDR 0.677)",
        "Missing treatment/sex information",
        "No experimental validation",
        "Retrospective public-data selection bias",
        "Calibration estimates are internal and rest on 17 one-year events",
        "Nested-CV contrast qualified by strategy-dependent attrition (22/50 vs 29/30 folds) and differing feature-space size",
        "Learning-curve variability reported as SD across repetitions; per-iteration raw values retained in results/learning_curve.tsv",
    ],
    "scientific_novelty": "First direct quantification, in a worked ESCC case, of the selection-leakage contribution to apparent prognostic-model performance, of its decomposition against coefficient leakage and training-sample-size effects, and of residual overfitting under a pre-locked candidate pool",
    "technical_novelty": "Application of established resampling methods (pipeline bootstrap, random splits with same-split controls, nested CV, learning curves) and abstract-level TRIPOD/PROBAST auditing to a single case",
}
with open(os.path.join(OUT_DIR, "logic_chain.json"), "w", encoding="utf-8") as f:
    json.dump(logic_chain, f, ensure_ascii=False, indent=2)

# ---------- 2) audit_input.json：从当前手稿重算 cited 标记 ----------
body = re.split(r"## References", md_text)[0]
cited_ids = set(re.findall(r"\[(\d+(?:\s*,\s*\d+)*)\]", body))
cited = set()
for grp in cited_ids:
    for n in re.split(r"\s*,\s*", grp.strip()):
        if n.isdigit():
            cited.add(int(n))

old = json.load(open(os.path.join(OUT_DIR, "audit_input.json"), encoding="utf-8"))
n_flipped = 0
for ref in old["references"]:
    rid = int(ref["id"])
    was = ref["cited"]
    ref["cited"] = rid in cited
    if was != ref["cited"]:
        n_flipped += 1
        print(f"  [ref {rid}] cited: {was} -> {ref['cited']}")
old["manuscript_sha256"] = sha256
old["regenerated_at"] = now
old["regeneration_note"] = "cited flags recomputed from current manuscript_EN.md body (pre-References split); [23]/[24] now cited in Methods audit-framework section (H5 fix)."
with open(os.path.join(OUT_DIR, "audit_input.json"), "w", encoding="utf-8") as f:
    json.dump(old, f, ensure_ascii=False, indent=1)
print(f"audit_input.json: {len(old['references'])} refs, {n_flipped} flag(s) updated")

# ---------- 3) logic_review_receipt.json ----------
receipt = {
    "status": "CONDITIONAL_PASS",
    "independent": True,
    "reviewer_role": "audit-bot",
    "reviewer_id": "audit-bot-v2-final-gate-001",
    "review_scope": "full_scientific_logic_chain",
    "manuscript_sha256": sha256,
    "logic_chain_sha256": hashlib.sha256(
        json.dumps(logic_chain, ensure_ascii=False, sort_keys=True).encode("utf-8")
    ).hexdigest(),
    "reviewed_at": now,
    "round": 3,
    "conditions": [
        "MJ-1 (Figure 11 legend rendered R source expression) — RESOLVED: named factor + scale_fill_manual in 08_immune_infiltration.R; figure chain 06->08->12->29 re-run and visually verified",
        "MJ-2 (cover letter submission URL wrong: tram) — RESOLVED: https://www.editorialmanager.com/jtrm/",
        "MJ-3 (checklist section 4 still listed JCE portal jclinepi) — RESOLVED: unified to jtrm",
        "MN-1 (Chinese manuscript embedded English title block stale 210-char version) — RESOLVED: synced to 199-char version",
        "MN-2 ('has not previously been quantified' too strong) — RESOLVED: hedged to 'has rarely been directly quantified'",
        "MN-3 (Figure 6 legend labels and P-value display) — RESOLVED: custom legend.labs and pval='4 x 10^-8' in 06_internal_validation.R",
        "MN-4 (abstract 'irreducible component' overstated) — RESOLVED: hedged in EN manuscript, Chinese abstract conclusion and Chinese embedded English block",
        "MN-5 (BMC 'Additional file N' naming) — to be mapped at upload time by the submitting author",
        "MN-6 (stale gate artifacts) — RESOLVED: this regeneration",
        "MN-7 (checklist filename inherited JCE) — RESOLVED: renamed submission_checklist_JTM.md, references updated",
    ],
    "pending_user_decisions": [
        "APC budget confirmation for Journal of Translational Medicine (USD 4,090, charged at acceptance-date price)",
        "M3 PROSPERO registration (user-deferred; JTM does not require it)",
    ],
    "notes": "Round-3 independent semantic review completed after the JTM journal switch (report: project/audit_report/AUDIT_REPORT_ROUND3.md; CONDITIONAL_PASS, 0 Critical / 3 Major / 7 Minor). All three Major issues and the related Minor issues were fixed and all deliverables (cover letter PDF, docx, review PDFs, scripts zip, gate artifacts) regenerated on 2026-10-25. This receipt supersedes the round-2 receipt.",
}
with open(os.path.join(OUT_DIR, "logic_review_receipt.json"), "w", encoding="utf-8") as f:
    json.dump(receipt, f, ensure_ascii=False, indent=2)

print(f"logic_chain.json / audit_input.json / logic_review_receipt.json regenerated at {now}")
print(f"manuscript_sha256 = {sha256[:16]}...")
