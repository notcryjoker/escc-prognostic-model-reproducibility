#!/usr/bin/env python3
"""
生成补充材料 Additional file 1/2：文献审计数据集（205 条记录 + 提取的分类）与列说明。

数据来源（两者均已通过独立复算验证）：
  - results/literature_audit/tripod_audit.csv     —— 最终 TRIPOD 条目分类（script 21，权威版本）
  - results/literature_audit/audit_records.csv    —— 文献元数据（journal / title / max_auc / max_n）

输出：
  - manuscript/supplementary/Additional file 1.tsv   数据表
  - manuscript/supplementary/Additional file 2.md               列说明与规则注释

合并键为 PMID；205/205 全部匹配。
"""
import csv
import os

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

with open(os.path.join(ROOT, "results/literature_audit/tripod_audit.csv")) as f:
    tripod = list(csv.DictReader(f))
with open(os.path.join(ROOT, "results/literature_audit/audit_records.csv")) as f:
    records = {r["pmid"]: r for r in csv.DictReader(f)}

assert len(tripod) == 205, f"预期 205 条记录，实际 {len(tripod)}"

# 列名映射：内部代码 -> 可读列名（0 = 否，1 = 是）
COL_MAP = [
    ("T1_ext",             "T1a_External_validation_reported"),
    ("T1_int",             "T1b_Internal_split_or_CV_reported"),
    ("T1_none",            "T1c_No_validation_mentioned"),
    ("T2_rerun",           "T2_Pipeline_rerun_stated"),
    ("T3_discrimination",  "T3_Discrimination_AUC_Cindex_reported"),
    ("T4_calibration",     "T4_Calibration_reported"),
    ("T5_full_model",      "T5_Full_model_reported"),
    ("T6_sample_size",     "T6_Sample_size_reported"),
    ("T7_multisource",     "T7_Multicentre_or_multisource"),
]

out_rows = []
for r in tripod:
    meta = records.get(r["pmid"], {})
    row = {
        "PMID": r["pmid"],
        "Year": r["year"],
        "Journal": meta.get("journal", ""),
        "Title": meta.get("title", "").replace("\t", " ").replace("\n", " ").strip(),
    }
    for src, dst in COL_MAP:
        row[dst] = "1" if r[src] == "1" else "0"
    auc = meta.get("max_auc", "")
    row["Max_reported_AUC_or_Cindex"] = auc if auc else ""
    row["Max_sample_size"] = meta.get("max_n", "") if meta.get("max_n", "") else ""
    out_rows.append(row)

# 按年份 + PMID 排序，便于查阅
out_rows.sort(key=lambda x: (x["Year"], x["PMID"]))

fieldnames = ["PMID", "Year", "Journal", "Title"] + [d for _, d in COL_MAP] + \
             ["Max_reported_AUC_or_Cindex", "Max_sample_size"]

tsv_path = os.path.join(ROOT, "manuscript/supplementary/Additional file 1.tsv")
with open(tsv_path, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                       extrasaction="raise", lineterminator="\n")
    w.writeheader()
    w.writerows(out_rows)

# ---- 自检：写出的文件必须复算出与手稿一致的数字 ----
n = len(out_rows)
checks = {
    "外部验证": sum(int(r["T1a_External_validation_reported"]) for r in out_rows),
    "内部验证": sum(int(r["T1b_Internal_split_or_CV_reported"]) for r in out_rows),
    "无验证": sum(int(r["T1c_No_validation_mentioned"]) for r in out_rows),
    "重跑流程": sum(int(r["T2_Pipeline_rerun_stated"]) for r in out_rows),
    "判别效能": sum(int(r["T3_Discrimination_AUC_Cindex_reported"]) for r in out_rows),
    "校准": sum(int(r["T4_Calibration_reported"]) for r in out_rows),
    "完整模型": sum(int(r["T5_Full_model_reported"]) for r in out_rows),
    "样本量": sum(int(r["T6_Sample_size_reported"]) for r in out_rows),
}
expected = {"外部验证": 55, "内部验证": 22, "无验证": 142, "重跑流程": 0,
            "判别效能": 81, "校准": 15, "完整模型": 85, "样本量": 35}
for k, v in checks.items():
    status = "OK" if v == expected[k] else "!!! 不一致 !!!"
    print(f"{k}: {v}/{n} ({v/n*100:.1f}%)  预期 {expected[k]}  {status}")
    assert v == expected[k], f"{k} 复算失败"

print(f"\nAdditional file 1 已生成: {tsv_path}  ({n} 行 x {len(fieldnames)} 列)")

# ---- 列说明文件 ----
notes = """# Additional file 2 — Literature-audit dataset: notes and column dictionary

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
"""
notes_path = os.path.join(ROOT, "manuscript/supplementary/Additional file 2.md")
with open(notes_path, "w") as f:
    f.write(notes)
print(f"列说明已生成: {notes_path}")
