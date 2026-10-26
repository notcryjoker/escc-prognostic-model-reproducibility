#!/usr/bin/env python3
"""
TRIPOD 报告规范导向的审计（基于标题+摘要）
定位：评估"研究报告了什么"(reporting practices)，而非"研究做了什么"。
      这是无需全文即可严格完成的定位，且比 PROBAST 偏倚风险评估更贴合数据来源。

评估条目（选取可从摘要可靠判定的 TRIPOD 相关条目）：
  T1  是否报告验证类型（外部/内部/无）            -- TRIPOD 16
  T2  是否在验证中重跑包括基因筛选的建模流程        -- TRIPOD 16 / 分析可重复性
  T3  是否报告判别效能(AUC/C-index)               -- TRIPOD 13b
  T4  是否报告校准(calibration)                   -- TRIPOD 13c
  T5  是否报告完整模型(系数/公式/nomogram)         -- TRIPOD 15
  T6  是否报告样本量与事件数                       -- TRIPOD 13a / 17
  T7  是否报告多中心/外部数据来源                  -- TRIPOD 4
输出：条目报告率 + PRISMA 计数 + 结果表
"""
import csv, re, os, json

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OUT = os.path.join(ROOT, "results", "literature_audit")
rows = list(csv.DictReader(open(f"{OUT}/audit_records.csv", encoding="utf-8")))
mod = [r for r in rows if r["is_model_study"] == "1" and r["has_abstract"] == "1"]
n = len(mod)
print(f"纳入 {n} 篇（检索命中 {len(rows)} 篇）\n")

P_EXT = [r"external validation", r"externally validated", r"external cohort",
         r"validation cohort", r"independent cohort", r"independent dataset",
         r"independent validation", r"validation set", r"independent queue"]
P_INT = [r"internal validation", r"training set", r"testing set", r"test set",
         r"train set", r"randomly (divided|split)", r"(divided|split) into a train",
         r"cross-validation", r"cross validation"]
P_RERUN = [r"re-r[au]n", r"rerunning", r"entire pipeline", r"whole pipeline",
           r"re-perform", r"nested cross-validation", r"nested cv"]
P_CAL = [r"calibration", r"hosmer-lemeshow", r"calibration (curve|plot|slope)"]
P_COEF = [r"coefficient", r"risk score\s*=", r"formula", r"nomogram",
          r"prognostic (index|score)\s*=", r"risk model\s*="]
P_EVT = [r"\b\d+\s+(?:deaths|events)\b", r"events?\s*(?:=|:)\s*\d+",
         r"death[s]? occurred in", r"follow-up"]
P_MULTI = [r"multi[- ]?center", r"multicenter", r"multi[- ]?institution"]

res = []
for r in mod:
    t = (r["title"] + " " + r["abstract"]).lower()
    def prob(pats):
        return int(any(re.search(p, t) for p in pats))
    ext, intv = prob(P_EXT), prob(P_INT)
    has_n = r["max_n"] != ""
    res.append({
        "pmid": r["pmid"], "year": r["year"],
        "T1_ext": ext, "T1_int": intv,
        "T1_none": int(not ext and not intv),
        "T2_rerun": prob(P_RERUN),
        "T3_discrimination": int(r["reports_auc"] == "1"),
        "T4_calibration": prob(P_CAL),
        "T5_full_model": prob(P_COEF),
        "T6_sample_size": int(has_n),
        "T7_multisource": prob(P_MULTI),
    })

with open(f"{OUT}/tripod_audit.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(res[0].keys()))
    w.writeheader(); w.writerows(res)

print("===== TRIPOD 条目报告率（基于摘要，n = %d）=====" % n)
items = [("T1_ext", "报告独立外部验证", None),
         ("T1_int", "报告内部拆分/交叉验证", None),
         ("T1_none", "未提及任何验证", None),
         ("T2_rerun", "明确重跑建模流程", None),
         ("T3_discrimination", "报告判别效能(AUC/C-index)", None),
         ("T4_calibration", "报告校准(calibration)", None),
         ("T5_full_model", "报告完整模型(系数/公式)", None),
         ("T6_sample_size", "报告样本量", None),
         ("T7_multisource", "报告多中心/多来源", None)]
out = {}
for key, lab, _ in items:
    v = sum(x[key] for x in res)
    out[lab] = {"n": v, "pct": round(100 * v / n, 1)}
    print(f"  {lab:28s} {v:4d}  ({100*v/n:5.1f}%)")

# 组合指标
both = sum(1 for x in res if x["T3_discrimination"] and x["T4_calibration"])
print(f"\n  同时报告判别与校准: {both} ({100*both/n:.1f}%)")
fullmod = sum(1 for x in res if x["T5_full_model"] and x["T6_sample_size"])
print(f"  同时报告模型与样本量: {fullmod} ({100*fullmod/n:.1f}%)")

# 年份趋势
yrs = {}
for x in res:
    if x["year"]:
        yrs[x["year"]] = yrs.get(x["year"], 0) + 1
print("\n年份分布:", dict(sorted(yrs.items())))

# PRISMA 计数
prisma = {"records_identified": len(rows),
          "records_no_abstract": sum(1 for r in rows if r["has_abstract"] == "0"),
          "records_not_model_study": sum(1 for r in rows if r["is_model_study"] == "0"),
          "included": n}
print("\nPRISMA:", prisma)

with open(f"{OUT}/tripod_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== TRIPOD 报告规范导向审计 =====\n\n")
    f.write(f"PRISMA: 检索 {prisma['records_identified']} 篇; "
            f"排除非预后建模 {prisma['records_not_model_study']} 篇; "
            f"无摘要 {prisma['records_no_abstract']} 篇; 纳入 {n} 篇\n\n")
    f.write("条目报告率:\n")
    for lab, v in out.items():
        f.write(f"  {lab}: {v['n']} ({v['pct']}%)\n")
    f.write(f"\n同时报告判别与校准: {both} ({100*both/n:.1f}%)\n")
    f.write(f"同时报告模型与样本量: {fullmod} ({100*fullmod/n:.1f}%)\n")
    f.write("\n年份分布: " + str(dict(sorted(yrs.items()))) + "\n")
print(f"\n已保存 -> tripod_audit.csv / tripod_summary.txt")
