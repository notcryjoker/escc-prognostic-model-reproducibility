#!/usr/bin/env python3
"""
PROBAST 偏倚风险评估（基于已获取的全文）
评估 4 个域中可从全文自动判定的关键信号问题，并给出总体偏倚风险判定。

PROBAST 判定规则（Wolff et al. 2019）：
  总体高风险  = 任一域为高风险
  总体低风险  = 所有域均为低风险
  总体不清楚  = 至少一域不清楚且无高风险域

关键条目（Analysis 域）：
  4.1 事件数/候选变量比 (EPV) 是否足够 (>=10)
  4.5 是否避免仅依据单变量分析选择预测因子
  4.7 是否报告模型性能(判别+校准)
  4.8 是否进行了恰当验证(外部验证 vs 仅内部/无)
"""
import csv, re, os, json
from collections import Counter

OUT = "/workspace/escc_bioinformatics/results/literature_audit"
FT = f"{OUT}/fulltext"

recs = list(csv.DictReader(open(f"{OUT}/audit_records.csv", encoding="utf-8")))
mod = [r for r in recs if r["is_model_study"] == "1" and r["has_abstract"] == "1"]

# ---------- 载入全文 ----------
ft = {}
if os.path.isdir(FT):
    for fn in os.listdir(FT):
        if fn.endswith(".txt"):
            p = fn.replace(".txt", "")
            try:
                t = open(f"{FT}/{fn}", encoding="utf-8", errors="replace").read()
                if len(t) > 3000:
                    ft[p] = t.lower()
            except Exception:
                pass
print(f"可用全文: {len(ft)} 篇", flush=True)

P_EXT = [r"external validation", r"externally validated", r"external cohort",
         r"validation cohort", r"independent cohort", r"independent dataset",
         r"validation set", r"independent queue", r"external dataset"]
P_INT = [r"internal validation", r"training set", r"testing set", r"test set",
         r"train set", r"randomly (divided|split)", r"(divided|split) into a train",
         r"cross-validation", r"cross validation"]
P_RERUN = [r"re-r[au]n", r"rerunning", r"entire pipeline", r"whole pipeline",
           r"nested cross-validation", r"nested cv"]
P_CAL = [r"calibration", r"hosmer-lemeshow"]
P_DISCR = [r"\bauc\b", r"c-index", r"concordance index", r"roc"]
P_UNI = [r"univariate cox", r"univariable cox", r"univariate analysis",
         r"univariable analysis", r"univariate regression",
         r"univariate logistic"]
P_PENAL = [r"lasso", r"ridge", r"elastic net", r"penaliz", r"penalis",
           r"shrinkage", r"overfitting", r"over-fitting"]
P_MISS = [r"missing data", r"imputation", r"complete case", r"missing values"]
P_EVT = [r"\b(\d{2,4})\s+(?:deaths?|events?)\b",
         r"(?:deaths?|events?)\s*(?:=|:)\s*(\d{2,4})",
         r"(\d{2,4})\s+patients?\s+(?:died|experienced)"]
P_N = [r"\bn\s*=\s*(\d{2,4})\b", r"(\d{2,4})\s+(?:patients|cases|samples|tumors)"]


def probe(t, pats):
    return any(re.search(p, t) for p in pats)


rows = []
for r in mod:
    pmid = r["pmid"]
    has_ft = pmid in ft
    t = ft.get(pmid, (r["title"] + " " + r["abstract"]).lower())
    src = "fulltext" if has_ft else "abstract"

    ext = probe(t, P_EXT)
    intv = probe(t, P_INT)
    rerun = probe(t, P_RERUN)
    cal = probe(t, P_CAL)
    disc = probe(t, P_DISCR)
    uni = probe(t, P_UNI)
    penal = probe(t, P_PENAL)
    miss = probe(t, P_MISS)

    # 事件数与样本量
    ev = [int(x) for x in re.findall(P_EVT[0], t)]
    ev += [int(x) for x in re.findall(P_EVT[1], t)]
    ev = [e for e in ev if 5 <= e <= 2000]
    ns = [int(x) for x in re.findall(P_N[0], t)]
    ns += [int(x) for x in re.findall(P_N[1], t)]
    ns = [n for n in ns if 20 <= n <= 5000]
    max_ev = max(ev) if ev else ""
    max_n = max(ns) if ns else ""

    # ---- PROBAST Analysis 域条目判定 ----
    # 4.8 验证: 外部+重跑 = 低; 仅外部 = 不清楚; 仅内部/无 = 高
    if ext and rerun:
        v48, rob48 = "Low", "Low"
    elif ext:
        v48, rob48 = "Unclear", "Unclear"
    elif intv or not (ext or intv):
        v48, rob48 = "High", "High"

    # 4.7 性能报告: 判别+校准 = 低; 仅判别 = 不清楚; 均无 = 高
    if disc and cal:
        rob47 = "Low"
    elif disc or cal:
        rob47 = "Unclear"
    else:
        rob47 = "High"

    # 4.5 单变量筛选: 未仅用单变量 = 低; 用单变量但未提及惩罚 = 高
    if not uni:
        rob45 = "Low"
    elif penal:
        rob45 = "Unclear"
    else:
        rob45 = "High"

    # 4.1 EPV: 若有事件数与变量数(近似用 10) 判定
    if max_ev != "" and isinstance(max_ev, int):
        epv = max_ev / 10.0     # 假定最终模型约 10 个变量（保守）
        rob41 = "Low" if epv >= 10 else ("Unclear" if epv >= 5 else "High")
        epv_val = round(epv, 1)
    else:
        rob41, epv_val = "Unclear", ""

    # 4.4 缺失数据: 提及 = 低; 未提及 = 不清楚
    rob44 = "Low" if miss else "Unclear"

    domain_risks = [rob48, rob47, rob45, rob41, rob44]
    if "High" in domain_risks:
        overall = "High"
    elif "Unclear" in domain_risks:
        overall = "Unclear"
    else:
        overall = "Low"

    rows.append({
        "pmid": pmid, "source": src, "year": r["year"],
        "external_validation": int(ext), "internal_validation": int(intv),
        "rerun_pipeline": int(rerun),
        "discrimination": int(disc), "calibration": int(cal),
        "univariate_screening": int(uni), "penalization": int(penal),
        "missing_data": int(miss),
        "max_n": max_n, "max_events": max_ev, "epv_approx": epv_val,
        "probast_4_1_epv": rob41, "probast_4_4_missing": rob44,
        "probast_4_5_selection": rob45, "probast_4_7_performance": rob47,
        "probast_4_8_validation": rob48,
        "probast_overall": overall,
    })

with open(f"{OUT}/probst_assessment.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)

n = len(rows)
nft = sum(1 for x in rows if x["source"] == "fulltext")
print(f"\n===== PROBAST 评估汇总 (n = {n}; 其中基于全文 {nft} 篇) =====", flush=True)

c = Counter(x["probast_overall"] for x in rows)
for k in ["Low", "Unclear", "High"]:
    v = c.get(k, 0)
    print(f"  总体偏倚风险 {k:8s}: {v:4d} ({100*v/n:5.1f}%)")

print("\n各条目判定:", flush=True)
for key, lab in [("probast_4_8_validation", "4.8 验证充分性"),
                 ("probast_4_7_performance", "4.7 性能报告(判别+校准)"),
                 ("probast_4_5_selection", "4.5 预测因子选择"),
                 ("probast_4_1_epv", "4.1 事件数充足性"),
                 ("probast_4_4_missing", "4.4 缺失数据处理")]:
    cc = Counter(x[key] for x in rows)
    tot = sum(cc.values())
    hi = cc.get("High", 0)
    print(f"  {lab:24s} 高风险 {hi:4d} ({100*hi/tot:5.1f}%)")

# 全文子集（更可靠）
sub = [x for x in rows if x["source"] == "fulltext"]
if sub:
    print(f"\n[仅全文子集 n = {len(sub)}]", flush=True)
    c2 = Counter(x["probast_overall"] for x in sub)
    for k in ["Low", "Unclear", "High"]:
        v = c2.get(k, 0)
        print(f"  总体 {k:8s}: {v:4d} ({100*v/len(sub):5.1f}%)")
    for key, lab in [("probast_4_8_validation", "4.8 验证"),
                     ("probast_4_7_performance", "4.7 性能"),
                     ("probast_4_5_selection", "4.5 选择")]:
        cc = Counter(x[key] for x in sub)
        print(f"  {lab:12s} 高风险 {cc.get('High',0):4d} "
              f"({100*cc.get('High',0)/len(sub):5.1f}%)")

with open(f"{OUT}/probst_assessment_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== PROBAST 偏倚风险评估汇总 =====\n\n")
    f.write(f"评估 {n} 篇（其中基于全文 {nft} 篇，其余基于摘要）\n\n")
    f.write("[总体偏倚风险]\n")
    for k in ["Low", "Unclear", "High"]:
        v = c.get(k, 0)
        f.write(f"  {k}: {v} ({100*v/n:.1f}%)\n")
    f.write("\n[各条目高风险比例]\n")
    for key, lab in [("probast_4_8_validation", "4.8 验证充分性"),
                     ("probast_4_7_performance", "4.7 性能报告"),
                     ("probast_4_5_selection", "4.5 预测因子选择"),
                     ("probast_4_1_epv", "4.1 事件数充足性"),
                     ("probast_4_4_missing", "4.4 缺失数据处理")]:
        cc = Counter(x[key] for x in rows)
        f.write(f"  {lab}: {cc.get('High',0)}/{sum(cc.values())} "
                f"({100*cc.get('High',0)/sum(cc.values()):.1f}%)\n")
    if sub:
        f.write(f"\n[仅全文子集 n={len(sub)}]\n")
        c2 = Counter(x["probast_overall"] for x in sub)
        for k in ["Low", "Unclear", "High"]:
            f.write(f"  {k}: {c2.get(k,0)} ({100*c2.get(k,0)/len(sub):.1f}%)\n")
print("\n已保存 probst_assessment.csv / probst_assessment_summary.txt", flush=True)
