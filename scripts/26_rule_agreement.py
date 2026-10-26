#!/usr/bin/env python3
"""
规则间一致性验证（回应"非双人独立筛选"的质疑）：
  用第二套独立定义的关键词规则重新分类 205 篇文献，
  与主分析规则（脚本 21）逐条目计算 Cohen's κ 与观察一致率。

设计原则：第二套规则由不同的关键词组合定义，边界与主规则不同
（例如对"外部验证"采用更严格的定义，不含 "validation set"/"validation cohort"
 这类边界模糊的词），以体现规则敏感性对结论的影响。
"""
import csv, re, os
from collections import Counter

OUT = "/workspace/escc_bioinformatics/results/literature_audit"
rows = list(csv.DictReader(open(f"{OUT}/audit_records.csv", encoding="utf-8")))
mod = [r for r in rows if r["is_model_study"] == "1" and r["has_abstract"] == "1"]
n = len(mod)
print(f"纳入 {n} 篇\n")


def probe(t, pats):
    return int(any(re.search(p, t) for p in pats))


# ===== 主规则（与脚本 21 一致，作为参考） =====
R1_EXT = [r"external validation", r"externally validated", r"external cohort",
          r"validation cohort", r"independent cohort", r"independent dataset",
          r"validation set", r"independent queue", r"external dataset"]
R1_INT = [r"internal validation", r"training set", r"testing set", r"test set",
          r"train set", r"randomly (divided|split)", r"(divided|split) into a train",
          r"cross-validation", r"cross validation"]
R1_RERUN = [r"re-r[au]n", r"rerunning", r"entire pipeline", r"whole pipeline",
            r"re-perform", r"nested cross-validation", r"nested cv"]
R1_CAL = [r"calibration", r"hosmer-lemeshow"]
R1_COEF = [r"coefficient", r"risk score\s*=", r"formula", r"nomogram",
           r"prognostic (index|score)\s*=", r"risk model\s*="]
R1_EVT = [r"\b\d+\s+(?:deaths|events)\b", r"events?\s*(?:=|:)\s*\d+",
          r"death[s]? occurred in", r"follow-up"]
R1_N = [r"\bn\s*=\s*\d{2,4}\b", r"\d{2,4}\s+(?:patients|cases|samples)"]

# ===== 第二套规则（独立定义，不同边界） =====
# 外部验证：严格定义，要求出现 "external" 或 "independent" 明确修饰
R2_EXT = [r"external (validation|cohort|dataset|independent)",
          r"externally validated", r"independent (validation|cohort|dataset)",
          r"validated in an independent", r"independent external"]
# 内部验证：要求同时出现 train 与 test，或明确的 CV 表述
R2_INT = [r"training (set|cohort).{0,40}test", r"test (set|cohort).{0,40}train",
          r"(k-fold|leave-one-out|10-fold|5-fold).{0,20}(cross-?validation|cv)",
          r"cross-?validation", r"randomly (divided|split|partitioned)"]
R2_RERUN = [r"re-r[au]n", r"re-?perform.{0,30}(selection|screen|pipeline)",
            r"entire (modelling|pipeline)", r"nested cross-?validation"]
R2_CAL = [r"calibration (curve|plot|slope|test)", r"hosmer-lemeshow",
          r"calibration was (assessed|evaluated|performed)"]
R2_COEF = [r"regression coefficient", r"risk score\s*=\s*[\d.+\- ]",
           r"nomogram", r"prognostic index\s*=", r"formula\s*="]
R2_N = [r"\bn\s*=\s*\d{2,4}\b", r"cohort of \d{2,4}",
        r"\d{2,4}\s+(?:patients|cases|samples|tumors|specimens)"]


def classify(text, ext_r, int_r, rerun_r, cal_r, coef_r, n_r):
    t = text.lower()
    ext = probe(t, ext_r)
    intv = probe(t, int_r)
    rerun = probe(t, rerun_r)
    cal = probe(t, cal_r)
    coef = probe(t, coef_r)
    ns = [int(x) for x in re.findall(r"n\s*=\s*(\d{2,4})", t)]
    ns += [int(x) for x in re.findall(r"(\d{2,4})\s+(?:patients|cases|samples|tumors|specimens)", t)]
    has_n = any(20 <= v <= 5000 for v in ns)
    return ext, intv, rerun, cal, coef, int(has_n)


def kappa(a, b):
    n = len(a)
    if n == 0:
        return float("nan")
    po = sum(1 for i in range(n) if a[i] == b[i]) / n
    pa = sum(a) / n
    pb = sum(b) / n
    pe = pa * pb + (1 - pa) * (1 - pb)
    if (1 - pe) < 1e-9:
        return 1.0 if po == 1.0 else 0.0
    return (po - pe) / (1 - pe)


# 双规则分类
R1, R2 = [], []
for r in mod:
    text = r["title"] + " " + r["abstract"]
    R1.append(classify(text, R1_EXT, R1_INT, R1_RERUN, R1_CAL, R1_COEF, R1_N))
    R2.append(classify(text, R2_EXT, R2_INT, R2_RERUN, R2_CAL, R2_COEF, R2_N))

items = ["external_validation", "internal_validation", "rerun_pipeline",
         "calibration", "full_model", "sample_size"]

print("===== 规则间一致性（Cohen's κ）=====")
print(f"{'条目':<22} {'R1 阳性':<10} {'R2 阳性':<10} {'观察一致率':<12} {'κ':<8}")
print("-" * 64)
res = {}
for k, idx in zip(items, range(6)):
    a = [r1[idx] for r1 in R1]
    b = [r2[idx] for r2 in R2]
    k_val = kappa(a, b)
    po = sum(1 for i in range(n) if a[i] == b[i]) / n
    res[k] = {"R1_n": sum(a), "R2_n": sum(b), "po": po, "kappa": k_val}
    lab = {"external_validation": "外部验证", "internal_validation": "内部验证",
           "rerun_pipeline": "重跑建模流程", "calibration": "报告校准",
           "full_model": "报告完整模型", "sample_size": "报告样本量"}[k]
    print(f"{lab:<22} {sum(a):<10} {sum(b):<10} {po:<12.3f} {k_val:<8.3f}")

# 核心结论的稳健性：未提及任何验证的比例
none1 = sum(1 for i in range(n) if R1[i][0] == 0 and R1[i][1] == 0)
none2 = sum(1 for i in range(n) if R2[i][0] == 0 and R2[i][1] == 0)
a_none = [int(R1[i][0] == 0 and R1[i][1] == 0) for i in range(n)]
b_none = [int(R2[i][0] == 0 and R2[i][1] == 0) for i in range(n)]
k_none = kappa(a_none, b_none)
print(f"\n--- 核心结论稳健性 ---")
print(f"未提及任何验证: R1={none1} ({100*none1/n:.1f}%)  R2={none2} ({100*none2/n:.1f}%)")
print(f"  κ = {k_none:.3f}  观察一致率 = {sum(1 for i in range(n) if a_none[i]==b_none[i])/n:.3f}")

# 总体判定
mean_k = sum(res[k]["kappa"] for k in items) / len(items)
print(f"\n六条目平均 κ = {mean_k:.3f}")
if mean_k >= 0.81:
    print("→ 几乎完美一致（almost perfect, κ ≥ 0.81）")
elif mean_k >= 0.61:
    print("→ 实质一致（substantial, 0.61 ≤ κ < 0.81）")
elif mean_k >= 0.41:
    print("→ 中等一致（moderate, 0.41 ≤ κ < 0.61）")

# 保存
with open(f"{OUT}/rule_agreement.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["item", "R1_positive", "R2_positive", "observed_agreement", "kappa"])
    for k in items:
        w.writerow([k, res[k]["R1_n"], res[k]["R2_n"],
                    round(res[k]["po"], 4), round(res[k]["kappa"], 4)])
    w.writerow(["any_validation_none", none1, none2,
                round(sum(1 for i in range(n) if a_none[i]==b_none[i])/n, 4),
                round(k_none, 4)])

with open(f"{OUT}/rule_agreement_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== 规则间一致性验证（Cohen's κ）=====\n\n")
    f.write(f"两套独立关键词规则对 {n} 篇文献分别分类\n\n")
    f.write(f"{'条目':<22} {'R1':<8} {'R2':<8} {'一致率':<8} {'κ':<8}\n")
    for k in items:
        lab = {"external_validation": "外部验证", "internal_validation": "内部验证",
               "rerun_pipeline": "重跑流程", "calibration": "校准",
               "full_model": "完整模型", "sample_size": "样本量"}[k]
        f.write(f"{lab:<22} {res[k]['R1_n']:<8} {res[k]['R2_n']:<8} "
                f"{res[k]['po']:<8.3f} {res[k]['kappa']:<8.3f}\n")
    f.write(f"\n未提及任何验证: R1={none1} ({100*none1/n:.1f}%) "
            f"R2={none2} ({100*none2/n:.1f}%) κ={k_none:.3f}\n")
    f.write(f"\n六条目平均 κ = {mean_k:.3f}\n")
print(f"\n已保存 rule_agreement.csv / rule_agreement_summary.txt")
