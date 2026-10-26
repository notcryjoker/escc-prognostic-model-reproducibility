#!/usr/bin/env python3
"""
系统评价升级：
  1. 将纳入的 205 篇文献映射到 PMC，获取开放获取全文
  2. 基于全文执行 PROBAST(精简版) 偏倚风险条目评估
  3. 以全文判定为参考，评估"摘要层面判定"的敏感度/特异度/一致率
输出：PRISMA 与 PROBAST 结果
"""
import urllib.request, urllib.parse, json, csv, re, time, os, sys

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "literature_audit")
os.makedirs(f"{OUT}/fulltext", exist_ok=True)


def fetch(url, tries=4):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research)"})
            with urllib.request.urlopen(req, timeout=180) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"   重试{i+1} {type(e).__name__}", flush=True)
            time.sleep(6)
    return None


recs = list(csv.DictReader(open(f"{OUT}/audit_records.csv")))
recs = [r for r in recs if r["is_model_study"] == "1" and r["has_abstract"] == "1"]
print(f"纳入文献 {len(recs)} 篇", flush=True)

# ---------- 1. PMID -> PMCID ----------
print("\n[1] 映射 PMID -> PMCID ...", flush=True)
pmc_map = {}
pmids = [r["pmid"] for r in recs]
for i in range(0, len(pmids), 100):
    chunk = pmids[i:i + 100]
    u = f"{BASE}/elink.fcgi?dbfrom=pubmed&db=pmc&id={','.join(chunk)}&retmode=json"
    raw = fetch(u)
    if not raw:
        continue
    try:
        j = json.loads(raw)
    except Exception:
        continue
    for ls in j.get("linksets", []):
        pid = str(ls.get("ids", [None])[0])
        for db in ls.get("linksetdbs", []):
            if db.get("dbto") == "pmc":
                for l in db.get("links", []):
                    pmc_map[pid] = str(l.get("id")) if isinstance(l, dict) else str(l)
    print(f"   已映射 {len(pmc_map)}/{len(pmids)}", flush=True)
    time.sleep(3)
print(f"PMC 可获取: {len(pmc_map)} 篇", flush=True)

# ---------- 2. 获取全文 ----------
print("\n[2] 下载 PMC 全文 ...", flush=True)
def xml2text(xml):
    # 去掉 back matter（参考文献）与表格标签，保留正文文本
    xml = re.sub(r"<back>.*?</back>", " ", xml, flags=re.S)
    xml = re.sub(r"<(ref-list|xref|table-wrap|inline-formula|disp-formula)[^>]*>.*?</\1>",
                 " ", xml, flags=re.S)
    xml = re.sub(r"<[^>]+>", " ", xml)
    return re.sub(r"\s+", " ", xml).strip()


full = {}
pmc_ids = sorted(set(pmc_map.values()))
for i in range(0, len(pmc_ids), 50):
    chunk = pmc_ids[i:i + 50]
    u = f"{BASE}/efetch.fcgi?db=pmc&id={','.join(chunk)}&retmode=xml"
    xml = fetch(u)
    if not xml:
        continue
    arts = re.findall(r"<article\b.*?</article>", xml, re.S)
    for a in arts:
        m = re.search(r"<article-id[^>]*pub-id-type=\"pmc\">\s*(\d+)\s*</article-id>", a)
        if not m:
            continue
        full[m.group(1)] = xml2text(a)
    print(f"   已下载 {len(full)}/{len(pmc_ids)}", flush=True)
    time.sleep(3)
print(f"全文获取成功: {len(full)} 篇", flush=True)

# ---------- 3. PROBAST 条目（基于全文） ----------
print("\n[3] PROBAST 条目评估 ...", flush=True)

def probe(t, pats):
    return int(any(re.search(p, t) for p in pats))

P_EXTERNAL = [r"external validation", r"externally validated", r"external cohort",
              r"validation cohort", r"independent cohort", r"independent dataset",
              r"independent validation", r"external independent", r"validation set",
              r"independent external", r"validated in.{0,30}independent"]
P_INTERNAL = [r"internal validation", r"training set", r"testing set", r"test set",
              r"train set", r"randomly (divided|split)", r"were (divided|split) into",
              r"cross-validation", r"cross validation"]
P_RERUN   = [r"re-r[au]n", r"rerunning", r"entire pipeline", r"whole pipeline",
             r"re-perform", r"nested cross-validation", r"nested cv",
             r"repeated.{0,30}(selection|screening)", r"re-ran"]
P_CALIB   = [r"calibration", r"calibration curve", r"calibration plot",
             r"hosmer-lemeshow", r"calibration slope"]
P_COEF    = [r"coefficient", r"risk score\s*=", r"formula", r"nomogram",
             r"risk model\s*=", r"beta coefficient", r"ß"]
P_UNI     = [r"univariate cox", r"univariable cox", r"univariate analysis",
             r"univariable analysis", r"univariate regression"]
P_MULTI   = [r"multivariate cox", r"multivariable cox", r"multivariate analysis",
             r"multivariable regression"]
P_EPMISS  = [r"missing data", r"imputation", r"multiple imputation",
             r"complete case", r"missing values"]

rows = []
for r in recs:
    pmcid = pmc_map.get(r["pmid"], "")
    ft = full.get(pmcid, "")
    abst = (r["title"] + " " + r["abstract"]).lower()
    full_l = ft.lower()

    # 全文判定（若有全文），否则回退到摘要
    src = full_l if ft else abst
    has_ft = int(bool(ft))

    val_ext = probe(src, P_EXTERNAL)
    val_int = probe(src, P_INTERNAL)
    rerun   = probe(src, P_RERUN)
    calib   = probe(src, P_CALIB)
    coef    = probe(src, P_COEF)
    uni     = probe(src, P_UNI)
    multi   = probe(src, P_MULTI)
    miss    = probe(src, P_EPMISS)

    # PROBAST 验证域分级
    if val_ext and rerun:
        rob = "Low"
    elif val_ext:
        rob = "Unclear-External"
    elif val_int:
        rob = "High (internal only)"
    else:
        rob = "High (no validation)"

    # 摘要层面判定（用于与全文对比）
    a_ext = probe(abst, P_EXTERNAL)
    a_int = probe(abst, P_INTERNAL)

    # 事件数/样本量（全文优先）
    ns = [int(x) for x in re.findall(r"\bn\s*=\s*(\d{2,4})\b", src)]
    ns += [int(x) for x in re.findall(r"(\d{2,4})\s+(?:patients|cases|samples)", src)]
    ns = [n for n in ns if 20 <= n <= 5000]

    rows.append({
        "pmid": r["pmid"], "pmcid": pmcid, "has_fulltext": has_ft,
        "year": r["year"],
        "ft_external": val_ext, "ft_internal": val_int, "ft_rerun": rerun,
        "ft_calibration": calib, "ft_coefficients": coef,
        "ft_univariate_screening": uni, "ft_multivariable": multi,
        "ft_missing_data": miss,
        "probst_validation_risk": rob,
        "abst_external": a_ext, "abst_internal": a_int,
        "max_n": max(ns) if ns else "",
    })

with open(f"{OUT}/probst_fulltext.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader(); w.writerows(rows)
print(f"已保存 -> probst_fulltext.csv ({len(rows)} 行)", flush=True)

# ---------- 4. 摘要 vs 全文判定一致性 ----------
print("\n[4] 摘要判定 vs 全文判定一致性（以全文为参考）", flush=True)
sub = [r for r in rows if r["has_fulltext"] == 1]
n = len(sub)
tp = sum(1 for r in sub if r["ft_external"] and r["abst_external"])
fp = sum(1 for r in sub if not r["ft_external"] and r["abst_external"])
fn = sum(1 for r in sub if r["ft_external"] and not r["abst_external"])
tn = sum(1 for r in sub if not r["ft_external"] and not r["abst_external"])
sens = tp / (tp + fn) if (tp + fn) else float("nan")
spec = tn / (tn + fp) if (tn + fp) else float("nan")
agree = (tp + tn) / n if n else float("nan")

print(f"  有全文可用于核对: {n} 篇")
print(f"  真阳性 {tp} | 假阳性 {fp} | 假阴性 {fn} | 真阴性 {tn}")
print(f"  敏感度(摘要识别外部验证的能力): {sens:.3f}")
print(f"  特异度: {spec:.3f}")
print(f"  总一致率: {agree:.3f}")

# ---------- 5. PROBAST 汇总 ----------
print("\n[5] PROBAST 验证域分级分布（全文，n = %d）" % n, flush=True)
from collections import Counter
c = Counter(r["probst_validation_risk"] for r in sub)
for k, v in c.most_common():
    print(f"  {k:24s} {v:4d} ({100*v/n:.1f}%)")

print("\n[6] 其他 PROBAST 相关条目（全文，n = %d）" % n, flush=True)
for key, lab in [("ft_rerun", "明确重跑建模流程"),
                 ("ft_calibration", "报告校准"),
                 ("ft_coefficients", "报告模型系数/公式"),
                 ("ft_univariate_screening", "使用单变量筛选"),
                 ("ft_multivariable", "使用多变量分析"),
                 ("ft_missing_data", "提及缺失数据处理")]:
    v = sum(r[key] for r in sub)
    print(f"  {lab:20s} {v:4d} ({100*v/n:.1f}%)")

with open(f"{OUT}/probst_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== PROBAST(精简) 验证域评估 & 摘要-全文一致性 =====\n\n")
    f.write(f"纳入 {len(rows)} 篇; 有全文 {n} 篇\n\n")
    f.write(f"[一致性] 敏感度 {sens:.3f} | 特异度 {spec:.3f} | 一致率 {agree:.3f}\n")
    f.write(f"  TP {tp} | FP {fp} | FN {fn} | TN {tn}\n\n")
    f.write("[验证域分级]\n")
    for k, v in c.most_common():
        f.write(f"  {k}: {v} ({100*v/n:.1f}%)\n")
    f.write("\n[其他条目]\n")
    for key, lab in [("ft_rerun", "重跑建模流程"), ("ft_calibration", "报告校准"),
                     ("ft_coefficients", "报告系数"), ("ft_univariate_screening", "单变量筛选"),
                     ("ft_multivariable", "多变量分析"), ("ft_missing_data", "缺失数据处理")]:
        v = sum(r[key] for r in sub)
        f.write(f"  {lab}: {v} ({100*v/n:.1f}%)\n")
print(f"\n汇总已保存 -> {OUT}/probst_summary.txt", flush=True)
