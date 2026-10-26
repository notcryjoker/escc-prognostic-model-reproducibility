#!/usr/bin/env python3
"""
文献审计：已发表的食管鳞癌(ESCC)转录组预后模型研究的验证实践
目的：系统评估领域内"验证"的实际做法，统计
  - 独立外部验证的比例
  - 仅内部拆分/交叉验证"验证"的比例
  - 明确说明在验证中重跑基因筛选环节的比例
  - 报告的判别效能(AUC/C-index)与样本量分布
数据源：PubMed E-utilities
"""
import urllib.request, urllib.parse, time, json, re, csv, sys, os

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "literature_audit")
os.makedirs(OUT, exist_ok=True)


def fetch(url):
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research)"})
            with urllib.request.urlopen(req, timeout=120) as r:
                return r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"  重试 {attempt+1}: {type(e).__name__}", flush=True)
            time.sleep(6)
    return None


# ---------- 1. 检索 ----------
TERM = ('("esophageal squamous cell carcinoma"[Title/Abstract] OR "oesophageal squamous cell carcinoma"[Title/Abstract] '
        'OR "ESCC"[Title/Abstract]) '
        'AND ("prognostic signature"[Title/Abstract] OR "prognostic model"[Title/Abstract] '
        'OR "prognostic gene"[Title/Abstract] OR "prognosis-related gene"[Title/Abstract] '
        'OR "risk score"[Title/Abstract] OR "risk model"[Title/Abstract] '
        'OR "gene signature"[Title/Abstract] OR "nomogram"[Title/Abstract]) '
        'AND ("TCGA"[Title/Abstract] OR "GEO"[Title/Abstract] OR "Gene Expression Omnibus"[Title/Abstract] '
        'OR "bioinformatics"[Title/Abstract] OR "transcriptome"[Title/Abstract] '
        'OR "gene expression"[Title/Abstract] OR "RNA-seq"[Title/Abstract]) '
        'AND 2015:2026[pdat] AND English[Language]')

print("检索 PubMed ...", flush=True)
u = f"{BASE}/esearch.fcgi?db=pubmed&term={urllib.parse.quote(TERM)}&retmax=1000&retmode=json"
raw = fetch(u)
pmids = json.loads(raw)["esearchresult"]["idlist"]
print(f"命中 {len(pmids)} 篇", flush=True)
time.sleep(3)

# ---------- 2. 获取摘要 ----------
print("下载摘要 ...", flush=True)
records = []
BATCH = 100
for i in range(0, len(pmids), BATCH):
    chunk = pmids[i:i + BATCH]
    u = (f"{BASE}/efetch.fcgi?db=pubmed&id={','.join(chunk)}"
         f"&retmode=xml&rettype=abstract")
    xml = fetch(u)
    if not xml:
        continue
    arts = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml, re.S)
    for a in arts:
        pmid = re.search(r"<PMID[^>]*>(\d+)</PMID>", a)
        title = re.search(r"<ArticleTitle>(.*?)</ArticleTitle>", a, re.S)
        abst = re.search(r"<Abstract>(.*?)</Abstract>", a, re.S)
        year = re.search(r"<PubDate>.*?<Year>(\d{4})</Year>", a, re.S)
        jour = re.search(r"<Journal>.*?<Title>(.*?)</Title>", a, re.S)
        if not (pmid and title):
            continue

        def clean(m):
            if not m:
                return ""
            t = re.sub(r"<[^>]+>", " ", m.group(1))
            return re.sub(r"\s+", " ", t).strip()

        records.append({
            "pmid": pmid.group(1),
            "year": year.group(1) if year else "",
            "journal": clean(jour),
            "title": clean(title),
            "abstract": clean(abst),
        })
    print(f"  已获取 {len(records)}/{len(pmids)}", flush=True)
    time.sleep(3)

print(f"共获取 {len(records)} 篇记录", flush=True)

# ---------- 3. 规则化提取 ----------
def has(text, pats):
    t = text.lower()
    return any(re.search(p, t) for p in pats)

EXTERNAL = [r"external validation", r"externally validated", r"external cohort",
            r"independent validation", r"independent external", r"validation cohort",
            r"independent cohort", r"external independent", r"validated in an independent",
            r"independent dataset", r"independent queue"]
INTERNAL = [r"internal validation", r"training set", r"testing set", r"test set",
            r"split", r"cross-validation", r"cross validation", r"train set"]
RERUN   = [r"re-r[au]n", r"rerunning", r"entire pipeline", r"whole pipeline",
           r"re-performed", r"repeat(ed)? the (gene )?screen",
           r"gene selection .{0,40}(within|inside|repeated|rerun)",
           r"nested cross-validation", r"nested cv"]
AUC     = [r"\bauc\b", r"c-index", r"concordance index", r"roc"]
SAMPLE  = [r"\bn\s*=\s*\d{2,4}\b", r"\d{2,4}\s+patients", r"\d{2,4}\s+cases",
           r"\d{2,4}\s+samples", r"cohort of \d{2,4}"]

rows = []
for r in records:
    txt = r["title"] + " " + r["abstract"]
    low = txt.lower()
    # 是否为转录组预后建模研究
    is_model = has(low, [r"prognostic signature", r"prognostic model", r"risk score",
                         r"risk model", r"gene signature", r"prognosis-related",
                         r"prognostic gene", r"nomogram"])
    gse = set(re.findall(r"GSE\d{4,6}", txt))
    # AUC 数值
    aucs = [float("0." + x) for x in
            re.findall(r"(?:auc|c-index|concordance)[^0-9]{0,30}?0\.(\d{2,3})", low)]
    aucs = [a for a in aucs if 0.4 <= a <= 1.0]
    # 样本量
    ns = [int(x) for x in re.findall(r"\bn\s*=\s*(\d{2,4})\b", low)]
    ns += [int(x) for x in re.findall(r"(\d{2,4})\s+(?:patients|cases|samples)", low)]
    ns = [n for n in ns if 20 <= n <= 5000]

    rows.append({
        "pmid": r["pmid"], "year": r["year"], "journal": r["journal"],
        "title": r["title"][:250],
        "abstract": r["abstract"][:4000],
        "is_model_study": int(is_model),
        "n_geo_cohorts": len(gse),
        "external_validation": int(has(low, EXTERNAL)),
        "internal_validation": int(has(low, INTERNAL)),
        "rerun_pipeline": int(has(low, RERUN)),
        "reports_auc": int(has(low, AUC)),
        "max_auc": max(aucs) if aucs else "",
        "max_n": max(ns) if ns else "",
        "has_abstract": int(bool(r["abstract"])),
    })

with open(f"{OUT}/audit_records.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)
print(f"已保存 {len(rows)} 条记录 -> audit_records.csv", flush=True)

# ---------- 4. 汇总统计 ----------
mod = [r for r in rows if r["is_model_study"] and r["has_abstract"]]
n = len(mod)
print(f"\n===== 文献审计汇总 =====")
print(f"检索命中 {len(rows)} 篇；其中预后建模研究且有摘要: {n} 篇\n")

ext = sum(r["external_validation"] for r in mod)
intv = sum(r["internal_validation"] for r in mod)
rerun = sum(r["rerun_pipeline"] for r in mod)
auc_rep = sum(r["reports_auc"] for r in mod)
both = sum(1 for r in mod if r["external_validation"] and r["internal_validation"])
neither = sum(1 for r in mod if not r["external_validation"] and not r["internal_validation"])

print(f"报告独立外部验证:                {ext:4d} ({100*ext/n:.1f}%)")
print(f"报告内部拆分/交叉验证:          {intv:4d} ({100*intv/n:.1f}%)")
print(f"两者均报告:                      {both:4d} ({100*both/n:.1f}%)")
print(f"两者均未提及:                    {neither:4d} ({100*neither/n:.1f}%)")
print(f"明确说明在验证中重跑建模流程:    {rerun:4d} ({100*rerun/n:.1f}%)")
print(f"报告 AUC / C-index:              {auc_rep:4d} ({100*auc_rep/n:.1f}%)")

aucs_all = [float(r["max_auc"]) for r in mod if r["max_auc"] != ""]
if aucs_all:
    aucs_all.sort()
    med = aucs_all[len(aucs_all) // 2]
    print(f"\n报告的最大 AUC/C-index: n={len(aucs_all)}, 中位 {med:.3f}, "
          f"范围 {min(aucs_all):.3f}–{max(aucs_all):.3f}")
    print(f"  ≥0.80 的比例: {100*sum(1 for a in aucs_all if a>=0.80)/len(aucs_all):.1f}%")
    print(f"  ≥0.85 的比例: {100*sum(1 for a in aucs_all if a>=0.85)/len(aucs_all):.1f}%")

ns_all = [int(r["max_n"]) for r in mod if r["max_n"] != ""]
if ns_all:
    ns_all.sort()
    print(f"\n报告的样本量: n={len(ns_all)}, 中位 {ns_all[len(ns_all)//2]}, "
          f"范围 {min(ns_all)}–{max(ns_all)}")
    print(f"  <100 例的比例: {100*sum(1 for x in ns_all if x<100)/len(ns_all):.1f}%")

yrs = {}
for r in mod:
    if r["year"]:
        yrs[r["year"]] = yrs.get(r["year"], 0) + 1
print("\n年份分布:", dict(sorted(yrs.items())))

# 保存汇总
with open(f"{OUT}/audit_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== ESCC 预后模型研究验证实践文献审计 =====\n\n")
    f.write(f"检索式: PubMed, 2015-2026, 英文\n")
    f.write(f"命中 {len(rows)} 篇; 预后建模研究且有摘要 {n} 篇\n\n")
    f.write(f"报告独立外部验证:             {ext:4d} ({100*ext/n:.1f}%)\n")
    f.write(f"报告内部拆分/交叉验证:        {intv:4d} ({100*intv/n:.1f}%)\n")
    f.write(f"两者均报告:                   {both:4d} ({100*both/n:.1f}%)\n")
    f.write(f"两者均未提及:                 {neither:4d} ({100*neither/n:.1f}%)\n")
    f.write(f"明确说明验证中重跑建模流程:   {rerun:4d} ({100*rerun/n:.1f}%)\n")
    f.write(f"报告 AUC/C-index:             {auc_rep:4d} ({100*auc_rep/n:.1f}%)\n")
    if aucs_all:
        f.write(f"\nAUC/C-index: 中位 {med:.3f}, 范围 {min(aucs_all):.3f}-{max(aucs_all):.3f}\n")
        f.write(f"  >=0.80: {100*sum(1 for a in aucs_all if a>=0.80)/len(aucs_all):.1f}%\n")
    if ns_all:
        f.write(f"\n样本量: 中位 {ns_all[len(ns_all)//2]}, 范围 {min(ns_all)}-{max(ns_all)}\n")
    f.write("\n年份分布: " + str(dict(sorted(yrs.items()))) + "\n")
print(f"\n汇总已保存 -> {OUT}/audit_summary.txt", flush=True)
