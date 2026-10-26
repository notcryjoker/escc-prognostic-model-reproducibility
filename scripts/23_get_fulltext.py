#!/usr/bin/env python3
"""
为冲 JCE 补强：获取全文以执行 PROBAST 偏倚风险评估
  1. 从 PubMed 获取 205 篇文献的 DOI
  2. 通过 Unpaywall 查询合法开放获取副本链接
  3. 下载 HTML 全文并提取正文文本
输出：fulltext/ 目录 + 可获取率报告
"""
import urllib.request, urllib.parse, json, csv, re, time, os, sys, subprocess

BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "literature_audit")
FT = f"{OUT}/fulltext"
os.makedirs(FT, exist_ok=True)
EMAIL = "openaccess@example.com"


def fetch(url, tries=3, timeout=90):
    for i in range(tries):
        try:
            req = urllib.request.Request(url, headers={
                "User-Agent": "Mozilla/5.0 (compatible; research-audit)"})
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return r.read()
        except Exception as e:
            if i == tries - 1:
                return None
            time.sleep(4)
    return None


recs = list(csv.DictReader(open(f"{OUT}/audit_records.csv", encoding="utf-8")))
recs = [r for r in recs if r["is_model_study"] == "1" and r["has_abstract"] == "1"]
pmids = [r["pmid"] for r in recs]
print(f"纳入 {len(recs)} 篇", flush=True)

# ---------- 1. 获取 DOI ----------
print("\n[1] 获取 DOI ...", flush=True)
doi_map = {}
for i in range(0, len(pmids), 100):
    chunk = pmids[i:i + 100]
    xml = fetch(f"{BASE}/efetch.fcgi?db=pubmed&id={','.join(chunk)}&retmode=xml")
    if not xml:
        continue
    xml = xml.decode("utf-8", errors="replace")
    arts = re.findall(r"<PubmedArticle>(.*?)</PubmedArticle>", xml, re.S)
    for a in arts:
        pid = re.search(r"<PMID[^>]*>(\d+)</PMID>", a)
        if not pid:
            continue
        dois = re.findall(r'<ArticleId[^>]*IdType="doi"[^>]*>([^<]+)</ArticleId>', a)
        if dois:
            doi_map[pid.group(1)] = dois[0].strip()
    print(f"   已获取 DOI {len(doi_map)}/{len(pmids)}", flush=True)
    time.sleep(3)
print(f"有 DOI: {len(doi_map)} 篇", flush=True)

# ---------- 2. Unpaywall 查询 OA 链接 ----------
print("\n[2] 查询 Unpaywall 开放获取链接 ...", flush=True)
oa = {}
no_doi = 0
for k, pmid in enumerate(pmids, 1):
    doi = doi_map.get(pmid)
    if not doi:
        no_doi += 1
        continue
    url = f"https://api.unpaywall.org/v2/{urllib.parse.quote(doi)}?email={EMAIL}"
    raw = fetch(url, tries=2, timeout=45)
    if not raw:
        continue
    try:
        j = json.loads(raw.decode("utf-8", errors="replace"))
    except Exception:
        continue
    if j.get("is_oa"):
        best = j.get("best_oa_location") or {}
        link = best.get("url_for_pdf") or best.get("url_for_landing_page") or best.get("url")
        if link:
            oa[pmid] = {"doi": doi, "url": link,
                        "host": best.get("host_type", ""),
                        "evidence": best.get("evidence", "")}
    if k % 40 == 0:
        print(f"   已查询 {k}/{len(pmids)}, 找到 OA {len(oa)}", flush=True)
    time.sleep(1.1)          # Unpaywall 礼貌限速
print(f"\n开放获取可获取: {len(oa)}/{len(pmids)} 篇 ({100*len(oa)/len(pmids):.1f}%)", flush=True)
hosts = {}
for v in oa.values():
    hosts[v["host"]] = hosts.get(v["host"], 0) + 1
print("来源分布:", hosts, flush=True)

with open(f"{OUT}/oa_links.csv", "w", newline="", encoding="utf-8") as f:
    w = csv.writer(f)
    w.writerow(["pmid", "doi", "url", "host_type", "evidence"])
    for pmid, v in oa.items():
        w.writerow([pmid, v["doi"], v["url"], v["host"], v["evidence"]])

# ---------- 3. 下载全文 ----------
print("\n[3] 下载全文 ...", flush=True)
def html2text(h):
    h = re.sub(r"<(script|style)[^>]*>.*?</\1>", " ", h, flags=re.S | re.I)
    h = re.sub(r"<[^>]+>", " ", h)
    return re.sub(r"\s+", " ", h).strip()


def pdf2text(p):
    try:
        r = subprocess.run(["pdftotext", "-", "-"], input=p,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=90)
        return r.stdout.decode("utf-8", errors="replace")
    except Exception:
        return ""


ok = 0
for k, (pmid, v) in enumerate(oa.items(), 1):
    dest = f"{FT}/{pmid}.txt"
    if os.path.exists(dest) and os.path.getsize(dest) > 2000:
        ok += 1
        continue
    raw = fetch(v["url"], tries=2, timeout=90)
    if not raw:
        continue
    if raw[:4] == b"%PDF":
        txt = pdf2text(raw)
    else:
        txt = html2text(raw.decode("utf-8", errors="replace"))
    if len(txt) > 3000:
        open(dest, "w", encoding="utf-8").write(txt)
        ok += 1
    if k % 25 == 0:
        print(f"   已下载 {ok}/{len(oa)}", flush=True)
    time.sleep(0.8)

print(f"\n全文下载成功: {ok} 篇 ({100*ok/len(pmids):.1f}% of 205)", flush=True)

with open(f"{OUT}/fulltext_summary.txt", "w", encoding="utf-8") as f:
    f.write("===== 全文获取汇总 =====\n\n")
    f.write(f"纳入 {len(pmids)} 篇; 有 DOI {len(doi_map)} 篇; "
            f"Unpaywall 报告 OA {len(oa)} 篇; 全文下载成功 {ok} 篇\n")
    f.write(f"来源分布: {hosts}\n")
print("汇总已保存 -> fulltext_summary.txt", flush=True)
