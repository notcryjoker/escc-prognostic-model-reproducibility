#!/usr/bin/env python3
"""
从 Ensembl BioMart 获取 Agilent Whole Human Genome 4x44K v2 (GPL18109) 探针 -> 基因符号映射。
用于 GSE53624 外部队列的探针注释。
"""
import urllib.parse, urllib.request, sys

out_path = sys.argv[1] if len(sys.argv) > 1 else "GPL18109_probe2gene_biomart.tsv"

xml = '''<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="agilent_wholegenome_4x44k_v2"/>
    <Attribute name="hgnc_symbol"/>
  </Dataset>
</Query>'''

url = "https://www.ensembl.org/biomart/martservice?query=" + urllib.parse.quote(xml)
print("正在查询 Ensembl BioMart ...", flush=True)
r = urllib.request.urlopen(url, timeout=300)
txt = r.read().decode("utf-8", errors="replace")
print(f"原始返回: {len(txt)} 字符", flush=True)

pairs = {}
for line in txt.strip().split("\n"):
    if not line.strip():
        continue
    parts = line.split("\t")
    if len(parts) < 2:
        continue
    probe, sym = parts[0].strip(), parts[1].strip()
    if not probe or not sym or sym == "":
        continue
    pairs.setdefault(probe, set()).add(sym)

n = 0
with open(out_path, "w") as fh:
    fh.write("probe\tgene\n")
    for probe, syms in pairs.items():
        for s in syms:
            fh.write(f"{probe}\t{s}\n")
            n += 1

print(f"唯一探针数: {len(pairs)}")
print(f"映射行数: {n}")
print(f"已保存: {out_path}", flush=True)
