#!/usr/bin/env python3
"""
生成带图全文审读版 PDF（manuscript_EN_review_with_figures.pdf）。

用途：供作者/合作者通读、核对图文对应关系；不用于 JCE 投稿
（JCE 按 Elsevier 规范：正文纯文本 + 图单独上传 figures/submission/Figure*.png）。

做法：在每张图首次被正文引用的段落之后插入对应图（Figure 1 无正文引用，
置于 ## Methods 标题之后），再经 pandoc -> HTML -> weasyprint -> PDF。
"""
import os
import re
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "manuscript/manuscript_EN.md")
PDF_OUT = os.path.join(ROOT, "manuscript/manuscript_EN_review_with_figures.pdf")

FIG_DIR_REL = "../figures/submission"  # 相对 manuscript/ 目录


def build_md_with_figures():
    lines = open(MD).read().split("\n")
    # 每个图的首次引用行号
    first_cite = {}
    for i, ln in enumerate(lines):
        for n in range(1, 14):
            if n in first_cite:
                continue
            # 匹配 "Fig. N" 或 "Fig. NA"（A-E 面板字母），但不匹配 Fig. 10 中的 1 等
            if re.search(rf"Fig\.\s*{n}(?:[A-E])?(?![0-9A-Za-z])", ln):
                first_cite[n] = i
    # Figure 1 无正文引用 -> 置于 ## Methods 之后
    if 1 not in first_cite:
        for i, ln in enumerate(lines):
            if ln.strip() == "## Methods":
                first_cite[1] = i
                break

    # 按行号收集插入内容（同一行多处插入时按图号排序）
    inserts = {}
    for n, i in first_cite.items():
        inserts.setdefault(i, []).append(n)

    out = []
    for i, ln in enumerate(lines):
        out.append(ln)
        for n in sorted(inserts.get(i, [])):
            out.append("")
            out.append(f"![Figure {n}]({FIG_DIR_REL}/Figure{n}.png)")
    missing = [n for n in range(1, 14) if n not in first_cite]
    if missing:
        raise RuntimeError(f"以下图未找到插入位置: {missing}")
    return "\n".join(out)


CSS = """
@page { size: A4; margin: 2.0cm 2.0cm 2.2cm 2.0cm;
        @bottom-center { content: counter(page) " / " counter(pages);
                         font-size: 8.5pt; color: #666; } }
body { font-family: "Times New Roman", "Noto Serif CJK SC", serif;
       font-size: 10.5pt; line-height: 1.5; color: #111; }
h1 { font-family: "Noto Sans CJK SC", "Arial", sans-serif; font-size: 15pt;
     border-bottom: 2px solid #333; padding-bottom: 4px; margin: 0 0 12px 0; }
h2 { font-family: "Noto Sans CJK SC", "Arial", sans-serif; font-size: 12.5pt;
     margin: 16px 0 8px 0; color: #1a1a2e; border-bottom: 1px solid #ccc;
     padding-bottom: 2px; }
h3 { font-family: "Noto Sans CJK SC", "Arial", sans-serif; font-size: 11pt;
     margin: 12px 0 6px 0; }
p { margin: 6px 0; text-align: justify; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 9pt; }
th, td { border: 1px solid #999; padding: 4px 6px; text-align: left; }
th { background: #eef1f5; font-family: "Noto Sans CJK SC", "Arial", sans-serif; }
figure { margin: 14px 0; text-align: center; page-break-inside: avoid; }
figure img { max-width: 100%; border: 1px solid #d5d5d5; padding: 4px;
             background: #fff; }
figcaption { font-family: "Noto Sans CJK SC", "Arial", sans-serif;
             font-size: 8.5pt; color: #555; margin-top: 4px; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 9pt;
       background: #f3f3f3; padding: 0 3px; }
hr { border: none; border-top: 1px solid #bbb; margin: 14px 0; }
ul, ol { margin: 6px 0; padding-left: 22px; }
strong { color: #000; }
"""


def main():
    md_text = build_md_with_figures()
    with tempfile.NamedTemporaryFile("w", suffix=".md", delete=False,
                                     dir=os.path.join(ROOT, "manuscript")) as tmp:
        tmp.write(md_text)
        tmp_md = tmp.name
    tmp_html = tmp_md.replace(".md", ".html")
    try:
        subprocess.run(["pandoc", tmp_md, "-s", "--metadata", "pagetitle=Review"],
                       check=True, stdout=open(tmp_html, "w"))
        html = open(tmp_html).read()
        html = html.replace("</head>", f"<style>{CSS}</style></head>")
        open(tmp_html, "w").write(html)
        subprocess.run(["weasyprint", tmp_html, PDF_OUT], check=True,
                       capture_output=True)
        n_fig = md_text.count("![Figure")
        print(f"审读版 PDF 已生成: {PDF_OUT}")
        print(f"  嵌入图: {n_fig} 张 | 大小: "
              f"{os.path.getsize(PDF_OUT)/1024/1024:.1f} MB")
    finally:
        for p in (tmp_md, tmp_html):
            if os.path.exists(p):
                os.remove(p)


if __name__ == "__main__":
    main()
