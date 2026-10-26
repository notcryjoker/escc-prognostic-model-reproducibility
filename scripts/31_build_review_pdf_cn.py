#!/usr/bin/env python3
"""
生成中文带图审读版 PDF（manuscript_CN_with_figures.pdf）。

在每张图首次被正文引用的段落之后插入对应投稿图（引用格式为"图 N"/"图 NX"）；
未被正文引用的图（如图 1）置于"## 2 材料与方法"标题之后。
经 pandoc -> HTML -> weasyprint -> PDF。
"""
import os
import re
import subprocess
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MD = os.path.join(ROOT, "manuscript/manuscript.md")
PDF_OUT = os.path.join(ROOT, "manuscript/manuscript_CN_with_figures.pdf")
FIG_DIR_REL = "../figures/submission"

CITE = {n: re.compile(rf"图\s*{n}(?![0-9])") for n in range(1, 14)}


def build_md_with_figures():
    lines = open(MD).read().split("\n")
    first_cite = {}
    for i, ln in enumerate(lines):
        for n, pat in CITE.items():
            if n not in first_cite and pat.search(ln):
                first_cite[n] = i
    # 未引用的图：置于 "## 2 材料与方法" 之后
    for n in range(1, 14):
        if n not in first_cite:
            for i, ln in enumerate(lines):
                if ln.strip().startswith("## 2 "):
                    first_cite[n] = i
                    break
    inserts = {}
    for n, i in first_cite.items():
        inserts.setdefault(i, []).append(n)

    out = []
    for i, ln in enumerate(lines):
        out.append(ln)
        for n in sorted(inserts.get(i, [])):
            out.append("")
            out.append(f"![图 {n}]({FIG_DIR_REL}/Figure{n}.png)")
    missing = [n for n in range(1, 14) if n not in first_cite]
    if missing:
        raise RuntimeError(f"以下图未找到插入位置: {missing}")
    return "\n".join(out)


CSS = """
@page { size: A4; margin: 2.0cm 2.0cm 2.2cm 2.0cm;
        @bottom-center { content: counter(page) " / " counter(pages);
                         font-size: 8.5pt; color: #666; } }
body { font-family: "Noto Serif CJK SC", "Times New Roman", serif;
       font-size: 10.5pt; line-height: 1.65; color: #111; }
h1 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 16pt;
     border-bottom: 2px solid #333; padding-bottom: 4px; margin: 0 0 12px 0;
     line-height: 1.4; }
h2 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 13pt;
     margin: 18px 0 8px 0; color: #1a1a2e; border-bottom: 1px solid #ccc;
     padding-bottom: 2px; }
h3 { font-family: "Noto Sans CJK SC", sans-serif; font-size: 11.5pt;
     margin: 12px 0 6px 0; }
p { margin: 6px 0; text-align: justify; }
table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 9pt; }
th, td { border: 1px solid #999; padding: 4px 6px; text-align: left; }
th { background: #eef1f5; font-family: "Noto Sans CJK SC", sans-serif; }
figure { margin: 14px 0; text-align: center; page-break-inside: avoid; }
figure img { max-width: 100%; border: 1px solid #d5d5d5; padding: 4px;
             background: #fff; }
figcaption { font-family: "Noto Sans CJK SC", sans-serif;
             font-size: 8.5pt; color: #555; margin-top: 4px; }
blockquote { border-left: 3px solid #aaa; margin: 8px 0; padding: 2px 12px;
             color: #333; background: #fafafa; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 9pt;
       background: #f3f3f3; padding: 0 3px; }
hr { border: none; border-top: 1px solid #bbb; margin: 14px 0; }
ul, ol { margin: 6px 0; padding-left: 22px; }
li { margin-bottom: 3px; }
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
        subprocess.run(["pandoc", tmp_md, "-s", "--metadata", "pagetitle=中文审读版"],
                       check=True, stdout=open(tmp_html, "w"))
        html = open(tmp_html).read()
        html = html.replace("</head>", f"<style>{CSS}</style></head>")
        open(tmp_html, "w").write(html)
        subprocess.run(["weasyprint", tmp_html, PDF_OUT], check=True,
                       capture_output=True)
        n_fig = md_text.count("![图")
        print(f"中文审读版 PDF 已生成: {PDF_OUT}")
        print(f"  嵌入图: {n_fig} 张 | 大小: "
              f"{os.path.getsize(PDF_OUT)/1024/1024:.1f} MB")
    finally:
        for p in (tmp_md, tmp_html):
            if os.path.exists(p):
                os.remove(p)


if __name__ == "__main__":
    main()
