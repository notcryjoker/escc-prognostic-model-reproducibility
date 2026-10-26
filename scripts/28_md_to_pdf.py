#!/usr/bin/env python3
"""
将投稿用 Markdown 文书转换为 PDF（pandoc md->html + weasyprint html->pdf）。
用法: python3 scripts/28_md_to_pdf.py <file1.md> [file2.md ...]
输出: 同名 .pdf，位于源文件同目录。
"""
import subprocess
import sys
import os
import tempfile

CSS = """
@page { size: A4; margin: 2.2cm 2.2cm 2.4cm 2.2cm;
        @bottom-center { content: counter(page) " / " counter(pages);
                         font-size: 8.5pt; color: #666; } }
body { font-family: "Times New Roman", "Noto Serif CJK SC", serif;
       font-size: 10.5pt; line-height: 1.5; color: #111; }
h1 { font-family: "Noto Sans CJK SC", "Arial", sans-serif;
     font-size: 15pt; border-bottom: 2px solid #333; padding-bottom: 4px;
     margin: 0 0 12px 0; }
h2 { font-family: "Noto Sans CJK SC", "Arial", sans-serif;
     font-size: 12pt; margin: 16px 0 8px 0; color: #222; }
h3 { font-family: "Noto Sans CJK SC", "Arial", sans-serif;
     font-size: 11pt; margin: 12px 0 6px 0; }
table { border-collapse: collapse; width: 100%; margin: 10px 0;
        font-size: 9pt; page-break-inside: auto; }
th, td { border: 1px solid #999; padding: 4px 6px; text-align: left;
         vertical-align: top; }
th { background: #eef1f5; font-family: "Noto Sans CJK SC", "Arial", sans-serif; }
tr { page-break-inside: avoid; }
blockquote { border-left: 3px solid #aaa; margin: 8px 0; padding: 2px 12px;
             color: #333; background: #fafafa; }
code { font-family: "DejaVu Sans Mono", monospace; font-size: 9pt;
       background: #f3f3f3; padding: 0 3px; }
hr { border: none; border-top: 1px solid #bbb; margin: 14px 0; }
ul, ol { margin: 6px 0 6px 0; padding-left: 22px; }
li { margin-bottom: 3px; }
strong { color: #000; }
"""


def convert(md_path: str) -> str:
    pdf_path = os.path.splitext(md_path)[0] + ".pdf"
    with tempfile.NamedTemporaryFile("w", suffix=".html", delete=False) as tmp:
        html_path = tmp.name
    try:
        subprocess.run(
            ["pandoc", md_path, "-s", "--metadata", "title=",
             "-c", "/dev/null"],
            check=True, stdout=open(html_path, "w"))
        # 注入内联 CSS（避免依赖外部文件）
        html = open(html_path).read()
        html = html.replace("</head>", f"<style>{CSS}</style></head>")
        open(html_path, "w").write(html)
        subprocess.run(["weasyprint", html_path, pdf_path], check=True,
                       capture_output=True)
        return pdf_path
    finally:
        if os.path.exists(html_path):
            os.remove(html_path)


if __name__ == "__main__":
    for md in sys.argv[1:]:
        pdf = convert(md)
        size = os.path.getsize(pdf)
        print(f"OK  {pdf}  ({size/1024:.0f} KB)")
