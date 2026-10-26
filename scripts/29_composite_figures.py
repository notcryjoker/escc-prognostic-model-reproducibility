#!/usr/bin/env python3
"""
合成投稿用图文件：一张图一个文件（Elsevier Editorial Manager 要求），
多面板图加 A/B/C... 面板字母（与正文 Fig. XA / Fig. 7D 等引用对应）。

输出: figures/submission/Figure1.png ... Figure13.png（共 13 个，即上传文件）

布局说明:
  Figure 3  在既有双面板图上叠加 A/B 字母（正文引用 Fig. 3A, B）
  Figure 5  A=风险热图  B=风险评分 vs AJCC 分期（横向双联）
  Figure 7  A=ROC  B=全流程重跑拆分  C=固定模型拆分  D=学习曲线  E=嵌套CV（3列x2行）
  Figure 9  A=GO  B=KEGG ORA（横向双联）
  Figure 10 A=GSEA KEGG  B=GSEA Hallmark（横向双联）
  Figure 11 A=ssGSEA箱线  B=免疫细胞-风险相关  C=浸润热图  D=检查点相关  E=PD-L1散点（2+3双行）
  其余图（1,2,4,6,8,12,13）为单文件或已含面板字母，原样复制
"""
import os
import shutil

import numpy as np
from PIL import Image, ImageDraw, ImageFont

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")
OUT = os.path.join(FIG, "submission")
os.makedirs(OUT, exist_ok=True)

FONT_PATH = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf"
BG = (255, 255, 255)
MARGIN = 60          # 画布外边距
GAP = 90             # 面板间距
LETTER_H = 150       # 每面板顶部预留的字母区高度
MAX_SIDE = 4800      # 合成图最大边长（超出则等比缩小）


def load(name):
    return Image.open(os.path.join(FIG, name)).convert("RGB")


def fit(img, box_w, box_h):
    """等比缩放放入 box，返回缩放后的图（不放大超过原图 1.2 倍）。"""
    w, h = img.size
    s = min(box_w / w, box_h / h, 1.2)
    if s < 1:
        img = img.resize((max(1, int(w * s)), max(1, int(h * s))), Image.LANCZOS)
    return img


def letter_font(px=110):
    return ImageFont.truetype(FONT_PATH, px)


def paste_panel(canvas, img, letter, cx, cy_top):
    """在画布上放置一个面板：字母区 + 居中面板。cx 为面板区域中心 x。"""
    dr = ImageDraw.Draw(canvas)
    f = letter_font(110)
    # 字母绘制在面板区域左上角
    bbox = dr.textbbox((0, 0), letter, font=f)
    lw, lh = bbox[2] - bbox[0], bbox[3] - bbox[1]
    dr.text((cx - img.width // 2, cy_top + (LETTER_H - lh) // 2 - bbox[1]),
            letter, font=f, fill=(0, 0, 0))
    canvas.paste(img, (cx - img.width // 2, cy_top + LETTER_H))


def compose_rows(rows, cell_h, fname):
    """rows: [[(letter, img), ...], ...] 每行面板数可不同；各行在画布上水平居中。"""
    # 计算每行的实际宽度
    row_imgs = []
    for row in rows:
        fitted = [(letter, fit(img, 10 ** 9, cell_h - LETTER_H)) for letter, img in row]
        row_imgs.append(fitted)
    row_widths = [sum(im.width for _, im in r) + GAP * (len(r) - 1) for r in row_imgs]
    W = max(row_widths) + 2 * MARGIN
    H = sum(cell_h for _ in rows) + GAP * (len(rows) - 1) + 2 * MARGIN

    canvas = Image.new("RGB", (W, H), BG)
    y = MARGIN
    for r, rw in zip(row_imgs, row_widths):
        x = MARGIN + (W - 2 * MARGIN - rw) // 2
        for letter, im in r:
            paste_panel(canvas, im, letter, x + im.width // 2, y)
            x += im.width + GAP
        y += cell_h + GAP
    save(canvas, fname)


def compose_grid(panels, ncols, nrows, fname):
    """均匀网格：panels 按行优先填入，(letter, img) 列表。"""
    cell_h = None  # 由下方统一计算
    fitted = []
    # 先按统一 cell 高度拟合
    return fitted  # placeholder（未使用）


def save(canvas, fname):
    w, h = canvas.size
    if max(w, h) > MAX_SIDE:
        s = MAX_SIDE / max(w, h)
        canvas = canvas.resize((int(w * s), int(h * s)), Image.LANCZOS)
    path = os.path.join(OUT, fname)
    canvas.save(path, optimize=True)
    print(f"  {fname}: {canvas.size[0]}x{canvas.size[1]}  "
          f"{os.path.getsize(path)/1024/1024:.2f} MB")


def overlay_letters_3(fname):
    """Figure 3：在既有双面板图上按白色分隔带定位，叠加 A / B。"""
    img = load("fig3_lasso.png")
    arr = np.asarray(img.convert("L"))
    white_frac = (arr > 245).mean(axis=1)
    bands, start = [], None
    for i, v in enumerate(white_frac):
        if v > 0.985:
            if start is None:
                start = i
        elif start is not None:
            if i - start > 40:
                bands.append((start, i))
            start = None
    if start is not None and len(white_frac) - start > 40:
        bands.append((start, len(white_frac)))
    if not bands:
        raise RuntimeError("未找到 Figure 3 的白色分隔带")
    top, bottom = max(bands, key=lambda b: b[1] - b[0])
    print(f"  Figure 3 分隔带: 行 {top}-{bottom}")
    dr = ImageDraw.Draw(img)
    f = letter_font(130)
    dr.text((MARGIN, top - 190), "A", font=f, fill=(0, 0, 0))
    dr.text((MARGIN, bottom + 40), "B", font=f, fill=(0, 0, 0))
    save(img, fname)


print("合成投稿用图 ->", OUT)

# ---------- Figure 3：叠加字母 ----------
overlay_letters_3("Figure3.png")

# ---------- Figure 5：A 热图 + B 分期（横向） ----------
compose_rows([[("A", load("fig5_risk_heatmap.png")),
               ("B", load("fig5b_risk_by_stage.png"))]], 2100, "Figure5.png")

# ---------- Figure 7：3 列 x 2 行（末格留空） ----------
compose_rows([[("A", load("fig7_roc.png")),
               ("B", load("fig7b_split_validation.png")),
               ("C", load("fig7c_split_fixedmodel.png"))],
              [("D", load("fig7d_learning_curve.png")),
               ("E", load("fig7e_nested_cv_compare.png"))]], 1800, "Figure7.png")

# ---------- Figure 9 / 10：横向双联 ----------
compose_rows([[("A", load("fig9_go.png")),
               ("B", load("fig9b_kegg.png"))]], 2200, "Figure9.png")
compose_rows([[("A", load("fig10_gsea_kegg.png")),
               ("B", load("fig10b_gsea_hallmark.png"))]], 2200, "Figure10.png")

# ---------- Figure 11：2+3 双行 ----------
compose_rows([[("A", load("fig11_immune_box.png")),
               ("B", load("fig11b_immune_cor.png"))],
              [("C", load("fig11c_immune_heatmap.png")),
               ("D", load("fig11d_immune_checkpoint.png")),
               ("E", load("fig11e_pdl1_scatter.png"))]], 2000, "Figure11.png")

# ---------- 单文件图：原样复制 ----------
for n in [1, 2, 4, 6, 8, 12, 13]:
    src = load(f"fig{n}_" + {
        1: "flowchart", 2: "deg_panel", 4: "forest", 6: "km",
        8: "nomogram", 12: "tripod_audit", 13: "probst"}[n] + ".png")
    save(src, f"Figure{n}.png")

print("\n完成。投稿上传用图共", len(os.listdir(OUT)), "个：")
for f in sorted(os.listdir(OUT)):
    print(" ", f)
