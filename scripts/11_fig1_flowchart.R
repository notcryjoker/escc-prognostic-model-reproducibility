## ============================================================
## Fig 1: Study workflow (grid drawing) - main trunk + right-hand branches
## ============================================================
suppressMessages(library(grid))

png("figures/fig1_flowchart.png", width = 2200, height = 2600, res = 200)
pushViewport(viewport(w = 1, h = 1, x = 0.5, y = 0.5))

# ---- Palette ----
c_data  <- "#DCE9F7"; b_data  <- "#2E6DA4"
c_proc  <- "#FFF3CD"; b_proc  <- "#B8860B"
c_deg   <- "#D8F0DC"; b_deg   <- "#2E7D46"
c_mod   <- "#F8D7DA"; b_mod   <- "#B02A37"
c_val   <- "#E2D9F3"; b_val   <- "#6F42C1"
c_func  <- "#D1ECF1"; b_func  <- "#0C7C94"

box <- function(xc, yc, w, h, label, fill, border, fs = 12, fw = "bold") {
  grid.roundrect(x = xc, y = yc, width = w, height = h, r = unit(3, "mm"),
                 gp = gpar(fill = fill, col = border, lwd = 2.2))
  grid.text(label, x = xc, y = yc, gp = gpar(fontsize = fs, fontface = fw))
}
arrow <- function(x1, y1, x2, y2, lwd = 1.8) {
  grid.segments(x1, y1, x2, y2,
                arrow = grid::arrow(length = unit(3, "mm"), type = "closed"),
                gp = gpar(col = "#4A4A4A", lwd = lwd))
}
elbow <- function(x1, y1, x2, y2, lwd = 1.8) {  # horizontal first, then vertical (from the right edge of the left box to above the target)
  grid.segments(x1, y1, x2, y1, gp = gpar(col = "#4A4A4A", lwd = lwd))
  arrow(x2, y1, x2, y2, lwd = lwd)
}

XM <- 0.36   # main trunk centre
WM <- 0.52   # main trunk width
XR <- 0.80   # right-branch centre
WR <- 0.30   # right-branch width

# ================= Title =================
grid.text("Study workflow", x = 0.5, y = 0.968,
          gp = gpar(fontsize = 18, fontface = "bold", col = "#1A1A1A"))

# ================= L1: Data acquisition =================
grid.text("① Public data acquisition", x = 0.5, y = 0.938, gp = gpar(fontsize = 12, col = "#555555", fontface = "bold"))
bw <- 0.215; gp_ <- 0.017
xs <- c(0.5 - 1.5 * (bw + gp_), 0.5 - 0.5 * (bw + gp_),
        0.5 + 0.5 * (bw + gp_), 0.5 + 1.5 * (bw + gp_))
labs <- c("TCGA-ESCC\n95 tumour samples\n(discovery cohort)", "GSE20347\n17 paired samples\n(GPL571)",
          "GSE45670\n28 tumour / 10 normal\n(GPL570)", "GSE53624\n119 paired samples\n(GPL18109)")
y1 <- 0.882
for (i in 1:4) box(xs[i], y1, bw, 0.082, labs[i], c_data, b_data, fs = 10.5)
arrow(0.5, y1 - 0.041, 0.5, y1 - 0.062)

# ================= L2: Preprocessing =================
y2 <- 0.808
box(XM, y2, WM, 0.058, "② Data preprocessing and annotation harmonisation\nSTAR counts → TPM; probe/ENSG → gene symbols;\nlog2 transformation; sample pairing",
    c_proc, b_proc, fs = 10.5)
arrow(0.5, y2 - 0.029, 0.5, y2 - 0.052)

# ================= L3: DEG =================
y3 <- 0.738
box(XM, y3, WM, 0.056, "③ Differential expression analysis (limma)\nthresholds: |log2FC| ≥ 1, FDR < 0.05", c_deg, b_deg, fs = 10.5)
arrow(0.5, y3 - 0.028, 0.5, y3 - 0.050)

# ---- Consensus note + consensus box ----
y3b <- 0.682
box(XM, y3b, WM, 0.040, "Three-source consensus: significant in ≥2 data sources\nwith consistent direction (including ≥1 GEO cohort)",
    c_deg, b_deg, fs = 9.5, fw = "plain")
arrow(0.5, y3b - 0.020, 0.5, y3b - 0.040)
y3c <- 0.632
box(XM, y3c, 0.40, 0.040, "Consensus DEGs: 719 genes\n(406 up- / 313 down-regulated)", "#FFFFFF", b_deg, fs = 10)
arrow(0.5, y3c - 0.020, 0.5, y3c - 0.042)

# ---- Right branch: functional enrichment (from consensus DEGs) ----
yR1 <- 0.632
box(XR, yR1, WR, 0.070, "Functional enrichment analysis\n\n• GO (BP/MF/CC)\n• KEGG ORA + GSEA\n• Hallmark GSEA",
    c_func, b_func, fs = 9.5, fw = "plain")
elbow(XM + 0.20, y3c, XR, yR1 + 0.035)

# ================= L4: Screening + LASSO =================
y4 <- 0.556
box(XM, y4, WM, 0.046, "④ Univariable Cox regression screening (TCGA, P < 0.05)\n→ 48 candidate genes",
    "#FFFFFF", b_mod, fs = 10)
arrow(0.5, y4 - 0.023, 0.5, y4 - 0.044)
y4b <- 0.494
box(XM, y4b, WM, 0.046, "⑤ LASSO-Cox penalised regression (10-fold cross-validation)\n→ 17-gene risk-score model",
    "#FFFFFF", b_mod, fs = 10)

# ---- Right branch: immune infiltration (from the LASSO model) ----
yR2 <- 0.494
box(XR, yR2, WR, 0.070, "Tumour immune infiltration analysis\n\n• ssGSEA (27 immune cell types)\n• high- vs low-risk group comparison\n• correlation with risk score",
    c_func, b_func, fs = 9.5, fw = "plain")
elbow(XM + 0.26, y4b, XR, yR2 + 0.035)

# ================= L5: Validation =================
arrow(0.5, y4b - 0.023, 0.5, y4b - 0.046)
y5 <- 0.416
box(0.27, y5, 0.40, 0.084, "⑥ Internal validation (TCGA-ESCC, n = 94)\n\n• Kaplan–Meier analysis\n• timeROC (1/3/5 years)\n• multivariable Cox (age, stage)\n• nomogram",
    c_val, b_val, fs = 9.5, fw = "plain")
box(0.73, y5, 0.40, 0.084, "⑦ External validation (GSE53624)\n\n• apply TCGA model coefficients\n• grouping at the median cut-off\n• Kaplan–Meier + timeROC",
    c_val, b_val, fs = 9.5, fw = "plain")
arrow(0.5, y4b - 0.023, 0.5, y5 + 0.042)

# ================= L6: Output =================
arrow(0.5, y5 - 0.042, 0.5, y5 - 0.066)
y6 <- 0.330
box(0.5, y6, 0.70, 0.052, "⑧ Performance evaluation and biological interpretation of the prognostic model\n(KM curves / ROC / nomogram / enrichment maps / immune maps)",
    "#EDEDED", "#444444", fs = 10.5)

# ================= Legend =================
ly <- 0.245
lx <- c(0.13, 0.37, 0.60, 0.82)
lt <- c("Data acquisition", "Preprocessing", "Differential expression", "Model construction")
lf <- c(c_data, c_proc, c_deg, c_mod); lb <- c(b_data, b_proc, b_deg, b_mod)
for (i in 1:4) {
  grid.roundrect(x = lx[i], y = ly, width = 0.028, height = 0.024, r = unit(1.5, "mm"),
                 gp = gpar(fill = lf[i], col = lb[i], lwd = 1.8))
  grid.text(lt[i], x = lx[i] + 0.048, y = ly, gp = gpar(fontsize = 10, just = "left"))
}
ly <- 0.205
lt2 <- c("Validation analysis", "Functional annotation"); lf2 <- c(c_val, c_func); lb2 <- c(b_val, b_func)
for (i in 1:2) {
  grid.roundrect(x = lx[i], y = ly, width = 0.028, height = 0.024, r = unit(1.5, "mm"),
                 gp = gpar(fill = lf2[i], col = lb2[i], lwd = 1.8))
  grid.text(lt2[i], x = lx[i] + 0.048, y = ly, gp = gpar(fontsize = 10, just = "left"))
}

# ================= Footnote =================
grid.text("ESCC: oesophageal squamous cell carcinoma; TCGA: The Cancer Genome Atlas; GEO: Gene Expression Omnibus; DEG: differentially expressed gene;\nORA: over-representation analysis; GSEA: gene set enrichment analysis; ssGSEA: single-sample gene set enrichment analysis.",
          x = 0.5, y = 0.145, gp = gpar(fontsize = 9, col = "#666666"))

dev.off()
cat("fig1_flowchart.png generated\n")
