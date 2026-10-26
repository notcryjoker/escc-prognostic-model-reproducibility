## ============================================================
## Figure 12 re-made: PRISMA flow + TRIPOD reporting-guideline audit
## ============================================================
suppressMessages({ library(ggplot2); library(patchwork); library(grid) })

d <- read.csv("results/literature_audit/tripod_audit.csv")
n <- nrow(d)
cat("Included:", n, "records\n")

## ---------- A: PRISMA flow (grid drawing) ----------
prisma_box <- function(x, y, w, h, lab, sub, fill, border) {
  grid.roundrect(x = x, y = y, width = w, height = h, r = unit(2, "mm"),
                 gp = gpar(fill = fill, col = border, lwd = 2))
  grid.text(lab, x = x, y = y + h * 0.16, gp = gpar(fontsize = 10.5, fontface = "bold"))
  grid.text(sub, x = x, y = y - h * 0.16, gp = gpar(fontsize = 9.5))
}
pA <- ggplot() + theme_void() +
  annotation_custom(grob = grobTree(
    rectGrob(gp = gpar(fill = "white", col = NA)),
    # three boxes
    roundrectGrob(x = 0.5, y = 0.86, width = 0.62, height = 0.14, r = unit(2, "mm"),
                  gp = gpar(fill = "#DCE9F7", col = "#2E6DA4", lwd = 2)),
    textGrob("PubMed search: 208 records", x = 0.5, y = 0.895,
             gp = gpar(fontsize = 11, fontface = "bold")),
    textGrob("2015–2026, English, ESCC + prognostic modelling + public data", x = 0.5, y = 0.845,
             gp = gpar(fontsize = 9)),
    roundrectGrob(x = 0.5, y = 0.60, width = 0.62, height = 0.14, r = unit(2, "mm"),
                  gp = gpar(fill = "#FFF3CD", col = "#B8860B", lwd = 2)),
    textGrob("3 records excluded", x = 0.5, y = 0.635, gp = gpar(fontsize = 11, fontface = "bold")),
    textGrob("not prognostic modelling studies (n = 3); no abstract (n = 0)", x = 0.5, y = 0.585,
             gp = gpar(fontsize = 9)),
    roundrectGrob(x = 0.5, y = 0.34, width = 0.62, height = 0.14, r = unit(2, "mm"),
                  gp = gpar(fill = "#D8F0DC", col = "#2E7D46", lwd = 2)),
    textGrob("205 records audited", x = 0.5, y = 0.375, gp = gpar(fontsize = 11, fontface = "bold")),
    textGrob("TRIPOD reporting assessment at the abstract level", x = 0.5, y = 0.325,
             gp = gpar(fontsize = 9)),
    # arrows
    segmentsGrob(x0 = 0.5, y0 = 0.79, x1 = 0.5, y1 = 0.675,
                 arrow = grid::arrow(length = unit(3, "mm"), type = "closed"),
                 gp = gpar(lwd = 1.6, col = "#555")),
    segmentsGrob(x0 = 0.5, y0 = 0.53, x1 = 0.5, y1 = 0.415,
                 arrow = grid::arrow(length = unit(3, "mm"), type = "closed"),
                 gp = gpar(lwd = 1.6, col = "#555"))
  )) +
  labs(title = "A  PRISMA flow") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5))

## ---------- B: TRIPOD item reporting rates ----------
it <- data.frame(
  item = c("Discrimination reported (AUC/C-index)", "Full model reported (coefficients/formula)",
           "Independent external validation reported", "Internal split/cross-validation reported", "Sample size reported",
           "Calibration reported", "Explicit re-run of modelling pipeline", "Multicentre/multi-source data reported"),
  pct = c(100 * sum(d$T3_discrimination) / n,
          100 * sum(d$T5_full_model) / n,
          100 * sum(d$T1_ext) / n,
          100 * sum(d$T1_int) / n,
          100 * sum(d$T6_sample_size) / n,
          100 * sum(d$T4_calibration) / n,
          100 * sum(d$T2_rerun) / n,
          100 * sum(d$T7_multisource) / n),
  cnt = c(sum(d$T3_discrimination), sum(d$T5_full_model), sum(d$T1_ext),
          sum(d$T1_int), sum(d$T6_sample_size), sum(d$T4_calibration),
          sum(d$T2_rerun), sum(d$T7_multisource))
)
it$item <- factor(it$item, levels = rev(it$item))
pB <- ggplot(it, aes(x = item, y = pct)) +
  geom_bar(stat = "identity", width = 0.7, fill = "#2E6DA4") +
  geom_text(aes(label = paste0(round(pct, 1), "% (", cnt, ")")),
            hjust = -0.1, size = 3.3) +
  coord_flip() + ylim(0, 60) +
  labs(x = NULL, y = "Reporting rate (%)",
       title = paste0("B  TRIPOD item reporting (n = ", n, ")")) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank())

## ---------- C: proportion meeting key item combinations ----------
comb <- data.frame(
  grp = c("Discrimination + calibration\n(TRIPOD 13b/13c)", "Model + sample size\n(TRIPOD 15/17)",
          "External validation + re-run pipeline", "Any validation"),
  pct = c(100 * sum(d$T3_discrimination & d$T4_calibration) / n,
          100 * sum(d$T5_full_model & d$T6_sample_size) / n,
          100 * sum(d$T1_ext & d$T2_rerun) / n,
          100 * sum(d$T1_ext | d$T1_int) / n)
)
comb$grp <- factor(comb$grp, levels = rev(comb$grp))
pC <- ggplot(comb, aes(x = grp, y = pct)) +
  geom_bar(stat = "identity", width = 0.62, fill = "#B2182B") +
  geom_text(aes(label = paste0(round(pct, 1), "%")), vjust = -0.5, size = 3.6) +
  ylim(0, 50) +
  labs(x = NULL, y = "Proportion meeting both (%)", title = "C  Key item combinations") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(size = 12, face = "bold"))

combined <- pA / (pB | pC) + plot_layout(heights = c(0.75, 1))
ggsave("figures/fig12_tripod_audit.png", combined, width = 12, height = 10.5, dpi = 300)
cat("Saved: figures/fig12_tripod_audit.png\n")
cat("Key combination: discrimination+calibration =", round(100 * sum(d$T3_discrimination & d$T4_calibration) / n, 1), "%\n")
