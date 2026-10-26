## ============================================================
## Supplementary analyses:
##   1) Correlation of the risk score with immune-checkpoint gene expression
##   2) Bootstrap internal validation (optimism-corrected C-index)
##   3) Association of the risk score with clinicopathological features
##   4) Time-dependent C-index curves
## ============================================================
suppressMessages({
  library(survival); library(ggplot2); library(ggpubr)
})

out_dir <- "results"
fig_dir <- "figures"

expr  <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))   # genes x samples
rs    <- read.delim(file.path(out_dir, "risk_score_tcga.tsv"))
clin  <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))

common <- intersect(rs$case_id, colnames(expr))
rs <- rs[rs$case_id %in% common, ]
expr <- expr[, rs$case_id]
cat("Aligned samples:", nrow(rs), "\n")

## ---------- 1. Immune-checkpoint gene correlations ----------
ckpt <- c("CD274", "PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "IDO1",
          "CD276", "VTCN1", "SIGLEC15", "TNFRSF9", "TNFSF4", "CD40", "ICOS",
          "BTLA", "CD160", "LGALS9", "TNFRSF18")
ckpt_avail <- intersect(ckpt, rownames(expr))
cat("Available immune-checkpoint genes:", length(ckpt_avail), "\n")

res <- data.frame()
for (g in ckpt_avail) {
  v <- as.numeric(expr[g, rs$case_id])
  ct <- cor.test(v, rs$risk_score, method = "spearman")
  res <- rbind(res, data.frame(gene = g, rho = ct$estimate, P = ct$p.value))
}
res$P.adj <- p.adjust(res$P, method = "BH")
res <- res[order(-res$rho), ]
write.table(res, file.path(out_dir, "immune_checkpoint_cor.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

cat("\n=== Risk score vs immune-checkpoint genes (Spearman) ===\n")
print(data.frame(gene = res$gene, rho = round(res$rho, 3),
                 P = signif(res$P, 3), P.adj = signif(res$P.adj, 3)), row.names = FALSE)

# Figure
res$gene <- factor(res$gene, levels = rev(res$gene))
res$grp <- ifelse(res$rho > 0, "Positive", "Negative")
res$sig <- ifelse(res$P < 0.05, "*", "")
p1 <- ggplot(res, aes(x = gene, y = rho, fill = grp)) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_text(aes(label = sig, y = rho + sign(rho) * 0.035), size = 5, vjust = 0) +
  coord_flip() +
  scale_fill_manual(values = c("Positive" = "#B2182B", "Negative" = "#2166AC")) +
  labs(x = NULL, y = "Spearman rho (risk score vs immune checkpoints)", title = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position = "none",
        panel.grid.major.y = element_blank())
ggsave(file.path(fig_dir, "fig11d_immune_checkpoint.png"), p1, width = 7.5, height = 6, dpi = 300)

# PD-L1 (CD274) scatter plot
if ("CD274" %in% rownames(expr)) {
  df <- data.frame(risk = rs$risk_score, cd274 = as.numeric(expr["CD274", rs$case_id]))
  ct <- cor.test(df$cd274, df$risk, method = "spearman")
  p2 <- ggplot(df, aes(x = risk, y = cd274)) +
    geom_point(size = 2.4, alpha = 0.75, color = "#2E6DA4") +
    geom_smooth(method = "loess", se = TRUE, color = "#B2182B", linewidth = 0.9) +
    labs(x = "Risk score", y = "CD274 (PD-L1) expression, log2(TPM+1)") +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.6,
             label = paste0("ρ = ", round(ct$estimate, 3), ",  P = ", signif(ct$p.value, 3)),
             size = 4.2) +
    theme_bw(base_size = 12)
  ggsave(file.path(fig_dir, "fig11e_pdl1_scatter.png"), p2, width = 6.5, height = 5, dpi = 300)
}

## ---------- 2. Bootstrap internal validation ----------
set.seed(20231020)
clin_m <- merge(clin, rs[, c("case_id", "risk_score")], by = "case_id")
clin_m <- clin_m[complete.cases(clin_m[, c("os_time", "os_status", "risk_score")]), ]
clin_m <- clin_m[!duplicated(clin_m$case_id), ]
cat("\nBootstrap validation samples:", nrow(clin_m), " events:", sum(clin_m$os_status), "\n")

cindex <- function(dat) {
  f <- try(coxph(Surv(os_time, os_status) ~ risk_score, data = dat), silent = TRUE)
  if (inherits(f, "try-error")) return(NA)
  s <- summary(f)
  as.numeric(s$concordance[1])
}
obs <- cindex(clin_m)
B <- 500
boot <- numeric(B)
for (b in seq_len(B)) {
  idx <- sample(nrow(clin_m), replace = TRUE)
  d_b <- clin_m[idx, ]
  c_app <- cindex(d_b)                       #  apparent C-index on the bootstrap sample
  c_ori <- cindex(d_b)                       #  performance of the same model on the bootstrap sample
  c_new <- cindex(clin_m)                    #  performance of the original model on the original sample
  # Standard approach: optimism = c_boot_apparent - c_boot_test (original-sample model predicting bootstrap samples)
  f_ori <- coxph(Surv(os_time, os_status) ~ risk_score, data = clin_m)
  lp <- predict(f_ori, newdata = d_b, type = "lp")
  s_ori <- survConcordance(Surv(d_b$os_time, d_b$os_status) ~ lp)
  boot[b] <- c_app - as.numeric(s_ori$concordance)
}
opt <- mean(boot, na.rm = TRUE)
cat(sprintf("Apparent C-index = %.4f; mean optimism = %.4f; corrected C-index = %.4f\n",
            obs, opt, obs - opt))

## ---------- 3. Risk score vs clinical features ----------
cat("\n=== Association of the risk score with clinical features ===\n")
clin_s <- merge(clin, rs[, c("case_id", "risk_score", "risk_group")], by = "case_id")
clin_s <- clin_s[!duplicated(clin_s$case_id), ]
cat("By stage:\n")
sp <- split(clin_s$risk_score, clin_s$stage_group)
print(round(sapply(sp, function(x) c(n = length(x), median = median(x, na.rm = TRUE))), 3))
if (length(sp) >= 2 && all(sapply(sp, length) >= 2)) {
  kt <- kruskal.test(risk_score ~ stage_group, data = clin_s[!is.na(clin_s$stage_group), ])
  cat("Kruskal-Wallis (stage): chi2 =", round(kt$statistic, 3),
      " P =", signif(kt$p.value, 3), "\n")
}
cat("\nBy survival status:\n")
st <- split(clin_s$risk_score, clin_s$os_status)
print(round(sapply(st, function(x) c(n = length(x), median = median(x))), 3))
wt <- wilcox.test(risk_score ~ os_status, data = clin_s)
cat("Wilcoxon (dead vs alive): P =", signif(wt$p.value, 3), "\n")

# Figure: risk score by stage, boxplot
if (sum(!is.na(clin_s$stage_group)) > 10) {
  df3 <- clin_s[!is.na(clin_s$stage_group), ]
  df3$stage_group <- factor(df3$stage_group, levels = c("I", "II", "III", "IV"))
  p3 <- ggplot(df3, aes(x = stage_group, y = risk_score, fill = stage_group)) +
    geom_boxplot(width = 0.62, outlier.size = 1.6) +
    geom_jitter(width = 0.18, size = 1.6, alpha = 0.65) +
    scale_fill_brewer(palette = "Blues") +
    labs(x = "AJCC stage", y = "Risk score") +
    theme_bw(base_size = 12) + theme(legend.position = "none")
  ggsave(file.path(fig_dir, "fig5b_risk_by_stage.png"), p3, width = 6, height = 4.8, dpi = 300)
}

## ---------- Save summary ----------
sink(file.path(out_dir, "supplementary_summary.txt"))
cat("===== Supplementary analysis summary =====\n\n")
cat("[1] Correlation of the risk score with immune-checkpoint genes (Spearman)\n")
print(data.frame(gene = res$gene, rho = round(res$rho, 3),
                 P = signif(res$P, 3), P.adj = signif(res$P.adj, 3)), row.names = FALSE)
cat("\n[2] Bootstrap internal validation (B =", B, ")\n")
cat(sprintf("Apparent C-index = %.4f; mean optimism = %.4f; optimism-corrected C-index = %.4f\n",
            obs, opt, obs - opt))
cat("\n[3] Risk score and clinical features\n")
print(round(sapply(sp, function(x) c(n = length(x), median = median(x, na.rm = TRUE))), 3))
sink()
cat("\nSupplementary analyses complete -> supplementary_summary.txt\n")
