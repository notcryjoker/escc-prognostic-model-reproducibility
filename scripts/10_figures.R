# ============================================================================
# 10_figures.R
# Purpose: assemble publication-grade figures (300 dpi, unified styling)
#   Fig1 study workflow + cohort overview
#   Fig2 DEG volcano plots (GSE20347/GSE45670/TCGA) + consensus-gene Venn
#   Fig3 LASSO coefficient path + CV curve
#   Fig4 model-gene forest plot (TCGA univariable)
#   Fig5 risk-score heatmap + survival status
#   Fig6 KM/ROC (TCGA + GSE53624 produced separately by 06/09)
#   Fig7 enrichment plots (GO dot plot + KEGG bar plot + GSEA)
#   Fig8 immune infiltration (heatmap + boxplot + correlation)
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(ggplot2); library(data.table); library(glmnet)
  library(VennDiagram); library(grid); library(patchwork)
})

base <- "/workspace/escc_bioinformatics"
out_dir <- file.path(base, "results")
fig_dir <- file.path(base, "figures")
dir.create(fig_dir, showWarnings = FALSE)

theme_pub <- theme_bw(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        axis.title = element_text(size = 12),
        legend.position = "right")

# ---------- Fig2a: DEG volcano plots ----------
plot_volcano <- function(deg_df, title, up_color = "#B2182B", down_color = "#2166AC") {
  d <- deg_df
  d$sig <- "NS"
  d$sig[d$logFC >= 1 & d$adj.P.Val < 0.05] <- "Up"
  d$sig[d$logFC <= -1 & d$adj.P.Val < 0.05] <- "Down"
  ggplot(d, aes(x = logFC, y = -log10(P.Value), color = sig)) +
    geom_point(size = 0.5, alpha = 0.6) +
    scale_color_manual(values = c(Down = down_color, NS = "grey70", Up = up_color),
                       name = NULL) +
    geom_vline(xintercept = c(-1, 1), linetype = 2, color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey40") +
    labs(x = expression(log[2]~fold~change), y = expression(-log[10]~P),
         title = title) +
    theme_pub
}

# Panel letters A/B/C/D correspond to the in-text citations (Fig. 2A-C, Fig. 2D)
fig2_letters <- c("GSE20347" = "A", "GSE45670" = "B", "TCGA" = "C")
for (cohort in c("GSE20347", "GSE45670", "TCGA")) {
  d <- fread(file.path(out_dir, paste0("DEG_", cohort, ".tsv")), data.table = FALSE)
  p <- plot_volcano(d, paste0(fig2_letters[[cohort]], " ", cohort, ": Tumor vs Normal"))
  ggsave(file.path(fig_dir, paste0("fig2_volcano_", cohort, ".png")), p,
         width = 7, height = 6, dpi = 300)
}

# ---------- Fig2b: consensus-DEG Venn ----------
d1 <- fread(file.path(out_dir, "DEG_GSE20347.tsv"), data.table = FALSE)
d2 <- fread(file.path(out_dir, "DEG_GSE45670.tsv"), data.table = FALSE)
d3 <- fread(file.path(out_dir, "DEG_TCGA.tsv"), data.table = FALSE)
sig1 <- d1$gene[abs(d1$logFC) >= 1 & d1$adj.P.Val < 0.05]
sig2 <- d2$gene[abs(d2$logFC) >= 1 & d2$adj.P.Val < 0.05]
sig3 <- d3$gene[abs(d3$logFC) >= 1 & d3$adj.P.Val < 0.05]

png(file.path(fig_dir, "fig2_venn.png"), width = 8, height = 7, units = "in", res = 300)
v <- venn.diagram(
  x = list("GSE20347" = sig1, "GSE45670" = sig2, "TCGA" = sig3),
  filename = NULL, fill = c("#2166AC", "#B2182B", "#1B7837"),
  alpha = 0.5, cex = 1.5, cat.cex = 1.5, main = "D Consensus DEGs across cohorts")
grid.draw(v)
dev.off()

# ---------- Fig3: LASSO coefficient path + CV ----------
cv_fit <- readRDS(file.path(out_dir, "lasso_cv.rds"))
model <- readRDS(file.path(out_dir, "lasso_model.rds"))

p1 <- ggplot(data.frame(lambda = cv_fit$lambda,
                        cvm = cv_fit$cvm, cvup = cv_fit$cvup, cvlo = cv_fit$cvlo),
             aes(x = log(lambda))) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup), alpha = 0.3, fill = "grey60") +
  geom_line(aes(y = cvm), color = "#B2182B", linewidth = 0.8) +
  geom_vline(xintercept = log(model$lambda_min), linetype = 2, color = "#2166AC") +
  annotate("text", x = log(model$lambda_min), y = max(cv_fit$cvup),
           label = "lambda.min", hjust = -0.1, color = "#2166AC") +
  labs(x = expression(log(lambda)), y = "Partial likelihood deviance",
       title = "Cross-validation curve") +
  theme_pub

p2 <- ggplot(data.frame(lambda = cv_fit$lambda, ncoef = cv_fit$nzero),
             aes(x = log(lambda), y = ncoef)) +
  geom_line(color = "#1B7837", linewidth = 0.8) +
  geom_vline(xintercept = log(model$lambda_min), linetype = 2, color = "#2166AC") +
  labs(x = expression(log(lambda)), y = "Number of non-zero coefficients",
       title = "Coefficient shrinkage path") +
  theme_pub

p_comb <- p1 / p2
ggsave(file.path(fig_dir, "fig3_lasso.png"), p_comb, width = 7, height = 9, dpi = 300)

# ---------- Fig4: model-gene forest plot (TCGA univariable Cox) ----------
prog <- fread(file.path(out_dir, "prognostic_uni_cox.tsv"), data.table = FALSE)
mg <- model$genes
d <- prog[prog$gene %in% mg, ]
d <- d[order(d$HR), ]
d$gene <- factor(d$gene, levels = d$gene)
p <- ggplot(d, aes(x = HR, y = gene)) +
  geom_vline(xintercept = 1, linetype = 2, color = "grey40") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.15,
                 color = "#2166AC") +
  geom_point(aes(size = -log10(P)), color = "#B2182B") +
  labs(x = "Hazard ratio (95% CI)", y = NULL,
       title = "Univariate Cox regression of model genes (TCGA)") +
  theme_pub
ggsave(file.path(fig_dir, "fig4_forest.png"), p, width = 8, height = 7, dpi = 300)

# ---------- Fig5: risk-score heatmap + survival status ----------
risk_df <- fread(file.path(out_dir, "risk_score_tcga.tsv"), data.table = FALSE)
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
model_genes <- model$genes
risk_df <- risk_df[order(risk_df$risk_score), ]
heat <- expr[model_genes, risk_df$case_id, drop = FALSE]
# Unify colour semantics by coefficient direction: high expression of positive-coefficient genes = high risk
coefs <- model$coef[, 1]
heat_signed <- sweep(heat, 1, coefs, "*")
anno <- data.frame(Risk = factor(risk_df$risk_group, levels = c("Low", "High")),
                   Survival = factor(risk_df$os_status, levels = c(0, 1),
                                     labels = c("Alive", "Dead")))
rownames(anno) <- risk_df$case_id
ann_colors <- list(Risk = c(Low = "#2166AC", High = "#B2182B"),
                   Survival = c(Alive = "#1B7837", Dead = "grey30"))

suppressPackageStartupMessages(library(pheatmap))
pheatmap(heat_signed, annotation_col = anno, annotation_colors = ann_colors,
         show_colnames = FALSE, show_rownames = TRUE,
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         main = "Risk score heatmap (TCGA)",
         filename = file.path(fig_dir, "fig5_risk_heatmap.png"),
         width = 10, height = 7, dpi = 300)

cat("10_figures done\n")
