# ============================================================================
# 08_immune_infiltration.R
# Purpose: marker-gene-based ssGSEA (GSVA) to assess infiltration of 28 immune cell types,
#       comparing immune profiles between high/low-risk groups and correlations with the risk score
# Input: results/tcga_expr_symbol.rds, results/risk_score_tcga.tsv
# Output: results/immune_ssgsea.tsv          per-sample x 28 immune-cell scores
#        results/immune_high_vs_low.tsv     between-group comparison
#        results/immune_cor_risk.tsv        correlation with risk score
#        results/immune_heatmap.png         heatmap
#        results/immune_boxplot.png         representative immune-cell boxplot
#        results/immune_cor_risk.png        correlation scatter plot
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(GSVA); library(data.table); library(ggplot2); library(pheatmap)
})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")

# ---------- 1. Marker gene sets for 28 immune cell types (Charoentong et al., Cancer Res 2012) ----------
markers <- list(
  "Activated CD4 T" = c("CD40LG", "IL2RB"),
  "Activated CD8 T" = c("CD8A", "CD8B", "CTLA4", "GZMK", "PRF1"),
  "Activated dendritic" = c("CD80", "CD83", "CD86", "NRP1", "ITGAX"),
  "CD56bright NK" = c("KIR2DL1", "KLRB1", "NCAM1", "NCR1"),
  "CD56dim NK" = c("CX3CR1", "GNLY", "KIR2DL3", "KLRD1", "NCAM1", "NCR1"),
  "Central memory CD4 T" = c("CCR7", "IL7R", "SELL"),
  "Central memory CD8 T" = c("CCR7", "IL7R", "SELL", "TCF7"),
  "Effector memory CD4 T" = c("CCR6", "IL4R", "TNFRSF17"),
  "Effector memory CD8 T" = c("CCR5", "CXCR3", "GZMH", "KLF3", "S1PR5"),
  "Eosinophil" = c("CCR3", "IL5RA", "RNASE2", "SEMG1"),
  "Gamma delta T" = c("GZMK", "KLRC1", "KLRC3", "KLRG1"),
  "Immature B" = c("CD22", "IGHD", "MS4A1", "TNFRSF17"),
  "Immature dendritic" = c("CD1C", "CLEC10A", "FCER1A", "GZMB"),
  "Macrophage" = c("CD68", "CD163", "CSF1R", "MSR1"),
  "Mast" = c("CPA3", "MS4A2", "TPSAB1", "TPSB2"),
  "MDSC" = c("ARG1", "CD33", "CD68", "ITGAM", "S100A8", "S100A9"),
  "Memory B" = c("AIM2", "BANK1", "CD27", "CR2", "TNFRSF13B"),
  "Monocyte" = c("CD14", "CD36", "CSF1R", "FCGR2A", "FCGR3A", "MS4A7"),
  "Natural killer" = c("FCGR3A", "GNLY", "KLRB1", "NCAM1", "NCR1"),
  "Neutrophil" = c("CSF3R", "CXCR1", "CXCR2", "FCGR3B", "S100A8"),
  "Plasmacytoid dendritic" = c("GZMB", "IL3RA", "NRP1"),
  "Regulatory T" = c("CCR4", "FOXP3", "IL2RA", "TGFBR1"),
  "T follicular helper" = c("BCL6", "CXCR5", "ICOS", "IL21"),
  "T gamma delta" = c("KLRC1", "KLRC3", "KLRD1", "TRDC"),
  "Type 1 T helper" = c("CCR5", "IFNG", "IL12RB2", "STAT4", "TBX21"),
  "Type 2 T helper" = c("CCR4", "GATA3", "IL4", "IL5", "IL13"),
  "Type 17 T helper" = c("CCR6", "IL17A", "IL17F", "IL21", "RORC"),
  "Regulatory B" = c("CD1D", "IL10", "TGFB1")
)

# Keep genes present in the expression matrix
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
avail <- rownames(expr)
markers <- lapply(markers, function(g) intersect(g, avail))
markers <- markers[sapply(markers, length) >= 2]
cat("Available immune cell types:", length(markers), "\n")

# ---------- 2. ssGSEA (GSVA 2.6 API) ----------
gs_param <- ssgseaParam(as.matrix(expr), markers)
gs <- gsva(gs_param)
cat("ssGSEA matrix:", dim(gs), "\n")
gs_t <- t(gs)   # samples x immune cells

# ---------- 3. Merge risk groups ----------
risk_df <- fread(file.path(out_dir, "risk_score_tcga.tsv"), data.table = FALSE)
m <- match(rownames(gs_t), risk_df$case_id)
keep <- !is.na(m)
gs_t <- gs_t[keep, , drop = FALSE]
risk_df <- risk_df[m[keep], ]
cat("Samples for immune analysis:", nrow(gs_t), "\n")

immune_df <- data.frame(case_id = rownames(gs_t), risk_group = risk_df$risk_group,
                        risk_score = risk_df$risk_score, gs_t, check.names = FALSE)
fwrite(immune_df, file.path(out_dir, "immune_ssgsea.tsv"), sep = "\t")

# ---------- 4. High vs low-risk group comparison (Wilcoxon) ----------
cmp <- data.frame(cell = colnames(gs_t), mean_low = NA_real_, mean_high = NA_real_,
                  p = NA_real_)
for (i in seq_along(gs_t |> colnames())) {
  ct <- colnames(gs_t)[i]
  lo <- gs_t[risk_df$risk_group == "Low", ct]
  hi <- gs_t[risk_df$risk_group == "High", ct]
  cmp$mean_low[i] <- mean(lo); cmp$mean_high[i] <- mean(hi)
  cmp$p[i] <- wilcox.test(lo, hi)$p.value
}
cmp$p.adj <- p.adjust(cmp$p, method = "BH")
cmp <- cmp[order(cmp$p), ]
fwrite(cmp, file.path(out_dir, "immune_high_vs_low.tsv"), sep = "\t")
cat("\nImmune-infiltration differences between high and low-risk groups (P<0.05):\n")
print(cmp[cmp$p < 0.05, ])

# ---------- 5. Correlation with risk score (Spearman) ----------
cor_res <- data.frame(cell = colnames(gs_t), rho = NA_real_, p = NA_real_)
for (i in seq_along(gs_t |> colnames())) {
  ct <- colnames(gs_t)[i]
  ct_res <- suppressWarnings(cor.test(gs_t[, ct], risk_df$risk_score,
                                      method = "spearman", exact = FALSE))
  cor_res$rho[i] <- ct_res$estimate
  cor_res$p[i] <- ct_res$p.value
}
cor_res$p.adj <- p.adjust(cor_res$p, method = "BH")
cor_res <- cor_res[order(cor_res$p), ]
fwrite(cor_res, file.path(out_dir, "immune_cor_risk.tsv"), sep = "\t")
cat("\nImmune cells correlated with the risk score (Spearman, P<0.05):\n")
print(cor_res[cor_res$p < 0.05, ])

# ---------- 6. Visualisation ----------
# 6.1 Heatmap: samples x immune cells ordered by risk score (groups annotated)
anno <- data.frame(Risk = factor(risk_df$risk_group, levels = c("Low", "High")))
rownames(anno) <- rownames(gs_t)
pheatmap(gs_t, annotation_col = anno,
         annotation_colors = list(Risk = c(Low = "#2166AC", High = "#B2182B")),
         show_colnames = FALSE, show_rownames = FALSE,
         clustering_method = "ward.D2",
         main = "Immune infiltration landscape (ssGSEA, 27 cell types)",
         filename = file.path(out_dir, "immune_heatmap.png"),
         width = 10, height = 8, dpi = 300)

# 6.2 Representative immune-cell boxplots (6 most different between groups)
top_cells <- cmp$cell[1:min(6, sum(cmp$p < 0.05))]
if (length(top_cells) >= 1) {
  df_long <- data.frame(
    cell = factor(rep(top_cells, each = nrow(gs_t)),
                  levels = top_cells),
    value = as.vector(gs_t[, top_cells]),
    group = rep(risk_df$risk_group, length(top_cells)))
  p <- ggplot(df_long, aes(x = group, y = value, fill = group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.85) +
    geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
    facet_wrap(~ cell, scales = "free_y") +
    scale_fill_manual(values = c(Low = "#2166AC", High = "#B2182B")) +
    labs(x = NULL, y = "ssGSEA score",
         title = "Immune infiltration differs between risk groups") +
    theme_bw(base_size = 11) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position = "bottom")
  ggsave(file.path(out_dir, "immune_boxplot.png"), p, width = 11, height = 7, dpi = 300)
}

# 6.3 Bar plot of correlations with the risk score (top 10 |rho|)
top_cor <- cor_res[order(abs(cor_res$rho), decreasing = TRUE), ][1:min(10, nrow(cor_res)), ]
top_cor$direction <- ifelse(top_cor$rho > 0, "Positive", "Negative")
p <- ggplot(top_cor, aes(x = reorder(cell, rho), y = rho, fill = direction)) +
  geom_col(alpha = 0.9) +
  coord_flip() +
  scale_fill_manual(values = c(Positive = "#B2182B", Negative = "#2166AC")) +
  labs(x = NULL, y = "Spearman rho", fill = "Correlation direction",
       title = "Correlation of immune cells with risk score (top 10)") +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"))
ggsave(file.path(out_dir, "immune_cor_risk.png"), p, width = 9, height = 6, dpi = 300)

cat("\n08_immune_infiltration done\n")
