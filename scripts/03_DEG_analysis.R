# ============================================================================
# 03_DEG_analysis.R
# Purpose: multi-cohort integrative ESCC tumour vs normal differential expression analysis
#   Cohorts: GSE20347 (17v17 paired), GSE45670 (28v10), TCGA-ESCC (95v3)
#   Method: limma; consensus DEG = significant in >=2 data sources with consistent direction
# Output: results/DEG_<cohort>.tsv        full DEG table per cohort
#        results/DEG_consensus.tsv       consensus DEGs (with per-cohort statistics)
#        results/DEG_summary.txt         summary statistics
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(limma); library(data.table); library(matrixStats)
})

base <- "/workspace/escc_bioinformatics"
out_dir <- file.path(base, "results")

deg_list <- list()   # one data.frame per data source

# ---------- 1. GSE20347 (paired design) ----------
g1 <- readRDS(file.path(out_dir, "geo_GSE20347.rds"))
e1 <- g1$expr
ph1 <- g1$pheno
# Extract patient IDs from titles to establish pairing
ph1$patient <- sub(".*patient (E[0-9]+).*", "\\1", ph1$title)
ph1$patient[!grepl("^E[0-9]+$", ph1$patient)] <- ph1$gsm[!grepl("^E[0-9]+$", ph1$patient)]
stopifnot(all(grepl("^E[0-9]+$", ph1$patient)))  # confirm pairing extraction succeeded
cat("GSE20347 paired patients:", length(unique(ph1$patient)), "\n")
design1 <- model.matrix(~ ph1$patient + ph1$tissue)
colnames(design1) <- make.names(colnames(design1))
coef1 <- grep("tissueTumor", colnames(design1), value = TRUE)
stopifnot(length(coef1) == 1)
fit1 <- lmFit(e1, design1)
contrast1 <- makeContrasts(get(coef1), levels = design1)
fit1c <- contrasts.fit(fit1, contrast1)
fit1c <- eBayes(fit1c)
tt1 <- topTable(fit1c, number = Inf, sort.by = "P")
deg_list$GSE20347 <- data.frame(gene = rownames(tt1), logFC = tt1$logFC,
                                P.Value = tt1$P.Value, adj.P.Val = tt1$adj.P.Val)
cat("GSE20347: genes", nrow(tt1), "| significant (|logFC|>1 & FDR<0.05):",
    sum(abs(tt1$logFC) > 1 & tt1$adj.P.Val < 0.05), "\n")

# ---------- 2. GSE45670 (unpaired) ----------
g2 <- readRDS(file.path(out_dir, "geo_GSE45670.rds"))
e2 <- g2$expr
ph2 <- g2$pheno
design2 <- model.matrix(~ ph2$tissue)
colnames(design2) <- make.names(colnames(design2))
coef2 <- grep("tissueTumor", colnames(design2), value = TRUE)
stopifnot(length(coef2) == 1)
fit2 <- lmFit(e2, design2)
fit2 <- eBayes(fit2)
tt2 <- topTable(fit2, coef = coef2, number = Inf, sort.by = "P")
deg_list$GSE45670 <- data.frame(gene = rownames(tt2), logFC = tt2$logFC,
                                P.Value = tt2$P.Value, adj.P.Val = tt2$adj.P.Val)
cat("GSE45670: genes", nrow(tt2), "| significant (|logFC|>1 & FDR<0.05):",
    sum(abs(tt2$logFC) > 1 & tt2$adj.P.Val < 0.05), "\n")

# ---------- 3. GSE75241 (if already processed) ----------
if (file.exists(file.path(out_dir, "geo_GSE75241.rds"))) {
  g3 <- readRDS(file.path(out_dir, "geo_GSE75241.rds"))
  e3 <- g3$expr; ph3 <- g3$pheno
  ph3$patient <- sub(".*patient ([0-9]+).*", "\\1", ph3$title)
  if (all(grepl("^[0-9]+$", ph3$patient)) && length(unique(ph3$patient)) > 10) {
    design3 <- model.matrix(~ ph3$patient + ph3$tissue)
  } else {
    design3 <- model.matrix(~ ph3$tissue)
  }
  fit3 <- lmFit(e3, design3)
  coef_name <- grep("tissueTumor", colnames(design3), value = TRUE)
  fit3 <- eBayes(fit3)
  tt3 <- topTable(fit3, coef = coef_name, number = Inf, sort.by = "P")
  deg_list$GSE75241 <- data.frame(gene = rownames(tt3), logFC = tt3$logFC,
                                  P.Value = tt3$P.Value, adj.P.Val = tt3$adj.P.Val)
  cat("GSE75241: genes", nrow(tt3), "| significant:",
      sum(abs(tt3$logFC) > 1 & tt3$adj.P.Val < 0.05), "\n")
} else {
  cat("GSE75241: not processed; skipped\n")
}

# ---------- 4. TCGA-ESCC (95 tumour vs 3 normal, supporting evidence) ----------
# Use the gene-symbol-space matrix (converted by 01b; harmonised with the GEO cohorts)
tc <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
ne <- readRDS(file.path(out_dir, "tcga_normal_symbol.rds"))
genes <- intersect(rownames(tc), rownames(ne))
et <- cbind(tc[genes, ], ne[genes, ])
grp <- factor(c(rep("Tumor", ncol(tc)), rep("Normal", ncol(ne))),
              levels = c("Normal", "Tumor"))
design_t <- model.matrix(~ grp)
fit_t <- lmFit(et, design_t)
fit_t <- eBayes(fit_t)
tt_t <- topTable(fit_t, coef = "grpTumor", number = Inf, sort.by = "P")
deg_list$TCGA <- data.frame(gene = rownames(tt_t), logFC = tt_t$logFC,
                            P.Value = tt_t$P.Value, adj.P.Val = tt_t$adj.P.Val)
cat("TCGA: genes", nrow(tt_t), "| significant:", sum(abs(tt_t$logFC) > 1 & tt_t$adj.P.Val < 0.05), "\n")

# ---------- 5. Consensus DEGs ----------
FC_CUT <- 1; FDR_CUT <- 0.05
sig_sets <- lapply(deg_list, function(d)
  d$gene[abs(d$logFC) >= FC_CUT & d$adj.P.Val < FDR_CUT])
names(sig_sets) <- names(deg_list)
all_genes <- unique(unlist(lapply(deg_list, function(d) d$gene)))
n_sources <- length(sig_sets)

consensus <- data.frame(gene = all_genes)
for (nm in names(deg_list)) {
  d <- deg_list[[nm]]
  consensus[[paste0(nm, "_logFC")]] <- d$logFC[match(consensus$gene, d$gene)]
  consensus[[paste0(nm, "_P")]] <- d$P.Value[match(consensus$gene, d$gene)]
  consensus[[paste0(nm, "_FDR")]] <- d$adj.P.Val[match(consensus$gene, d$gene)]
}
# Count of significant sources and direction consistency
fc_cols <- paste0(names(deg_list), "_logFC")
fdr_cols <- paste0(names(deg_list), "_FDR")
sig_mat <- mapply(function(f, q) abs(consensus[[f]]) >= FC_CUT & consensus[[q]] < FDR_CUT,
                  fc_cols, fdr_cols)
colnames(sig_mat) <- names(deg_list)
sig_mat[is.na(sig_mat)] <- FALSE   # a gene missing from a cohort counts as non-significant
consensus$n_sig <- rowSums(sig_mat)
consensus$mean_logFC <- rowMeans(consensus[, fc_cols], na.rm = TRUE)
# Direction consistency: all significant sources have the same logFC sign (significant sources only, so non-significant cohorts cannot interfere)
fc_mat <- as.matrix(consensus[, fc_cols])
sig_fc_mat <- fc_mat
sig_fc_mat[!sig_mat] <- NA
sign_mat <- sign(sig_fc_mat)
consensus$direction_consistent <- apply(sign_mat, 1, function(s) {
  s <- s[!is.na(s)]
  length(s) > 0 && all(s == s[1])
})
# Direction is set by the first significant source (all significant sources are confirmed same-sign at this point)
first_sig <- apply(sign_mat, 1, function(s) {
  s <- s[!is.na(s)]
  if (length(s) == 0) NA_real_ else s[1]
})
consensus$direction <- ifelse(first_sig > 0, "Up_in_Tumor",
                       ifelse(first_sig < 0, "Down_in_Tumor", NA_character_))

# Consensus criteria: significant in >=2 data sources, consistent direction, and >=1 GEO cohort
geo_names <- intersect(c("GSE20347", "GSE45670", "GSE75241"), names(deg_list))
consensus$is_consensus <- consensus$n_sig >= 2 & consensus$direction_consistent &
  rowSums(sig_mat[, geo_names, drop = FALSE]) >= 1

cons_degs <- consensus[consensus$is_consensus, ]
cat("\n=== Consensus DEGs ===\n")
cat("Total genes:", nrow(consensus), "| consensus DEGs:", nrow(cons_degs),
    "| up in tumour:", sum(cons_degs$direction == "Up_in_Tumor"),
    "| down in tumour:", sum(cons_degs$direction == "Down_in_Tumor"), "\n")

# ---------- 6. Save ----------
for (nm in names(deg_list)) {
  fwrite(deg_list[[nm]], file.path(out_dir, paste0("DEG_", nm, ".tsv")), sep = "\t")
}
fwrite(consensus, file.path(out_dir, "DEG_consensus.tsv"), sep = "\t")
fwrite(cons_degs, file.path(out_dir, "DEG_consensus_sig.tsv"), sep = "\t")
saveRDS(deg_list, file.path(out_dir, "DEG_list.rds"))
saveRDS(consensus, file.path(out_dir, "DEG_consensus.rds"))

summary_lines <- c(
  "===== ESCC differential expression analysis summary =====",
  paste0("Data sources: ", paste(names(deg_list), collapse = ", ")),
  paste0("Thresholds: |log2FC| >= ", FC_CUT, ", FDR < ", FDR_CUT),
  paste0("Consensus criteria: >=", 2, " data sources significant with consistent direction (including >=1 GEO cohort)"),
  ""
)
for (nm in names(deg_list)) {
  d <- deg_list[[nm]]
  summary_lines <- c(summary_lines,
    paste0(nm, ": significant ", sum(abs(d$logFC) >= FC_CUT & d$adj.P.Val < FDR_CUT),
           " (up ", sum(d$logFC >= FC_CUT & d$adj.P.Val < FDR_CUT),
           ", down ", sum(d$logFC <= -FC_CUT & d$adj.P.Val < FDR_CUT), ")"))
}
summary_lines <- c(summary_lines,
  paste0("Consensus DEGs: ", nrow(cons_degs), " (up in tumour ", sum(cons_degs$direction == "Up_in_Tumor"),
         ", down in tumour ", sum(cons_degs$direction == "Down_in_Tumor"), ")"))
writeLines(summary_lines, file.path(out_dir, "DEG_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")
cat("\n03_DEG_analysis done\n")
