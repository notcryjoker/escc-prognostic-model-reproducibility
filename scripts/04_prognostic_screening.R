# ============================================================================
# 04_prognostic_screening.R
# Purpose: univariable Cox regression of consensus DEGs in the TCGA-ESCC cohort to screen prognosis-related genes
# Input: results/tcga_expr_tpm.rds, tcga_clinical_surv.rds, DEG_consensus.rds
# Output: results/prognostic_uni_cox.tsv   univariable Cox results for all consensus DEGs
#        results/prognostic_genes.tsv     screened prognostic genes (P<0.05)
#       results/prognostic_summary.txt
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(survival); library(data.table)
})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")

# ---------- 1. Read data ----------
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))  # gene symbols x samples, log2(TPM+1)
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
consensus <- readRDS(file.path(out_dir, "DEG_consensus.rds"))
cons_degs <- consensus$gene[consensus$is_consensus]
cat("Consensus DEGs:", length(cons_degs), "\n")

# Align samples
samples <- intersect(colnames(expr), rownames(clin))
expr <- expr[, samples]
clin <- clin[samples, ]
surv_obj <- Surv(clin$os_time, clin$os_status)
cat("Samples for Cox analysis:", length(samples), "| death events:", sum(clin$os_status), "\n")

# ---------- 2. Univariable Cox (gene by gene) ----------
genes <- intersect(cons_degs, rownames(expr))
cat("Analysable genes (present in the expression matrix):", length(genes), "\n")

res <- data.frame(gene = genes, HR = NA_real_, CI_lower = NA_real_,
                  CI_upper = NA_real_, coef = NA_real_, se = NA_real_,
                  z = NA_real_, P = NA_real_)

for (i in seq_along(genes)) {
  g <- genes[i]
  x <- as.numeric(expr[g, ])
  # Remove zero-variance genes (not expected among consensus DEGs)
  if (sd(x) == 0) next
  fit <- try(coxph(surv_obj ~ x), silent = TRUE)
  if (!inherits(fit, "try-error")) {
    s <- summary(fit)
    # s$coef columns: coef, exp(coef), se(coef), z, Pr(>|z|)
    res$coef[i] <- s$coef[1, 1]
    res$HR[i]   <- s$coef[1, 2]
    res$se[i]   <- s$coef[1, 3]
    res$z[i]    <- s$coef[1, 4]
    res$P[i]    <- s$coef[1, 5]
    # s$conf.int columns: exp(coef), exp(-coef), lower .95, upper .95
    res$CI_lower[i] <- s$conf.int[1, 3]
    res$CI_upper[i] <- s$conf.int[1, 4]
  }
}

# BH adjustment
res$FDR <- p.adjust(res$P, method = "BH")
res$direction <- ifelse(res$HR > 1, "risk", "protective")

n_uni <- sum(res$P < 0.05, na.rm = TRUE)
n_uni_fdr <- sum(res$FDR < 0.05, na.rm = TRUE)
cat("Univariable Cox significant (P<0.05):", n_uni, "\n")
cat("Univariable Cox significant (FDR<0.05):", n_uni_fdr, "\n")

# ---------- 3. Screen prognostic genes ----------
# Use P<0.05 (the conventional univariable screening threshold) and also report FDR
P_CUT <- 0.05
prog <- res[!is.na(res$P) & res$P < P_CUT, ]
prog <- prog[order(prog$P), ]
cat("Prognostic genes after screening:", nrow(prog),
    " (risk genes:", sum(prog$direction == "risk"),
    "| protective genes:", sum(prog$direction == "protective"), ")\n")

# ---------- 4. Save ----------
fwrite(res, file.path(out_dir, "prognostic_uni_cox.tsv"), sep = "\t")
fwrite(prog, file.path(out_dir, "prognostic_genes.tsv"), sep = "\t")
saveRDS(prog, file.path(out_dir, "prognostic_genes.rds"))

lines <- c(
  "===== ESCC prognostic-gene univariable Cox screening =====",
  paste0("Cohort: TCGA-ESCC (n=", length(samples), ", events=", sum(clin$status), ")"),
  paste0("Candidate pool: consensus DEGs ", length(cons_degs)),
  paste0("Univariable Cox P<0.05: ", n_uni),
  paste0("Univariable Cox FDR<0.05: ", n_uni_fdr),
  paste0("Prognostic genes included (P<", P_CUT, "): ", nrow(prog)),
  paste0("  Risk genes HR>1: ", sum(prog$direction == "risk")),
  paste0("  Protective genes HR<1: ", sum(prog$direction == "protective")),
  "",
  "Top 15 prognostic genes (sorted by P):",
  paste(sprintf("  %-12s HR=%.2f (%.2f-%.2f)  P=%.2e  FDR=%.3f  %s",
                prog$gene[1:min(15, nrow(prog))], prog$HR[1:min(15, nrow(prog))],
                prog$CI_lower[1:min(15, nrow(prog))], prog$CI_upper[1:min(15, nrow(prog))],
                prog$P[1:min(15, nrow(prog))], prog$FDR[1:min(15, nrow(prog))],
                prog$direction[1:min(15, nrow(prog))]), collapse = "\n")
)
writeLines(lines, file.path(out_dir, "prognostic_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
cat("\n04_prognostic_screening done\n")
