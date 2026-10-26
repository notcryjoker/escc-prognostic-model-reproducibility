# ============================================================================
# 05_lasso_model.R
# Purpose: build an ESCC prognostic risk model with LASSO-Cox from univariable-screened prognostic genes
# Input: results/tcga_expr_symbol.rds, tcga_clinical_surv.rds,
#       results/prognostic_genes.rds
# Output: results/lasso_model.rds        model object (gene coefficients)
#        results/lasso_cv.rds          CV results
#        results/risk_score_tcga.tsv   TCGA cohort risk scores
#       results/lasso_summary.txt
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(data.table)
})

set.seed(123)
base <- "/workspace/escc_bioinformatics"
out_dir <- file.path(base, "results")

# ---------- 1. Data preparation ----------
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
prog <- readRDS(file.path(out_dir, "prognostic_genes.rds"))

genes <- intersect(prog$gene, rownames(expr))
cat("LASSO input genes:", length(genes), "\n")

samples <- intersect(colnames(expr), rownames(clin))
# Keep only samples with complete survival data
complete <- complete.cases(clin[samples, c("os_time", "os_status")])
samples <- samples[complete]
x <- t(expr[genes, samples, drop = FALSE])          # samples x genes
clin <- clin[samples, ]
surv_obj <- Surv(clin$os_time, clin$os_status)
cat("Modelling samples:", nrow(x), "| events:", sum(clin$os_status), "\n")

# ---------- 2. LASSO-Cox (10-fold cross-validation) ----------
cv_fit <- cv.glmnet(x, surv_obj, family = "cox", alpha = 1,
                    nfolds = 10, type.measure = "C", cox.ties = "breslow")
lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se
cat("lambda.min:", round(lambda_min, 4), "| lambda.1se:", round(lambda_1se, 4), "\n")

# Extract genes with non-zero coefficients (lambda.min)
coef_mat <- as.matrix(coef(cv_fit, s = lambda_min))
coef_mat <- coef_mat[coef_mat[, 1] != 0, , drop = FALSE]
model_genes <- rownames(coef_mat)
cat("Model genes selected (lambda.min):", length(model_genes), "\n")
print(round(coef_mat, 4))

# ---------- 3. Compute risk scores ----------
risk_score <- as.numeric(x[, model_genes, drop = FALSE] %*% coef_mat[, 1])
risk_df <- data.frame(case_id = rownames(x), risk_score = risk_score,
                      os_time = clin$os_time, os_status = clin$os_status,
                      stringsAsFactors = FALSE)
# Group at the median
med <- median(risk_score)
risk_df$risk_group <- ifelse(risk_score > med, "High", "Low")
cat("\nMedian risk score:", round(med, 3), "\n")
cat("High/low-risk group samples:", sum(risk_df$risk_group == "High"), "/",
    sum(risk_df$risk_group == "Low"), "\n")

# ---------- 4. Save ----------
saveRDS(list(coef = coef_mat, genes = model_genes, lambda_min = lambda_min,
             lambda_1se = lambda_1se, median = med),
        file.path(out_dir, "lasso_model.rds"))
saveRDS(cv_fit, file.path(out_dir, "lasso_cv.rds"))
fwrite(risk_df, file.path(out_dir, "risk_score_tcga.tsv"), sep = "\t")

lines <- c(
  "===== LASSO-Cox prognostic model construction =====",
  paste0("Input genes: ", length(genes), " (univariable Cox P<0.05)"),
  paste0("Modelling cohort: TCGA-ESCC, n=", nrow(x), ", events=", sum(clin$os_status)),
  paste0("lambda.min: ", round(lambda_min, 4)),
  paste0("lambda.1se: ", round(lambda_1se, 4)),
  paste0("Model genes selected: ", length(model_genes), " -> ",
         paste(model_genes, collapse = ", ")),
  "",
  "Gene coefficients:",
  paste(sprintf("  %-10s %.4f", model_genes, coef_mat[, 1]), collapse = "\n"),
  "",
  paste0("Median risk score: ", round(med, 3))
)
writeLines(lines, file.path(out_dir, "lasso_summary.txt"))
cat(paste(lines, collapse = "\n"), "\n")
cat("\n05_lasso_model done\n")
