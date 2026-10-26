# ============================================================================
# 09_external_validation.R
# Purpose: validate the 17-gene prognostic model in the independent GSE53624 cohort
#   Risk scores reuse the gene coefficients and median grouping cut-off from the TCGA training set
# Input: results/geo_GSE53624.rds, results/lasso_model.rds
# Output: results/validation_GSE53624.txt
#       results/km_GSE53624.png, results/roc_curves_GSE53624.png
#       results/risk_score_gse53624.tsv
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(survival); library(timeROC); library(data.table); library(survminer)
})

base <- "/workspace/escc_bioinformatics"
out_dir <- file.path(base, "results")

# ---------- 1. Data and model ----------
gse <- readRDS(file.path(out_dir, "geo_GSE53624.rds"))
model <- readRDS(file.path(out_dir, "lasso_model.rds"))
expr <- gse$expr
clin <- gse$pheno
model_genes <- model$genes
coefs <- model$coef[, 1]
med_cutoff <- model$median   # reuse the TCGA median cut-off

# ---------- 2. Compute risk scores ----------
genes <- intersect(model_genes, rownames(expr))
cat("Model genes available:", length(genes), "/", length(model_genes), "\n")
if (length(genes) < length(model_genes)) {
  cat("Warning: missing genes", paste(setdiff(model_genes, genes), collapse = ", "),
      " (scores computed from available genes)\n")
}
x <- t(expr[genes, , drop = FALSE])
coef_sub <- coefs[genes]
risk_score <- as.numeric(x %*% coef_sub)

risk_df <- data.frame(case_id = colnames(expr), risk_score = risk_score,
                      stringsAsFactors = FALSE)
# Align survival data
m <- match(risk_df$case_id, rownames(clin))
risk_df$os_time <- clin$os_time[m]
risk_df$os_status <- clin$os_status[m]
risk_df$stage <- clin$stage[m]
risk_df$age <- clin$Age[m]
complete <- complete.cases(risk_df[, c("os_time", "os_status")])
risk_df <- risk_df[complete, ]
risk_df$risk_group <- ifelse(risk_df$risk_score > med_cutoff, "High", "Low")
cat("External validation samples:", nrow(risk_df), "| events:", sum(risk_df$os_status),
    "| high/low:", sum(risk_df$risk_group == "High"), "/", sum(risk_df$risk_group == "Low"), "\n")

surv_obj <- with(risk_df, Surv(os_time, os_status))
val_lines <- c("===== GSE53624 external validation =====",
               paste0("n = ", nrow(risk_df), ", events = ", sum(risk_df$os_status)),
               paste0("Grouping cut-off (TCGA median): ", round(med_cutoff, 3)))

# ---------- 3. KM analysis ----------
surv_diff <- survdiff(surv_obj ~ risk_group, data = risk_df)
p_km <- pchisq(surv_diff$chisq, df = 1, lower.tail = FALSE)
val_lines <- c(val_lines,
  paste0("KM log-rank P = ", format.pval(p_km, digits = 3)))

png(file.path(out_dir, "km_GSE53624.png"), width = 7, height = 6, units = "in", res = 300)
km_fit <- survfit(surv_obj ~ risk_group, data = risk_df)
print(ggsurvplot(km_fit, data = risk_df, pval = TRUE, risk.table = TRUE,
                 palette = c("#2166AC", "#B2182B"),
                 xlab = "Overall survival time (days)", ylab = "Survival probability",
                 title = "Kaplan-Meier curves by risk group (GSE53624)",
                 legend.title = "Risk group", font.title = 12))
dev.off()

# ---------- 4. Time-dependent ROC ----------
roc_res <- timeROC(T = risk_df$os_time, delta = risk_df$os_status,
                   marker = risk_df$risk_score,
                   cause = 1, times = c(365, 1095, 1825))
aucs <- sapply(c("t=365", "t=1095", "t=1825"),
               function(t) round(roc_res$AUC[[t]], 3))
val_lines <- c(val_lines,
  paste0("Time-dependent AUC: 1-year = ", aucs[1], ", 3-year = ", aucs[2], ", 5-year = ", aucs[3]))

png(file.path(out_dir, "roc_curves_GSE53624.png"), width = 7, height = 6, units = "in", res = 300)
plot(roc_res$FP[, "t=365"], roc_res$TP[, "t=365"], type = "l", col = "#2166AC",
     lwd = 2, xlab = "False positive rate", ylab = "True positive rate",
     main = "Time-dependent ROC curves (GSE53624)", xlim = c(0, 1), ylim = c(0, 1))
lines(roc_res$FP[, "t=1095"], roc_res$TP[, "t=1095"], col = "#B2182B", lwd = 2)
lines(roc_res$FP[, "t=1825"], roc_res$TP[, "t=1825"], col = "#1B7837", lwd = 2)
abline(0, 1, lty = 2, col = "grey")
legend("bottomright", legend = paste0(c("1-year (AUC=", "3-year (AUC=", "5-year (AUC="),
       aucs, ")"), col = c("#2166AC", "#B2182B", "#1B7837"), lwd = 2, cex = 0.9)
dev.off()

# ---------- 5. Univariable + multivariable Cox ----------
uni <- coxph(surv_obj ~ risk_score, data = risk_df)
s_uni <- summary(uni)
val_lines <- c(val_lines,
  paste0("Univariable Cox (risk score): HR = ", round(s_uni$conf.int[1], 3),
         " (", round(s_uni$conf.int[3], 3), "-", round(s_uni$conf.int[4], 3),
         "), P = ", format.pval(s_uni$coefficients[1, 5], digits = 3)))

cox_df <- data.frame(risk_score = risk_df$risk_score,
                     age = risk_df$age,
                     stage = factor(risk_df$stage,
                                    levels = c("I", "II", "III", "IV")),
                     stringsAsFactors = FALSE)
cc <- complete.cases(cox_df)
cox_df <- cox_df[cc, , drop = FALSE]
cox_df$stage <- droplevels(cox_df$stage)
surv_c <- with(risk_df[cc, ], Surv(os_time, os_status))
multi <- coxph(surv_c ~ risk_score + age + stage, data = cox_df)
s_multi <- summary(multi)
val_lines <- c(val_lines,
  paste0("Multivariable Cox (risk score + age + stage): HR = ",
         round(s_multi$conf.int[1, 1], 3), " (",
         round(s_multi$conf.int[1, 3], 3), "-", round(s_multi$conf.int[1, 4], 3), "), P = ",
         format.pval(s_multi$coefficients[1, 5], digits = 3)))

# ---------- 6. Save ----------
fwrite(risk_df, file.path(out_dir, "risk_score_gse53624.tsv"), sep = "\t")
writeLines(val_lines, file.path(out_dir, "validation_GSE53624.txt"))
cat(paste(val_lines, collapse = "\n"), "\n")
cat("\n09_external_validation done\n")
