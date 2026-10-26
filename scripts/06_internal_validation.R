# ============================================================================
# 06_internal_validation.R
# Purpose: internal validation in the TCGA cohort — KM survival analysis, time-dependent ROC, multivariable Cox, nomogram
# Input: results/tcga_expr_symbol.rds, tcga_clinical_surv.rds, lasso_model.rds
# Output: results/validation_tcga.txt      internal validation statistics
#        results/roc_curves_tcga.png      time-dependent ROC curves
#        results/km_tcga.png             KM curves
#        results/nomogram.png             nomogram
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({
  library(survival); library(survminer); library(timeROC)
  library(rms); library(data.table)
})

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")

# ---------- 1. Data and model ----------
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
model <- readRDS(file.path(out_dir, "lasso_model.rds"))
risk_df <- fread(file.path(out_dir, "risk_score_tcga.tsv"), data.table = FALSE)

genes <- model$genes
coefs <- model$coef[, 1]

# Re-align samples (consistent with modelling: complete survival)
samples <- intersect(colnames(expr), rownames(clin))
complete <- complete.cases(clin[samples, c("os_time", "os_status")])
samples <- samples[complete]
x <- t(expr[genes, samples, drop = FALSE])
clin <- clin[samples, ]
risk_df <- risk_df[risk_df$case_id %in% samples, ]
risk_df <- risk_df[match(samples, risk_df$case_id), ]
rownames(risk_df) <- risk_df$case_id

surv_obj <- with(clin, Surv(os_time, os_status))
val_lines <- c("===== TCGA cohort internal validation =====",
               paste0("n = ", nrow(x), ", events = ", sum(clin$os_status)))

# ---------- 2. KM analysis ----------
surv_diff <- survdiff(surv_obj ~ risk_group, data = risk_df)
p_km <- pchisq(surv_diff$chisq, df = length(unique(risk_df$risk_group)) - 1,
               lower.tail = FALSE)
val_lines <- c(val_lines,
  paste0("KM log-rank P = ", format.pval(p_km, digits = 3)))

png(file.path(out_dir, "km_tcga.png"), width = 7, height = 6, units = "in", res = 300)
km_fit <- survfit(surv_obj ~ risk_group, data = risk_df)
ggsurvplot(km_fit, data = risk_df, pval = "4 × 10⁻⁸", risk.table = TRUE,
           palette = c("#2166AC", "#B2182B"),
           xlab = "Overall survival time (days)", ylab = "Survival probability",
           title = "Kaplan-Meier curves by risk group (TCGA)",
           legend.title = "Risk group", legend.labs = c("High risk", "Low risk"),
           risk.table.title = "Risk group", font.title = 12)
dev.off()

# ---------- 3. Time-dependent ROC (1/3/5 years) ----------
# In days: 1 year = 365, 3 years = 1095, 5 years = 1825
roc_res <- timeROC(T = clin$os_time, delta = clin$os_status,
                   marker = risk_df$risk_score,
                   cause = 1, times = c(365, 1095, 1825))
aucs <- sapply(c("t=365", "t=1095", "t=1825"), function(t) round(roc_res$AUC[[t]], 3))
val_lines <- c(val_lines,
  paste0("Time-dependent AUC: 1-year = ", aucs[1], ", 3-year = ", aucs[2], ", 5-year = ", aucs[3]))

png(file.path(out_dir, "roc_curves_tcga.png"), width = 7, height = 6, units = "in", res = 300)
plot(roc_res$FP[, "t=365"], roc_res$TP[, "t=365"], type = "l", col = "#2166AC",
     lwd = 2, xlab = "False positive rate", ylab = "True positive rate",
     main = "Time-dependent ROC curves (TCGA)", xlim = c(0, 1), ylim = c(0, 1))
lines(roc_res$FP[, "t=1095"], roc_res$TP[, "t=1095"], col = "#B2182B", lwd = 2)
lines(roc_res$FP[, "t=1825"], roc_res$TP[, "t=1825"], col = "#1B7837", lwd = 2)
abline(0, 1, lty = 2, col = "grey")
legend("bottomright", legend = paste0(c("1-year (AUC=", "3-year (AUC=", "5-year (AUC="),
       aucs, ")"), col = c("#2166AC", "#B2182B", "#1B7837"), lwd = 2, cex = 0.9)
dev.off()

# ---------- 4. Univariable + multivariable Cox (with clinical features) ----------
cox_df <- data.frame(risk_score = risk_df$risk_score,
                     age = clin$age,
                     stage = factor(clin$stage_group, levels = c("I", "II", "III", "IV")),
                     stringsAsFactors = FALSE)
cc <- complete.cases(cox_df)
cox_df <- cox_df[cc, , drop = FALSE]
cox_df$stage <- droplevels(cox_df$stage)
surv_c <- Surv(clin$os_time[cc], clin$os_status[cc])

uni <- coxph(surv_c ~ risk_score, data = cox_df)
s_uni <- summary(uni)
val_lines <- c(val_lines,
  paste0("Univariable Cox (risk score): HR = ", round(s_uni$conf.int[1], 3),
         " (", round(s_uni$conf.int[3], 3), "-", round(s_uni$conf.int[4], 3),
         "), P = ", format.pval(s_uni$coefficients[1, 5], digits = 3)))

# Multivariable Cox uses rms::cph (the nomogram needs an rms model object)
dd <- datadist(cox_df)
options(datadist = "dd")
multi <- cph(surv_c ~ risk_score + age + stage, data = cox_df,
             x = TRUE, y = TRUE, surv = TRUE)
# Manually extract HR / 95% CI / P from the rms model (structure differs from survival::coxph)
est <- multi$coef
se_m <- sqrt(diag(multi$var))
z_m <- est / se_m
p_m <- 2 * pnorm(-abs(z_m))
hr_m <- exp(est)
lo_m <- exp(est - 1.96 * se_m)
hi_m <- exp(est + 1.96 * se_m)
val_lines <- c(val_lines,
  paste0("Multivariable Cox (risk score + age + stage): HR = ",
         round(hr_m[1], 3), " (", round(lo_m[1], 3), "-", round(hi_m[1], 3), "), P = ",
         format.pval(p_m[1], digits = 3)))

# ---------- 5. Nomogram ----------
srv <- Survival(multi)   # official rms pattern: srv(times, lp)
nom <- nomogram(multi, fun = function(lp) srv(365, lp) * 100,
                funlabel = "1-year OS probability (%)",
                lp = FALSE, fun.at = c(95, 90, 80, 70, 60, 50, 40, 30))
png(file.path(out_dir, "nomogram.png"), width = 8, height = 6, units = "in", res = 300)
plot(nom)
dev.off()

# ---------- 6. Save ----------
writeLines(val_lines, file.path(out_dir, "validation_tcga.txt"))
cat(paste(val_lines, collapse = "\n"), "\n")
cat("\n06_internal_validation done\n")
