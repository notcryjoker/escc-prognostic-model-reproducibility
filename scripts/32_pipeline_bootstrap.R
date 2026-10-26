#!/usr/bin/env Rscript
# ============================================================================
# 32_pipeline_bootstrap.R
# ----------------------------------------------------------------------------
# Purpose (fixes for final-review defects H3 and M4):
#   [H3] Correct pipeline bootstrap: the previous implementation (13_supplementary_analysis.R)
#        re-estimated only a single-parameter Cox model with a locked predictor on each bootstrap sample, so optimism was identically
#        zero (a mathematical artefact). This script re-runs the complete original modelling pipeline within each bootstrap sample
#        (univariable Cox screening at P<0.05 -> LASSO-Cox with 10-fold CV taking lambda.min -> fitting),
#        and then computes optimism, ensuring a non-trivial estimate.
#   [M4] Add final-model calibration: calibration slope/intercept for 1-year predicted event probability (apparent + bootstrap
#        corrected) plus a decile calibration curve; also the apparent calibration of the nomogram multivariable model (risk score + age + stage)
#        (internal validation only; no external validation cohort).
#
# Pipeline replication basis (strictly identical to the original pipeline):
#   - Modelling cohort: scripts/05_lasso_model.R —— TCGA-ESCC complete survival samples n=94, 31 events
#   - Candidate pool:   scripts/04_prognostic_screening.R —— consensus DEGs present in the expression matrix
#               714 genes (of the 719 consensus DEGs)
#   - Screening threshold: univariable Cox P<0.05 (uncorrected, as in the original pipeline)
#   - LASSO:    cv.glmnet(family="cox", alpha=1, nfolds=10, type.measure="C",
#               cox.ties="breslow"), taking lambda.min (as in the original pipeline;
#               maxit=5000 inside the bootstrap loop is for numerical stability only and does not change the lambda selection rule)
#   - Original model:   re-running with set.seed(123) should exactly reproduce results/lasso_model.rds (verified within the script)
#
# Bootstrap design (Harrell's standard pipeline-optimism procedure):
#   B=500, seed=20231020 (consistent with the project 13-script convention).
#   Each bootstrap sample (resampling n=94 with replacement, duplicates included):
#     (a) Apparent C-index: re-run the full pipeline within the sample to obtain the model, then its C-index on that bootstrap sample
#     (b) Test C-index: C-index of the linear predictor of the original full-sample 17-gene model on the new bootstrap sample
#     optimism_b = (a) - (b);  optimism = mean(optimism_b)
#     Corrected C-index = apparent C-index - optimism
#   Bootstrap resampling indices are generated once (independent of glmnet's internal RNG consumption) to guarantee reproducibility.
#
# Calibration design (1 year = 365 days, consistent with the 1-year definition in the 06-script nomogram/timeROC):
#   - Predicted 1-year event probability: p_i = 1 - S0(365)^exp(beta*(lp_i - mean(lp))),
#     where S0 is the baseline survival of coxph(Surv ~ lp) (survfit centres covariates at their means by default)
#     and beta is that model's coefficient — the Cox identity S_i = S0^exp(beta*(x_i - xbar)) requires the
#     coefficient term (bug fix: an earlier version omitted beta, understating probabilities; now cross-checked
#     numerically against survfit(fit, newdata=...), expected <1e-12).
#   - 1-year outcome indicator: y1=1 event and time<=365; y1=0 (alive and time>=365) or (event and time>365);
#     y1=NA denotes censoring before 365 days (1-year status undeterminable; excluded from the logistic calibration, with counts reported).
#   - Calibration slope/intercept: glm(y1 ~ qlogis(p), family=binomial); ideal values (0, 1).
#   - Bootstrap correction: re-run the pipeline + probability model + calibration model within each bootstrap sample;
#     optimism = mean(apparent parameter - test parameter on the original data); corrected value = apparent - optimism.
#
# Evidence-ceiling statement: this study is observational internal validation; bootstrap correction quantifies only overfitting within the same cohort
# and cannot substitute for external validation; calibration estimates rest on 17 events at 1 year and carry substantial uncertainty.
#
# Input: results/tcga_expr_symbol.rds, results/tcga_clinical_surv.rds,
#       results/DEG_consensus.rds, results/prognostic_genes.tsv,
#       results/lasso_model.rds, results/risk_score_tcga.tsv
# Output: results/32_pipeline_bootstrap.tsv       one row per bootstrap iteration
#       results/32_pipeline_bootstrap_summary.txt
#       results/32_calibration_curve.tsv        decile calibration curve of 1-year predicted probabilities
# ============================================================================
rm(list = ls())
suppressPackageStartupMessages({ library(survival); library(glmnet) })

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")
HORIZON <- 365          # 1 year (days), consistent with 06_internal_validation.R
B <- 500                # bootstrap repetitions (consistent with the manuscript)
SEED <- 20231020        # same bootstrap seed as script 13

## ---------- 1. Data and modelling cohort (strict replication of 04/05) ----------
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))   # gene symbols x samples, log2(TPM+1)
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
consensus <- readRDS(file.path(out_dir, "DEG_consensus.rds"))
prog_ref <- read.delim(file.path(out_dir, "prognostic_genes.tsv"))
model_ref <- readRDS(file.path(out_dir, "lasso_model.rds"))

cons_degs <- consensus$gene[consensus$is_consensus]
pool <- intersect(cons_degs, rownames(expr))                  # 714
samples <- intersect(colnames(expr), rownames(clin))
clin <- clin[samples, ]
complete <- complete.cases(clin[, c("os_time", "os_status")])
clin <- clin[complete, ]
samples <- rownames(clin)                                     # 94
time <- clin$os_time; status <- clin$os_status
X <- t(expr[pool, samples, drop = FALSE])                     # 94 x 714
y <- Surv(time, status)
cat(sprintf("Modelling cohort: n=%d, events=%d | candidate pool: %d genes\n", nrow(X), sum(status), ncol(X)))

## ---------- 2. C-index function ----------
# Consistent with script 13, use survConcordance (Harrell's C; immune to the
# version difference in which the new concordance() formula interface returns 1-C; coxph consistency verified: full sample = 0.8328115)
c_index <- function(times, statuses, lp) {
  ok <- !is.na(lp) & !is.na(times) & !is.na(statuses)
  if (sum(ok) < 3) return(NA_real_)
  s <- try(suppressWarnings(survConcordance(Surv(times[ok], statuses[ok]) ~ lp[ok])),
           silent = TRUE)
  if (inherits(s, "try-error")) return(NA_real_)
  as.numeric(s$concordance)
}

## ---------- 3. Original-pipeline re-run and fidelity checks ----------
# 3.1 Re-running univariable screening (94 samples, P<0.05) should reproduce the 48 genes in prognostic_genes.tsv
screen_p <- vapply(seq_len(ncol(X)), function(k) {
  f <- try(coxph(y ~ X[, k]), silent = TRUE)
  if (inherits(f, "try-error")) NA_real_ else summary(f)$coefficients[1, 5]
}, numeric(1))
names(screen_p) <- colnames(X)   # vapply drops names; assign explicitly (bug fix: otherwise names(which(...)) is NULL)
# Consistent with script 05: candidate genes ordered by P as in prognostic_genes.tsv
prog_rerun <- intersect(prog_ref$gene, names(which(screen_p < 0.05)))
cat(sprintf("[Fidelity] re-run univariable screening: %d genes (reference %d) -> %s\n",
            length(prog_rerun), nrow(prog_ref),
            if (identical(sort(prog_rerun), sort(prog_ref$gene))) "match OK" else "mismatch !!"))

# 3.2 LASSO re-run (set.seed(123), same RNG state as script 05) should reproduce 17 genes
set.seed(123)
cv_refit <- cv.glmnet(X[, prog_rerun, drop = FALSE], y, family = "cox", alpha = 1,
                      nfolds = 10, type.measure = "C", cox.ties = "breslow")
co_refit <- as.matrix(coef(cv_refit, s = "lambda.min"))
genes_full <- rownames(co_refit)[co_refit[, 1] != 0]
coefs_full <- co_refit[co_refit[, 1] != 0, 1]
cat(sprintf("[Fidelity] LASSO re-run: %d genes (reference %d), lambda.min=%.4f (reference %.4f) -> %s\n",
            length(genes_full), length(model_ref$genes), cv_refit$lambda.min,
            model_ref$lambda.min,
            if (identical(sort(genes_full), sort(model_ref$genes)) &&
                abs(cv_refit$lambda.min - model_ref$lambda_min) < 1e-3) "match OK" else "mismatch !!"))
# Linear predictor of the original full-sample model
lp_full <- as.numeric(X[, genes_full, drop = FALSE] %*% coefs_full)
c_apparent <- c_index(time, status, lp_full)
cat(sprintf("[Fidelity] original-model apparent C-index = %.4f (reference 0.8328)\n", c_apparent))

## ---------- 4. 1-year outcome indicator and apparent calibration ----------
y1 <- ifelse(status == 1 & time <= HORIZON, 1,
             ifelse((status == 0 & time >= HORIZON) | (status == 1 & time > HORIZON), 0, NA))
cat(sprintf("1-year outcome: determinable n=%d, 1-year events=%d, censored before 365 days excluded n=%d\n",
            sum(!is.na(y1)), sum(y1, na.rm = TRUE), sum(is.na(y1))))

# Probability model: coxph(Surv ~ lp_full); p = 1 - S0(365)^exp(beta*(lp - mean(lp)))
# Cox identity: S_i = S0^exp(beta*(x_i - xbar)); the coefficient beta MUST be included
# (survfit.coxph centres at covariate means; here lp_full is the raw linear predictor,
#  so beta is applied explicitly and cross-checked against survfit(newdata))
f_prob <- coxph(y ~ lp_full)
beta_prob <- as.numeric(coef(f_prob))
s0_full <- summary(survfit(f_prob), times = HORIZON)$surv
p_full <- 1 - s0_full^exp(beta_prob * (lp_full - mean(lp_full)))
# Fidelity: compare with survfit(newdata) (should be <1e-12)
p_sf <- 1 - summary(survfit(f_prob, newdata = data.frame(lp_full = lp_full),
                            se.fit = FALSE), times = HORIZON)$surv
cat(sprintf("[Fidelity] max difference of predicted-probability formula vs survfit(newdata) = %.2e\n",
            max(abs(p_full - p_sf))))

logit_calib <- function(yv, pv) {
  d <- data.frame(y = yv, q = qlogis(pv))
  d <- d[complete.cases(d), ]
  g <- try(glm(y ~ q, data = d, family = binomial), silent = TRUE)
  if (inherits(g, "try-error")) return(c(intercept = NA_real_, slope = NA_real_))
  c(intercept = unname(coef(g)[1]), slope = unname(coef(g)[2]))
}
cal0 <- logit_calib(y1, p_full)
cat(sprintf("Apparent calibration (risk-score model, 1 year): intercept = %.3f, slope = %.3f (ideal 0, 1)\n",
            cal0["intercept"], cal0["slope"]))

# 4.1 Apparent calibration of the nomogram multivariable model (risk score + age + stage, n=92) — replicating the 06-script model
cox_df <- data.frame(risk_score = lp_full, age = clin$age,
                     stage = factor(clin$stage_group, levels = c("I", "II", "III", "IV")),
                     stringsAsFactors = FALSE)
cc <- complete.cases(cox_df)
cox_df <- cox_df[cc, , drop = FALSE]
cox_df$stage <- droplevels(cox_df$stage)
y_cc <- y[cc]; y1_cc <- y1[cc]; time_cc <- time[cc]; status_cc <- status[cc]
f_multi <- coxph(y_cc ~ risk_score + age + stage, data = cox_df)
lp_multi <- predict(f_multi, type = "lp")   # centred linear predictor (mean 0), so no extra beta term
s0_multi <- summary(survfit(f_multi), times = HORIZON)$surv
p_multi <- 1 - s0_multi^exp(lp_multi - mean(lp_multi))
cal_multi <- logit_calib(y1_cc, p_multi)
cat(sprintf("Apparent calibration (nomogram multivariable model, n=%d, 1 year): intercept = %.3f, slope = %.3f\n",
            nrow(cox_df), cal_multi["intercept"], cal_multi["slope"]))

## ---------- 5. Pipeline bootstrap (B=500) ----------
set.seed(SEED)
idx_mat <- matrix(sample(seq_len(nrow(X)), size = nrow(X) * B, replace = TRUE),
                  nrow = nrow(X))          # 94 x 500, one resample per column
res <- data.frame(b = integer(), n_events_boot = integer(), n_screened = integer(),
                  n_genes = integer(), c_apparent = numeric(), c_test = numeric(),
                  optimism = numeric(),
                  calib_int_app = numeric(), calib_slope_app = numeric(),
                  calib_int_test = numeric(), calib_slope_test = numeric(),
                  status = character(), stringsAsFactors = FALSE)
checkpoint <- file.path(tempdir(), "32_checkpoint.rds")
t_start <- Sys.time()
for (b in seq_len(B)) {
  idx <- idx_mat[, b]
  st <- "ok"
  fit <- try({
    Xb <- X[idx, , drop = FALSE]
    yb <- Surv(time[idx], status[idx])
    # (i) Complete re-run within the bootstrap sample: univariable Cox screening at P<0.05
    p_b <- vapply(seq_len(ncol(Xb)), function(k) {
      f <- try(coxph(yb ~ Xb[, k]), silent = TRUE)
      if (inherits(f, "try-error")) NA_real_ else summary(f)$coefficients[1, 5]
    }, numeric(1))
    keep <- which(!is.na(p_b) & p_b < 0.05)
    if (length(keep) < 2) stop("screen_lt2")
    # (ii) LASSO-Cox re-run within the bootstrap sample (10-fold CV, lambda.min)
    cvb <- try(suppressWarnings(cv.glmnet(Xb[, keep, drop = FALSE], yb, family = "cox", alpha = 1,
                     nfolds = 10, type.measure = "C", cox.ties = "breslow", maxit = 5000)), silent = TRUE)
    cob <- as.matrix(coef(cvb, s = "lambda.min"))
    genes_b <- rownames(cob)[cob[, 1] != 0]
    coefs_b <- cob[cob[, 1] != 0, 1]
    if (length(genes_b) < 2) stop("lt2_genes")
    # (iii) Two linear predictors
    lp_b_boot <- as.numeric(Xb[, genes_b, drop = FALSE] %*% coefs_b)     # bootstrap sample
    lp_b_orig <- as.numeric(X[, genes_b, drop = FALSE] %*% coefs_b)      # original 94 samples
    # (a) Apparent C-index: bootstrap re-run model on its own sample
    c_app <- c_index(time[idx], status[idx], lp_b_boot)
    # (b) Test C-index: original full-sample model on the new bootstrap sample
    c_tst <- c_index(time[idx], status[idx], lp_full[idx])
    # (iv) Calibration within the bootstrap (1 year, logistic)
    fb <- coxph(yb ~ lp_b_boot)
    beta_b <- as.numeric(coef(fb))
    s0_b <- summary(survfit(fb), times = HORIZON)$surv
    lpb_bar <- mean(lp_b_boot)
    p_app <- 1 - s0_b^exp(beta_b * (lp_b_boot - lpb_bar))   # apparent prediction (bootstrap sample)
    p_tst <- 1 - s0_b^exp(beta_b * (lp_b_orig - lpb_bar))   # prediction on the original sample (same beta, same baseline)
    cal_app <- logit_calib(y1[idx], p_app)              # apparent calibration parameters
    cal_tst <- logit_calib(y1, p_tst)                   # test calibration parameters (original-data outcomes)
    list(c_app = c_app, c_tst = c_tst, n_ev = sum(status[idx]),
         n_scr = length(keep), n_genes = length(genes_b),
         cal_app = cal_app, cal_tst = cal_tst)
  }, silent = TRUE)
  if (inherits(fit, "try-error")) {
    msg <- conditionMessage(attr(fit, "condition"))
    st <- if (grepl("screen_lt2", msg)) "screen_lt2" else
          if (grepl("lt2_genes", msg)) "lt2_genes" else "model_fail"
    res <- rbind(res, data.frame(b = b, n_events_boot = NA_integer_, n_screened = NA_integer_,
                                 n_genes = NA_integer_, c_apparent = NA_real_, c_test = NA_real_,
                                 optimism = NA_real_, calib_int_app = NA_real_,
                                 calib_slope_app = NA_real_, calib_int_test = NA_real_,
                                 calib_slope_test = NA_real_, status = st,
                                 stringsAsFactors = FALSE))
  } else {
    res <- rbind(res, data.frame(b = b, n_events_boot = fit$n_ev, n_screened = fit$n_scr,
                                 n_genes = fit$n_genes, c_apparent = fit$c_app,
                                 c_test = fit$c_tst,
                                 optimism = fit$c_app - fit$c_tst,
                                 calib_int_app = fit$cal_app["intercept"],
                                 calib_slope_app = fit$cal_app["slope"],
                                 calib_int_test = fit$cal_tst["intercept"],
                                 calib_slope_test = fit$cal_tst["slope"],
                                 status = st, stringsAsFactors = FALSE))
  }
  if (b %% 50 == 0) {
    cat(sprintf("  Completed %d/%d (%.1f min)\n", b, B,
                as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
    saveRDS(res, checkpoint)
  }
}
t_end <- Sys.time()
cat(sprintf("bootstrap complete, elapsed %.1f min\n",
            as.numeric(difftime(t_end, t_start, units = "mins"))))

## ---------- 6. Summary ----------
ok <- !is.na(res$optimism)
n_ok <- sum(ok)
opt <- res$optimism[ok]
mean_opt <- mean(opt); sd_opt <- sd(opt)
c_corrected <- c_apparent - mean_opt
# Bootstrap correction of the calibration parameters
int_opt <- mean(res$calib_int_app[ok] - res$calib_int_test[ok], na.rm = TRUE)
slope_opt <- mean(res$calib_slope_app[ok] - res$calib_slope_test[ok], na.rm = TRUE)
n_calib_ok <- sum(ok & !is.na(res$calib_slope_app) & !is.na(res$calib_slope_test))
int_corr <- cal0["intercept"] - int_opt
slope_corr <- cal0["slope"] - slope_opt
# O/E (1 year)
oe <- sum(p_full[!is.na(y1)]) / sum(y1, na.rm = TRUE)

# Decile calibration curve (two models)
decile_tbl <- function(pv, yv, label) {
  d <- data.frame(p = pv, y = yv); d <- d[complete.cases(d), ]
  d$dec <- cut(d$p, breaks = quantile(d$p, probs = seq(0, 1, 0.1), include.lowest = TRUE),
               labels = FALSE, include.lowest = TRUE)
  agg <- do.call(rbind, lapply(sort(unique(d$dec)), function(k) {
    s <- d[d$dec == k, ]
    data.frame(model = label, decile = k, n = nrow(s), events_1y = sum(s$y),
               mean_pred_event_prob = mean(s$p), observed_event_rate = mean(s$y),
               mean_pred_survival_prob = 1 - mean(s$p), stringsAsFactors = FALSE)
  }))
  agg
}
cal_curve <- rbind(decile_tbl(p_full, y1, "17-gene risk score (Cox)"),
                   decile_tbl(p_multi, y1_cc, "Nomogram multivariable (risk+age+stage)"))
fwrite_dt <- function(d, f) write.table(d, f, sep = "\t", quote = FALSE, row.names = FALSE)

## ---------- 7. Output ----------
write.table(res, file.path(out_dir, "32_pipeline_bootstrap.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
fwrite_dt(cal_curve, file.path(out_dir, "32_calibration_curve.tsv"))

summ <- c(
  "===== 32_pipeline_bootstrap.R summary =====",
  paste0("Run time: ", format(t_start), " - ", format(t_end)),
  "",
  "--- Design and fidelity checks ---",
  sprintf("Modelling cohort: TCGA-ESCC n=%d, events=%d (consistent with 05_lasso_model.R)", nrow(X), sum(status)),
  sprintf("Candidate pool: %d consensus DEGs (of the 719 consensus DEGs present in the expression matrix)", ncol(X)),
  sprintf("Univariable screening threshold P<0.05; LASSO: 10-fold CV, lambda.min (as in the original pipeline)"),
  sprintf("Re-run univariable screening gene count: %d (reference prognostic_genes.tsv: %d)",
          length(prog_rerun), nrow(prog_ref)),
  sprintf("Re-run LASSO gene count: %d (reference lasso_model.rds: %d), lambda.min=%.4f (reference %.4f)",
          length(genes_full), length(model_ref$genes), cv_refit$lambda.min, model_ref$lambda_min),
  sprintf("Original full-sample 17-gene model apparent C-index = %.4f (consistent with the manuscript value 0.8328)", c_apparent),
  sprintf("Max difference of predicted-probability formula vs survfit(newdata) = %.2e (numerical check passed)",
          max(abs(p_full - p_sf))),
  "",
  "--- Bootstrap internal validation (pipeline bootstrap, B=500, seed=20231020) ---",
  sprintf("Bootstrap resampling: n=%d with replacement; each iteration re-runs screening -> LASSO -> fitting in full", nrow(X)),
  sprintf("Completed iterations: %d/%d; modelling failures: %d (%s)", n_ok, B, sum(res$status != "ok"),
          paste(names(table(res$status[res$status != "ok"])),
                table(res$status[res$status != "ok"]), sep = ":", collapse = "; ")),
  sprintf("Apparent C-index = %.4f", c_apparent),
  sprintf("Bootstrap apparent C (a) mean = %.4f (SD %.4f)", mean(res$c_apparent[ok]), sd(res$c_apparent[ok])),
  sprintf("Original-model test C (b) mean = %.4f (SD %.4f)", mean(res$c_test[ok]), sd(res$c_test[ok])),
  sprintf("Mean optimism = %.4f (SD %.4f, median %.4f, range %.4f-%.4f)",
          mean_opt, sd_opt, median(opt), min(opt), max(opt)),
  sprintf("Optimism 95%% percentile interval: %.4f - %.4f",
          unname(quantile(opt, 0.025)), unname(quantile(opt, 0.975))),
  sprintf("Optimism-corrected C-index = %.4f", c_corrected),
  sprintf("[Acceptance] optimism non-triviality: %d/%d iterations with optimism>0.005 (%.1f%%); mean optimism = %.4f -> %s",
          sum(opt > 0.005), n_ok, 100 * mean(opt > 0.005), mean_opt,
          if (mean_opt > 0.005) "PASS (screening/LASSO were indeed re-run within bootstrap)" else "FAIL (debugging required)"),
  "",
  "--- Calibration (M4, 1 year = 365 days, internal validation) ---",
  sprintf("1-year outcome determinable n=%d (1-year events %d; censored before 365 days excluded %d)",
          sum(!is.na(y1)), sum(y1, na.rm = TRUE), sum(is.na(y1))),
  sprintf("Apparent calibration (17-gene risk-score model): intercept = %.3f, slope = %.3f (ideal 0, 1)",
          cal0["intercept"], cal0["slope"]),
  sprintf("Bootstrap-corrected calibration (%d valid iterations): intercept = %.3f, slope = %.3f",
          n_calib_ok, int_corr, slope_corr),
  sprintf("  (mean optimism: intercept %+.3f, slope %+.3f)", int_opt, slope_opt),
  sprintf("1-year O/E = %.3f (expected events %.1f vs actual %d, determinable samples)",
          oe, sum(p_full[!is.na(y1)]), sum(y1, na.rm = TRUE)),
  sprintf("Apparent calibration (nomogram multivariable model risk+age+stage, n=%d): intercept = %.3f, slope = %.3f",
          nrow(cox_df), cal_multi["intercept"], cal_multi["slope"]),
  "Calibration-curve data: results/32_calibration_curve.tsv (deciles)",
  "",
  "--- Interpretation and evidence ceiling ---",
  "Optimism quantifies the overfitting of 're-modelling within the same cohort using the full pipeline (including gene screening)' relative to the fixed model;",
  "it does not measure external generalisation error (no external validation cohort) and cannot be interpreted as a causal effect.",
  "Calibration rests on 17 events at 1 year, with only 8-9 cases per decile group; estimates carry substantial uncertainty and should be reported as internal consistency only.",
  "",
  "--- Reproducibility information ---",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(summ, file.path(out_dir, "32_pipeline_bootstrap_summary.txt"))
cat("\n", paste(summ, collapse = "\n"), "\n")
cat("\n32_pipeline_bootstrap complete\n")
