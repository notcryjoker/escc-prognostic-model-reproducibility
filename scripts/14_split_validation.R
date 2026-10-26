## ============================================================
## Internal robustness validation: 100 random 7:3 splits (model in the training set, evaluate in the test set)
## Purpose: in the absence of an external cohort, use repeated split validation to assess the reproducibility of the modelling pipeline
##       and the predictive performance of the model on samples not used in modelling
## ============================================================
suppressMessages({
  library(survival); library(glmnet); library(matrixStats)
})

set.seed(20231021)
out_dir <- "results"

expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))     # genes x samples
rs   <- read.delim(file.path(out_dir, "risk_score_tcga.tsv"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))

common <- intersect(rs$case_id, colnames(expr))
expr <- expr[, common]
clin <- clin[match(common, clin$case_id), ]
stopifnot(all(clin$case_id == common))

surv_df <- data.frame(rownames = common,
                      time = clin$os_time, status = clin$os_status)
cc <- complete.cases(surv_df$time, surv_df$status)
surv_df <- surv_df[cc, ]
expr <- expr[, surv_df$rownames]
cat("Analysis samples:", nrow(surv_df), " events:", sum(surv_df$status), "\n")

# Candidate gene pool: the 48 genes with univariable Cox P<0.05 (reused from the main analysis)
pg <- read.delim(file.path(out_dir, "prognostic_genes.tsv"))
cand <- intersect(pg$gene, rownames(expr))
cat("Candidate gene pool:", length(cand), "genes\n")

X <- t(expr[cand, , drop = FALSE])      # samples x genes
X <- X[surv_df$rownames, ]

cidx <- function(times, status, lp) {
  # Manually compute the C-index (Harrell's C)
  ok <- !is.na(lp) & !is.na(times) & !is.na(status)
  t <- times[ok]; s <- status[ok]; r <- lp[ok]
  n <- length(t)
  comp <- 0; concordant <- 0; ties <- 0
  for (i in seq_len(n)) {
    if (s[i] == 0) next
    for (j in seq_len(n)) {
      if (t[j] <= t[i]) next
      comp <- comp + 1
      if (r[j] < r[i]) concordant <- concordant + 1
      else if (r[j] == r[i]) ties <- ties + 0.5
    }
  }
  if (comp == 0) return(NA)
  (concordant + ties) / comp
}

B <- 100
res <- data.frame(rep = integer(), n_train = integer(), n_test = integer(),
                  n_events_test = integer(), n_genes = integer(),
                  c_test = numeric(), p_test = numeric())

for (b in seq_len(B)) {
  idx_tr <- sample(seq_len(nrow(X)), size = floor(0.7 * nrow(X)))
  idx_te <- setdiff(seq_len(nrow(X)), idx_tr)

  Xtr <- X[idx_tr, ]; Xte <- X[idx_te, ]
  ytr <- Surv(surv_df$time[idx_tr], surv_df$status[idx_tr])

  # Training set: univariable Cox screening (P<0.05)
  keep <- sapply(seq_len(ncol(Xtr)), function(k) {
    f <- try(coxph(ytr ~ Xtr[, k]), silent = TRUE)
    if (inherits(f, "try-error")) return(FALSE)
    summary(f)$coefficients[1, 5] < 0.05
  })
  if (sum(keep) < 3) next
  Xtr2 <- Xtr[, keep, drop = FALSE]

  # LASSO-Cox (10-fold CV)
  cvf <- try(cv.glmnet(Xtr2, ytr, family = "cox", nfolds = 10,
                       alpha = 1, ties = "breslow", maxit = 5000), silent = TRUE)
  if (inherits(cvf, "try-error")) next
  co <- coef(cvf, s = "lambda.min")
  nz <- which(as.numeric(co) != 0) - 1
  genes_sel <- rownames(co)[which(as.numeric(co) != 0)]
  if (length(genes_sel) < 2) next

  # Test-set risk score
  lp_te <- as.numeric(predict(cvf, Xte[, colnames(Xtr2), drop = FALSE],
                              s = "lambda.min", type = "link"))
  ct <- cidx(surv_df$time[idx_te], surv_df$status[idx_te], lp_te)

  # Test-set KM test (grouped at the training-set risk-score median cut-off)
  cut <- median(as.numeric(predict(cvf, Xtr2, s = "lambda.min", type = "link")))
  grp <- ifelse(lp_te > cut, "High", "Low")
  sd_ <- try(survdiff(Surv(surv_df$time[idx_te], surv_df$status[idx_te]) ~ grp),
             silent = TRUE)
  # survdiff objects have no $p.value; derive it from $chisq
  p_km <- if (inherits(sd_, "try-error")) NA_real_ else {
    ch <- suppressWarnings(as.numeric(sd_$chisq))
    if (length(ch) != 1 || is.na(ch)) NA_real_ else pchisq(ch, df = 1, lower.tail = FALSE)
  }
  if (length(p_km) != 1 || is.na(p_km)) p_km <- NA_real_

  res <- rbind(res, data.frame(rep = b, n_train = length(idx_tr), n_test = length(idx_te),
                               n_events_test = sum(surv_df$status[idx_te]),
                               n_genes = length(genes_sel),
                               c_test = ct, p_test = p_km))
  if (b %% 20 == 0) cat("Done", b, "/", B, "\n")
}

# Filter failed iterations
res <- res[!is.na(res$c_test), ]
cat("\n===== Split validation summary (", nrow(res), " successful runs) =====\n")
cat(sprintf("Test-set C-index: median = %.3f (IQR %.3f-%.3f), mean = %.3f, range %.3f-%.3f\n",
            median(res$c_test), quantile(res$c_test, 0.25), quantile(res$c_test, 0.75),
            mean(res$c_test), min(res$c_test), max(res$c_test)))
cat(sprintf("Proportion of test-set C-index > 0.5: %.1f%%\n", 100 * mean(res$c_test > 0.5)))
cat(sprintf("Proportion of test-set C-index > 0.7: %.1f%%\n", 100 * mean(res$c_test > 0.7)))
cat(sprintf("Proportion of test-set log-rank P < 0.05: %.1f%%\n",
            100 * mean(res$p_test < 0.05, na.rm = TRUE)))
cat(sprintf("Median number of genes selected per model: %.0f (range %.0f-%.0f)\n",
            median(res$n_genes), min(res$n_genes), max(res$n_genes)))

write.table(res, file.path(out_dir, "split_validation.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Figure: C-index distribution
suppressMessages(library(ggplot2))
p <- ggplot(res, aes(x = c_test)) +
  geom_histogram(bins = 18, fill = "#2E6DA4", color = "white", alpha = 0.9) +
  geom_vline(xintercept = median(res$c_test), color = "#B2182B", linetype = 2, linewidth = 1) +
  geom_vline(xintercept = 0.5, color = "grey40", linetype = 3) +
  annotate("text", x = median(res$c_test), y = Inf, vjust = 1.6, hjust = -0.1,
           label = paste0("Median ", round(median(res$c_test), 3)), size = 4, color = "#B2182B") +
  labs(x = "Test-set C-index", y = "Frequency",
       title = paste0("100 random 7:3 split validation (n = ", nrow(res), " runs)")) +
  theme_bw(base_size = 12)
ggsave("figures/fig7b_split_validation.png", p, width = 6.5, height = 4.8, dpi = 300)
cat("Figure saved: figures/fig7b_split_validation.png\n")
