## ============================================================
## 17_nested_cv.R
## Purpose: quantify whether a pre-locked candidate pool reduces selection leakage and improves generalizability
##   A) Locked-pool nested CV: candidate pool = 714 consensus DEGs, no univariable pre-screening
##   B) Original-pipeline nested CV: univariable Cox screening + LASSO within each training fold (control)
##   C) Learning curve: training-sample fraction -> test-set C-index (answers "how many samples are needed")
## ============================================================
suppressMessages({ library(survival); library(glmnet); library(matrixStats) })

set.seed(20231021)
out_dir <- "results"

expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
rs   <- read.delim(file.path(out_dir, "risk_score_tcga.tsv"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
deg  <- read.delim(file.path(out_dir, "DEG_consensus.tsv"))

samples <- rs$case_id
expr <- expr[, intersect(samples, colnames(expr))]
clin <- clin[match(samples, clin$case_id), ]
expr <- expr[, samples]
stopifnot(all(clin$case_id == samples))
time <- clin$os_time; status <- clin$os_status
cat("Samples:", length(samples), " events:", sum(status), "\n")

pool <- intersect(deg$gene[deg$is_consensus], rownames(expr))
Xall <- t(expr[pool, , drop = FALSE])[samples, ]
# Filter near-zero-variance genes
v <- apply(Xall, 2, sd, na.rm = TRUE)
Xall <- Xall[, v > 1e-6, drop = FALSE]
cat("Locked candidate pool:", ncol(Xall), "genes\n")

## ---------- C-index (vectorised) ----------
c_index <- function(time, status, lp) {
  ok <- !is.na(lp) & !is.na(time) & !is.na(status)
  t <- time[ok]; s <- status[ok]; r <- lp[ok]
  n <- length(t); conc <- 0; comp <- 0; ties <- 0
  for (i in seq_len(n)) {
    if (s[i] == 0) next
    idx <- which(t > t[i])
    if (!length(idx)) next
    comp <- comp + length(idx)
    conc <- conc + sum(r[idx] < r[i])
    ties <- ties + sum(r[idx] == r[i]) * 0.5
  }
  if (comp == 0) NA else (conc + ties) / comp
}

## ---------- Nested CV engine ----------
nested_cv <- function(X, time, status, n_folds = 10, n_rep = 5,
                      prescreen = FALSE, label = "") {
  res <- data.frame(rep = integer(), fold = integer(), n_test = integer(),
                    events = integer(), c_index = numeric(), n_genes = integer(),
                    genes = character(), stringsAsFactors = FALSE)
  n <- nrow(X)
  for (rep in seq_len(n_rep)) {
    fold_id <- sample(rep_len(seq_len(n_folds), n))
    for (f in seq_len(n_folds)) {
      te <- which(fold_id == f); tr <- which(fold_id != f)
      if (sum(status[te]) < 1 || sum(status[tr]) < 3) next
      Xtr <- X[tr, , drop = FALSE]; Xte <- X[te, , drop = FALSE]
      ytr <- Surv(time[tr], status[tr])

      Xtr_use <- Xtr
      if (prescreen) {                     # original pipeline: univariable Cox screening within the training set
        keep <- vapply(seq_len(ncol(Xtr)), function(k) {
          fit <- try(coxph(ytr ~ Xtr[, k]), silent = TRUE)
          if (inherits(fit, "try-error")) FALSE
          else summary(fit)$coefficients[1, 5] < 0.05
        }, logical(1))
        if (sum(keep) < 3) next
        Xtr_use <- Xtr[, keep, drop = FALSE]
      }

      cvf <- try(cv.glmnet(Xtr_use, ytr, family = "cox", nfolds = 5, alpha = 1,
                           ties = "breslow", maxit = 5000), silent = TRUE)
      if (inherits(cvf, "try-error")) next
      co <- coef(cvf, s = "lambda.min")
      nz <- which(as.numeric(co) != 0)
      genes <- rownames(co)[nz]
      if (length(genes) < 2) next

      lp <- as.numeric(predict(cvf, Xte[, colnames(Xtr_use), drop = FALSE],
                               s = "lambda.min", type = "link"))
      ct <- c_index(time[te], status[te], lp)
      res <- rbind(res, data.frame(rep = rep, fold = f, n_test = length(te),
                                   events = sum(status[te]),
                                   c_index = ct, n_genes = length(genes),
                                   genes = paste(genes, collapse = ";"),
                                   stringsAsFactors = FALSE))
    }
    cat(sprintf("[%s] repeat %d/%d complete, %d folds accumulated\n", label, rep, n_rep, nrow(res)))
  }
  res
}

## ---------- A) Locked pool ----------
cat("\n===== A) Locked pool (714 consensus DEGs, no univariable pre-screening) =====\n")
A <- nested_cv(Xall, time, status, n_folds = 10, n_rep = 5,
               prescreen = FALSE, label = "Locked pool")
A <- A[!is.na(A$c_index), ]
cat(sprintf("Folds: %d | test-set C-index median = %.3f (IQR %.3f-%.3f), mean = %.3f\n",
            nrow(A), median(A$c_index), quantile(A$c_index, .25),
            quantile(A$c_index, .75), mean(A$c_index)))
cat(sprintf("Median number of selected genes: %.0f (range %.0f-%.0f)\n",
            median(A$n_genes), min(A$n_genes), max(A$n_genes)))

# Gene stability (pairwise Jaccard)
gs <- strsplit(A$genes, ";")
jac <- numeric(0)
for (i in 1:(length(gs) - 1)) for (j in (i + 1):length(gs)) {
  jac <- c(jac, length(intersect(gs[[i]], gs[[j]])) / length(union(gs[[i]], gs[[j]])))
}
cat(sprintf("Gene-set stability (mean Jaccard): %.3f\n", mean(jac)))
# Frequent genes
tb <- sort(table(unlist(gs)), decreasing = TRUE)
cat("Top 10 genes by selection frequency:\n"); print(head(tb, 10))

## ---------- B) Original pipeline (control) ----------
cat("\n===== B) Original pipeline (univariable Cox screening + LASSO within each fold) =====\n")
B <- nested_cv(Xall, time, status, n_folds = 10, n_rep = 3,
               prescreen = TRUE, label = "Original pipeline")
B <- B[!is.na(B$c_index), ]
cat(sprintf("Folds: %d | test-set C-index median = %.3f (IQR %.3f-%.3f), mean = %.3f\n",
            nrow(B), median(B$c_index), quantile(B$c_index, .25),
            quantile(B$c_index, .75), mean(B$c_index)))
cat(sprintf("Median number of selected genes: %.0f (range %.0f-%.0f)\n",
            median(B$n_genes), min(B$n_genes), max(B$n_genes)))

## ---------- C) Learning curve ----------
cat("\n===== C) Learning curve (locked pool, training fraction -> generalizability) =====\n")
fracs <- c(0.3, 0.5, 0.7, 0.9)
lc <- data.frame(frac = numeric(), n_train = integer(), events_train = integer(),
                 c_mean = numeric(), c_median = numeric(), c_sd = numeric(),
                 n_rep = integer(), stringsAsFactors = FALSE)
for (fr in fracs) {
  cs <- numeric(0); ntr <- integer(0); ev <- integer(0)
  for (b in 1:20) {
    n_tr <- round(fr * length(samples))
    tr <- sample(seq_along(samples), n_tr)
    te <- setdiff(seq_along(samples), tr)
    if (sum(status[te]) < 2 || sum(status[tr]) < 5) next
    ytr <- Surv(time[tr], status[tr])
    cvf <- try(cv.glmnet(Xall[tr, ], ytr, family = "cox", nfolds = 5, alpha = 1,
                         ties = "breslow", maxit = 5000), silent = TRUE)
    if (inherits(cvf, "try-error")) next
    lp <- as.numeric(predict(cvf, Xall[te, , drop = FALSE], s = "lambda.min", type = "link"))
    ct <- c_index(time[te], status[te], lp)
    if (!is.na(ct)) { cs <- c(cs, ct); ntr <- c(ntr, n_tr); ev <- c(ev, sum(status[tr])) }
  }
  lc <- rbind(lc, data.frame(frac = fr, n_train = round(fr * length(samples)),
                             events_train = round(mean(ev)), c_mean = mean(cs),
                             c_median = median(cs), c_sd = sd(cs), n_rep = length(cs)))
  cat(sprintf("Training fraction %.0f%% (n=%d, events~%d): test C-index mean = %.3f (SD %.3f), median = %.3f [%d reps]\n",
              100 * fr, round(fr * length(samples)), round(mean(ev)), mean(cs), sd(cs), median(cs), length(cs)))
}

## ---------- Save ----------
write.table(A, file.path(out_dir, "nested_cv_lockedpool.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(B, file.path(out_dir, "nested_cv_original.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(lc, file.path(out_dir, "learning_curve.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

suppressMessages(library(ggplot2))
p <- ggplot(lc, aes(x = n_train, y = c_mean)) +
  geom_line(color = "#2E6DA4", linewidth = 1.1) +
  geom_point(size = 3, color = "#B2182B") +
  geom_errorbar(aes(ymin = c_mean - c_sd, ymax = c_mean + c_sd), width = 3, alpha = .5) +
  geom_hline(yintercept = 0.5, linetype = 3, color = "grey40") +
  geom_hline(yintercept = 0.7, linetype = 2, color = "#2E7D46") +
  annotate("text", x = max(lc$n_train), y = 0.71, label = "C-index = 0.7 (clinical-utility reference)",
           hjust = 1, size = 3.6, color = "#2E7D46") +
  labs(x = "Training-set size", y = "Independent test-set C-index (mean ± SD)",
       title = "Learning curve: generalizability versus training-set size") +
  theme_bw(base_size = 12)
ggsave("figures/fig7d_learning_curve.png", p, width = 6.8, height = 5, dpi = 300)

p2 <- ggplot(data.frame(grp = rep(c("Locked pool (714 DEGs)", "Original pipeline (univariable screening)"),
                                 c(nrow(A), nrow(B))),
                        c_index = c(A$c_index, B$c_index)),
             aes(x = grp, y = c_index, fill = grp)) +
  geom_boxplot(width = .6, outlier.size = 1.4) +
  geom_jitter(width = .15, size = 1.5, alpha = .6) +
  geom_hline(yintercept = .5, linetype = 3, color = "grey40") +
  scale_fill_manual(values = c("#2E6DA4", "#B2182B")) +
  labs(x = NULL, y = "Test-set C-index (nested CV)") + theme_bw(base_size = 12) +
  theme(legend.position = "none")
ggsave("figures/fig7e_nested_cv_compare.png", p2, width = 6, height = 5, dpi = 300)

cat("\n===== Summary =====\n")
cat(sprintf("Locked-pool nested CV: C-index median %.3f (n=%d folds)\n", median(A$c_index), nrow(A)))
cat(sprintf("Original-pipeline nested CV: C-index median %.3f (n=%d folds)\n", median(B$c_index), nrow(B)))
cat("Figures saved: fig7d_learning_curve.png, fig7e_nested_cv_compare.png\n")
