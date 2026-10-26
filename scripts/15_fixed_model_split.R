## ============================================================
## Performance of the fixed 17-gene model (TCGA full-data coefficients) across 100 7:3 splits
## on the test set — used to distinguish
##   (A) the generalizability of the model itself (this script; optimism from gene selection remains)
##   (B) the reproducibility of the full modelling pipeline (script 14; unbiased estimate)
## ============================================================
suppressMessages({ library(survival); library(ggplot2) })

set.seed(20231021)
out_dir <- "results"

expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))
rs   <- read.delim(file.path(out_dir, "risk_score_tcga.tsv"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
mdl  <- readRDS(file.path(out_dir, "lasso_model.rds"))

common <- intersect(rs$case_id, colnames(expr))
expr <- expr[, common]
clin <- clin[match(common, clin$case_id), ]
df <- data.frame(id = common, time = clin$os_time, status = clin$os_status)
df <- df[complete.cases(df$time, df$status), ]
expr <- expr[, df$id]

genes <- mdl$genes
coefs <- as.numeric(mdl$coef)
names(coefs) <- genes
gen_avail <- intersect(genes, rownames(expr))
cat("Available model genes:", length(gen_avail), "/", length(genes), "\n")

# Compute whole-cohort risk scores (fixed coefficients)
X <- t(expr[gen_avail, , drop = FALSE])[df$id, ]
risk <- as.numeric(as.matrix(X) %*% coefs[gen_avail])
cutoff <- median(risk)
cat("Whole-cohort risk-score median =", round(cutoff, 4), "\n")

cidx <- function(times, status, lp) {
  ok <- !is.na(lp); t <- times[ok]; s <- status[ok]; r <- lp[ok]
  n <- length(t); comp <- 0; conc <- 0; ties <- 0
  for (i in seq_len(n)) {
    if (s[i] == 0) next
    for (j in seq_len(n)) {
      if (t[j] <= t[i]) next
      comp <- comp + 1
      if (r[j] < r[i]) conc <- conc + 1 else if (r[j] == r[i]) ties <- ties + 0.5
    }
  }
  if (comp == 0) NA else (conc + ties) / comp
}

B <- 100
res <- data.frame(rep = integer(), n_test = integer(), events = integer(),
                  c_test = numeric(), p_km = numeric())
for (b in seq_len(B)) {
  ite <- sample(seq_len(nrow(df)), size = round(0.3 * nrow(df)))
  lp_te <- risk[ite]
  ct <- cidx(df$time[ite], df$status[ite], lp_te)
  grp <- ifelse(lp_te > cutoff, "High", "Low")
  sd_ <- try(survdiff(Surv(df$time[ite], df$status[ite]) ~ grp), silent = TRUE)
  p <- if (inherits(sd_, "try-error")) NA_real_ else {
    ch <- suppressWarnings(as.numeric(sd_$chisq))
    if (length(ch) != 1 || is.na(ch)) NA_real_ else pchisq(ch, df = 1, lower.tail = FALSE)
  }
  res <- rbind(res, data.frame(rep = b, n_test = length(ite),
                               events = sum(df$status[ite]), c_test = ct, p_km = p))
}
res <- res[!is.na(res$c_test), ]
cat("\n===== Fixed-model split validation (", nrow(res), " runs) =====\n")
cat(sprintf("Test-set C-index: median = %.3f (IQR %.3f-%.3f), mean = %.3f, range %.3f-%.3f\n",
            median(res$c_test), quantile(res$c_test, .25), quantile(res$c_test, .75),
            mean(res$c_test), min(res$c_test), max(res$c_test)))
cat(sprintf("C-index > 0.5: %.1f%% | > 0.6: %.1f%% | > 0.7: %.1f%%\n",
            100 * mean(res$c_test > .5), 100 * mean(res$c_test > .6),
            100 * mean(res$c_test > .7)))
cat(sprintf("log-rank P < 0.05: %.1f%%\n", 100 * mean(res$p_km < .05, na.rm = TRUE)))

write.table(res, file.path(out_dir, "split_validation_fixedmodel.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot(res, aes(x = c_test)) +
  geom_histogram(bins = 18, fill = "#0C7C94", color = "white", alpha = .9) +
  geom_vline(xintercept = median(res$c_test), color = "#B2182B", linetype = 2, linewidth = 1) +
  geom_vline(xintercept = .5, color = "grey40", linetype = 3) +
  annotate("text", x = median(res$c_test), y = Inf, vjust = 1.6, hjust = -.1,
           label = paste0("Median ", round(median(res$c_test), 3)), size = 4, color = "#B2182B") +
  labs(x = "Test-set C-index (fixed 17-gene model)", y = "Frequency") +
  theme_bw(base_size = 12)
ggsave("figures/fig7c_split_fixedmodel.png", p, width = 6.5, height = 4.8, dpi = 300)
cat("Figure saved\n")
