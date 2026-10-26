#!/usr/bin/env Rscript
# ============================================================================
# 33_control_split.R
# ----------------------------------------------------------------------------
# 目的（对应终审缺陷 H7 的修复）:
#   将"基因集泄漏"与"训练样本量/系数泄漏效应"分离。论文核心主张"选择泄漏约
#   0.29 C-index 单位"（0.850 固定策略 vs 0.556 重跑策略）的落差同时包含:
#     (i) 基因集用全部 94 样本（含测试样本）选定 -> 基因集选择泄漏;
#     (ii) 固定策略系数亦用全部 94 样本估计、重跑策略仅用 65 个训练样本
#          -> 系数泄漏 + 训练样本量效应（94 vs 65）。
#   本脚本增加对照拆分（基因集固定为全样本 17 基因，系数仅在训练集内重估）,
#   与两个既有策略构成同拆分下的三向分解:
#     A. 全重跑（脚本 14）:           基因集 + 系数均由训练集估计   -> 0.556
#     B. 对照（本脚本）:               基因集固定、系数训练集重估     -> 待计算
#     C. 全固定（脚本 15, 0.850）:     基因集 + 系数均固定           -> 0.850
#     B - A ≈ 纯基因集选择泄漏; C - B ≈ 系数泄漏 + 训练样本量效应
#
# 协议保真度（与脚本 14_split_validation.R 严格一致，保证可比）:
#   - set.seed(20231021); B=100; idx_tr <- sample(94, floor(0.7*94)=65);
#     idx_te <- setdiff(94, idx_tr)  -> 65 训练 / 29 测试
#     （注: 任务书与审计报告表述的"66 个训练样本"为约数; 脚本 14 实际协议为
#       floor(0.7*94)=65 训练 / 29 测试, 本脚本如实按 65/29 复现）
#   - 候选池: prognostic_genes.tsv 的 48 个基因（脚本 14 同）
#   - C-index: 脚本 14/15 的手动 Harrell's C 实现（保证与已发表数字同标尺）
#   - RNG 保真度: 脚本 14 中 cv.glmnet 的折分组会消耗 RNG, 故本脚本 Branch A
#     完整复跑"筛选 -> LASSO"以消耗相同 RNG, 保证 100 次拆分与脚本 14 逐次
#     一致; 并以 n_genes/c_test 与 results/split_validation.tsv 精确校验
#     （预期 max |diff| < 1e-12, 完成折数 98）。
#
# 设计披露:
#   - 脚本 15（全固定策略, 0.850）虽用同 seed 20231021, 但其抽样方式为
#     sample(94, round(0.3*94)=28) 直接抽测试集, 与脚本 14 的"抽 65 训练集"
#     消耗不同 RNG 流 -> 两者拆分实际不同。因此本脚本增设 Branch C
#     （全固定系数在脚本 14 同拆分上评估）, 给出同拆分锚点; 论文口径的
#     0.850（脚本 15）与 Branch C 的差异反映拆分协议间的随机波动。
#
# 输入: results/tcga_expr_symbol.rds, results/tcga_clinical_surv.rds,
#       results/risk_score_tcga.tsv, results/prognostic_genes.tsv,
#       results/lasso_model.rds, results/split_validation.tsv,
#       results/split_validation_fixedmodel.tsv
# 输出: results/33_control_split.tsv          每次拆分一行（三分支）
#       results/33_control_split_summary.txt
#
# 证据上限: 内部重抽样对照, 仅分离同一队列内的泄漏成分; 不构成外部泛化证据,
# 亦不解释为因果效应。
# ============================================================================
rm(list = ls())
suppressMessages({ library(survival); library(glmnet) })

base <- Sys.getenv("ESCC_ROOT", unset = "")

if (!nzchar(base)) {

  base <- tryCatch(dirname(dirname(normalizePath(sys.frame(1)$ofile))),

                error = function(e) getwd())

}
out_dir <- file.path(base, "results")
SEED <- 20231021   # 与脚本 14/15 一致
B <- 100

## ---------- 1. 数据（与脚本 14 完全一致的对齐方式） ----------
expr <- readRDS(file.path(out_dir, "tcga_expr_symbol.rds"))   # genes x samples
rs   <- read.delim(file.path(out_dir, "risk_score_tcga.tsv"))
clin <- readRDS(file.path(out_dir, "tcga_clinical_surv.rds"))
mdl  <- readRDS(file.path(out_dir, "lasso_model.rds"))
pg   <- read.delim(file.path(out_dir, "prognostic_genes.tsv"))
ref14 <- read.delim(file.path(out_dir, "split_validation.tsv"))          # 已发表: 全重跑
ref15 <- read.delim(file.path(out_dir, "split_validation_fixedmodel.tsv")) # 已发表: 全固定（不同拆分）

common <- intersect(rs$case_id, colnames(expr))
expr <- expr[, common]
clin <- clin[match(common, clin$case_id), ]
surv_df <- data.frame(rownames = common, time = clin$os_time, status = clin$os_status)
cc <- complete.cases(surv_df$time, surv_df$status)
surv_df <- surv_df[cc, ]
expr <- expr[, surv_df$rownames]
cat(sprintf("分析样本: n=%d, 事件=%d\n", nrow(surv_df), sum(surv_df$status)))

cand <- intersect(pg$gene, rownames(expr))                 # 48 基因候选池（脚本 14 同）
X48 <- t(expr[cand, , drop = FALSE])[surv_df$rownames, ]
genes17 <- mdl$genes                                        # 全样本 17 基因（固定基因集）
coefs17 <- as.numeric(mdl$coef); names(coefs17) <- genes17  # 全样本系数（固定系数）
X17 <- t(expr[genes17, , drop = FALSE])[surv_df$rownames, ]
lp_fixed <- as.numeric(X17 %*% coefs17)
cat(sprintf("候选池: %d 基因 | 固定基因集: %d 基因 | lp_fixed vs risk_score_tcga 最大差 = %.2e\n",
            ncol(X48), length(genes17), max(abs(lp_fixed - rs$risk_score))))

## ---------- 2. C-index（脚本 14/15 的手动 Harrell's C, 保证同标尺） ----------
cidx <- function(times, status, lp) {
  ok <- !is.na(lp) & !is.na(times) & !is.na(status)
  t <- times[ok]; s <- status[ok]; r <- lp[ok]
  n <- length(t); comp <- 0; concordant <- 0; ties <- 0
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

## ---------- 3. 100 次 7:3 拆分（复刻脚本 14 协议; 三分支并行） ----------
set.seed(SEED)
res <- data.frame(rep = integer(), n_train = integer(), n_test = integer(),
                  n_events_test = integer(),
                  # Branch A: 全重跑（复刻脚本 14, 用于 RNG 保真度与同拆分参照）
                  a_n_genes = integer(), a_c_test = numeric(), a_p_km = numeric(),
                  # Branch B: 固定 17 基因, 系数训练集 Cox 重估（本脚本对照）
                  b_c_test = numeric(), b_p_km = numeric(),
                  # Branch C: 全固定（同拆分锚点, 对应脚本 15 的 0.850 口径）
                  c_c_test = numeric(), c_p_km = numeric(),
                  stringsAsFactors = FALSE)
n_b_fail <- 0
for (b in seq_len(B)) {
  idx_tr <- sample(seq_len(nrow(X48)), size = floor(0.7 * nrow(X48)))   # 65
  idx_te <- setdiff(seq_len(nrow(X48)), idx_tr)                          # 29
  Xtr <- X48[idx_tr, , drop = FALSE]; Xte <- X48[idx_te, , drop = FALSE]
  ytr <- Surv(surv_df$time[idx_tr], surv_df$status[idx_tr])
  yte <- Surv(surv_df$time[idx_te], surv_df$status[idx_te])

  ## Branch A —— 与脚本 14 逐行一致的"筛选 -> LASSO"（含 RNG 消耗）
  a_genes <- NA_integer_; a_c <- NA_real_; a_p <- NA_real_
  keep <- sapply(seq_len(ncol(Xtr)), function(k) {
    f <- try(coxph(ytr ~ Xtr[, k]), silent = TRUE)
    if (inherits(f, "try-error")) return(FALSE)
    summary(f)$coefficients[1, 5] < 0.05
  })
  if (sum(keep) >= 3) {
    Xtr2 <- Xtr[, keep, drop = FALSE]
    cvf <- try(suppressWarnings(cv.glmnet(Xtr2, ytr, family = "cox", nfolds = 10,
                                          alpha = 1, ties = "breslow", maxit = 5000)),
               silent = TRUE)
    if (!inherits(cvf, "try-error")) {
      co <- coef(cvf, s = "lambda.min")
      gsel <- rownames(co)[which(as.numeric(co) != 0)]
      if (length(gsel) >= 2) {
        lp_te <- as.numeric(predict(cvf, Xte[, colnames(Xtr2), drop = FALSE],
                                    s = "lambda.min", type = "link"))
        a_genes <- length(gsel)
        a_c <- cidx(surv_df$time[idx_te], surv_df$status[idx_te], lp_te)
        cut <- median(as.numeric(predict(cvf, Xtr2, s = "lambda.min", type = "link")))
        grp <- ifelse(lp_te > cut, "High", "Low")
        sd_ <- try(survdiff(yte ~ grp), silent = TRUE)
        a_p <- if (inherits(sd_, "try-error")) NA_real_ else
               pchisq(as.numeric(sd_$chisq), df = 1, lower.tail = FALSE)
      }
    }
  }

  ## Branch B —— 固定 17 基因, 系数仅在 65 个训练样本内 Cox 重估
  b_c <- NA_real_; b_p <- NA_real_
  dtr <- data.frame(X17[idx_tr, , drop = FALSE])
  dtr$time <- surv_df$time[idx_tr]; dtr$status <- surv_df$status[idx_tr]
  fb <- try(suppressWarnings(coxph(Surv(time, status) ~ ., data = dtr)), silent = TRUE)
  if (!inherits(fb, "try-error")) {
    lp_te_b <- as.numeric(predict(fb, newdata = data.frame(X17[idx_te, , drop = FALSE]),
                                  type = "lp"))
    b_c <- cidx(surv_df$time[idx_te], surv_df$status[idx_te], lp_te_b)
    cut_b <- median(as.numeric(predict(fb, type = "lp")))
    grp_b <- ifelse(lp_te_b > cut_b, "High", "Low")
    sd_b <- try(survdiff(yte ~ grp_b), silent = TRUE)
    b_p <- if (inherits(sd_b, "try-error")) NA_real_ else
           pchisq(as.numeric(sd_b$chisq), df = 1, lower.tail = FALSE)
  } else n_b_fail <- n_b_fail + 1

  ## Branch C —— 全固定（基因 + 系数均用全样本 94 例的值）, 同拆分评估
  lp_te_c <- lp_fixed[idx_te]
  c_c <- cidx(surv_df$time[idx_te], surv_df$status[idx_te], lp_te_c)
  grp_c <- ifelse(lp_te_c > median(lp_fixed), "High", "Low")
  sd_c <- try(survdiff(yte ~ grp_c), silent = TRUE)
  c_p <- if (inherits(sd_c, "try-error")) NA_real_ else
         pchisq(as.numeric(sd_c$chisq), df = 1, lower.tail = FALSE)

  res <- rbind(res, data.frame(rep = b, n_train = length(idx_tr), n_test = length(idx_te),
                               n_events_test = sum(surv_df$status[idx_te]),
                               a_n_genes = a_genes, a_c_test = a_c, a_p_km = a_p,
                               b_c_test = b_c, b_p_km = b_p,
                               c_c_test = c_c, c_p_km = c_p,
                               stringsAsFactors = FALSE))
  if (b %% 20 == 0) cat(sprintf("  完成 %d/%d\n", b, B))
}

## ---------- 4. 保真度校验: Branch A 必须精确复现脚本 14 ----------
m <- merge(res[!is.na(res$a_c_test), ], ref14, by = "rep")
fid_genes <- identical(m$a_n_genes, m$n_genes)
fid_c <- max(abs(m$a_c_test - m$c_test))
cat(sprintf("\n[保真度] Branch A 完成 %d 次（参考 split_validation.tsv: %d）; 合并 %d 次\n",
            sum(!is.na(res$a_c_test)), nrow(ref14), nrow(m)))
cat(sprintf("[保真度] n_genes 完全一致: %s | c_test 最大绝对差 = %.3e\n",
            fid_genes, fid_c))

## ---------- 5. 汇总 ----------
summ_row <- function(x, label) {
  x <- x[!is.na(x)]
  sprintf("%s: 完成 %d/100 | 中位 = %.4f (IQR %.4f-%.4f) | 均值 = %.4f (SD %.4f) | 范围 %.4f-%.4f | >0.5: %.1f%% | >0.7: %.1f%%",
          label, length(x), median(x), quantile(x, .25), quantile(x, .75),
          mean(x), sd(x), min(x), max(x),
          100 * mean(x > 0.5), 100 * mean(x > 0.7))
}
cat("\n", summ_row(res$a_c_test, "A 全重跑（同拆分, 复刻脚本 14）"), "\n", sep = "")
cat("", summ_row(res$b_c_test, "B 固定基因集+系数训练集重估（本对照）"), "\n", sep = "")
cat("", summ_row(res$c_c_test, "C 全固定（同拆分锚点）"), "\n", sep = "")
cat("", summ_row(ref15$c_test, "参考: 脚本 15 全固定（已发表, 不同拆分）"), "\n", sep = "")

medA <- median(res$a_c_test, na.rm = TRUE)
medB <- median(res$b_c_test, na.rm = TRUE)
medC <- median(res$c_c_test, na.rm = TRUE)
med15 <- median(ref15$c_test)
cat(sprintf("\n三向分解（同拆分中位数）: A=%.4f -> B=%.4f -> C=%.4f\n", medA, medB, medC))
cat(sprintf("  C - B = %.4f ≈ 系数泄漏 + 训练样本量效应（94 vs 65）\n", medC - medB))
cat(sprintf("  B - A = %.4f ≈ 纯基因集选择泄漏（基因集用含测试样本的全样本选定）\n", medB - medA))
cat(sprintf("已发表口径: 0.850（脚本 15, 不同拆分） - B = %.4f ≈ 系数泄漏+样本量效应; B - 0.556（脚本 14） = %.4f ≈ 基因集泄漏\n",
            med15 - medB, medB - 0.5558))
cat(sprintf("Branch B Cox 拟合失败: %d/100\n", n_b_fail))
cat(sprintf("Branch B log-rank P<0.05 比例: %.1f%% | Branch C: %.1f%%\n",
            100 * mean(res$b_p_km < 0.05, na.rm = TRUE),
            100 * mean(res$c_p_km < 0.05, na.rm = TRUE)))

## ---------- 6. 输出 ----------
write.table(res, file.path(out_dir, "33_control_split.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

summ <- c(
  "===== 33_control_split.R 汇总 =====",
  paste0("运行时间: ", format(Sys.time())),
  "",
  "--- 设计与协议 ---",
  sprintf("队列: TCGA-ESCC n=%d, 事件=%d（与脚本 14 一致）", nrow(surv_df), sum(surv_df$status)),
  sprintf("协议: set.seed(%d), 100 次 7:3 拆分, floor(0.7*94)=65 训练 / 29 测试（复刻脚本 14）", SEED),
  "候选池: 48 基因（prognostic_genes.tsv, 脚本 14 同）; 固定基因集: 17 基因（lasso_model.rds 全样本 LASSO）",
  "C-index: 脚本 14/15 手动 Harrell's C（与已发表 0.556/0.850 同标尺）",
  "",
  "--- 保真度校验（Branch A vs results/split_validation.tsv） ---",
  sprintf("Branch A 完成 %d 次（参考 %d）; n_genes 完全一致: %s; c_test 最大绝对差 = %.3e",
          sum(!is.na(res$a_c_test)), nrow(ref14), fid_genes, fid_c),
  "",
  "--- 三分支测试集 C-index ---",
  summ_row(res$a_c_test, "A 全重跑（同拆分, 复刻脚本 14）"),
  summ_row(res$b_c_test, "B 固定基因集+系数训练集重估（本对照）"),
  summ_row(res$c_c_test, "C 全固定（同拆分锚点）"),
  summ_row(ref15$c_test, "参考: 脚本 15 全固定（已发表 0.850, 不同拆分）"),
  "",
  "--- 分解解读（中位数） ---",
  sprintf("同拆分三向分解: A=%.4f -> B=%.4f -> C=%.4f", medA, medB, medC),
  sprintf("C - B = %.4f ≈ 系数泄漏 + 训练样本量效应（94 vs 65 的估计精度差异）", medC - medB),
  sprintf("B - A = %.4f ≈ 纯基因集选择泄漏（基因集用全部 94 样本、含 29 个测试样本选定）", medB - medA),
  sprintf("已发表口径分解: 0.850 - B = %.4f ≈ 系数泄漏+样本量效应; B - 0.556 = %.4f ≈ 基因集泄漏",
          med15 - medB, medB - 0.5558),
  "",
  "--- 设计披露 ---",
  "1) 脚本 15（0.850）与脚本 14（0.556）虽同 seed 20231021, 但抽样方式不同（直接抽 28 个测试集 vs 抽 65 个训练集）,",
  "   RNG 流不同 -> 实际拆分不同; 故 0.850 与 0.556 之差含拆分间随机波动。Branch C 给出同拆分锚点。",
  "2) 任务书/审计报告中\"66 个训练样本\"为约数; 脚本 14 实际协议为 floor(0.7*94)=65 训练 / 29 测试。",
  "3) Branch B 的 Cox 为 17 基因 / 65 样本 / ~21 事件（EPV≈1.2）, 系数估计高度不稳定是设计意图——",
  "   这正是被量化的\"系数泄漏+样本量效应\"本身, 不是缺陷。",
  "",
  "--- 证据上限 ---",
  "内部重抽样对照: 仅分离同一队列内的泄漏成分, 不构成外部泛化证据, 不解释为因果效应。",
  "",
  "--- 复现信息 ---",
  paste(capture.output(sessionInfo()), collapse = "\n")
)
writeLines(summ, file.path(out_dir, "33_control_split_summary.txt"))
cat("\n", paste(summ, collapse = "\n"), "\n")
cat("\n33_control_split 完成\n")
