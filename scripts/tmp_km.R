suppressMessages(library(survival))
d <- read.delim("results/risk_score_tcga.tsv")
d$grp <- factor(ifelse(d$risk_score >= median(d$risk_score), "High", "Low"), levels = c("Low", "High"))
s <- survdiff(Surv(os_time, os_status) ~ grp, data = d)
cat("KM chisq =", as.numeric(s$statistic), " P =", signif(s$p.value, 3), "\n")
cat("High: n =", sum(d$grp == "High"), " deaths =", sum(d$os_status[d$grp == "High"]),
    " median OS =", median(d$os_time[d$grp == "High"]), "\n")
cat("Low:  n =", sum(d$grp == "Low"), " deaths =", sum(d$os_status[d$grp == "Low"]),
    " median OS =", median(d$os_time[d$grp == "Low"]), "\n")
fit <- survfit(Surv(os_time, os_status) ~ grp, data = d)
sf <- summary(fit, times = c(365, 1095, 1825))
st <- names(sf$strata)
for (i in seq_along(st)) {
  cat(st[i], " 1y =", round(sf$surv[(i - 1) * 3 + 1] * 100, 1),
      " 3y =", round(sf$surv[(i - 1) * 3 + 2] * 100, 1),
      " 5y =", round(sf$surv[(i - 1) * 3 + 3] * 100, 1), "\n")
}
clin <- readRDS("results/tcga_clinical_surv.rds")
cat("gender male =", sum(clin$gender == "male", na.rm = TRUE),
    " female =", sum(clin$gender == "female", na.rm = TRUE), "\n")
cat("age median =", median(clin$age_at_index, na.rm = TRUE),
    " range =", range(clin$age_at_index, na.rm = TRUE), "\n")
