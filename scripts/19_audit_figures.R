## ============================================================
## Visualisation of the literature-audit results
## ============================================================
suppressMessages({ library(ggplot2); library(patchwork) })

d <- read.csv("results/literature_audit/audit_records.csv", stringsAsFactors = FALSE)
mod <- d[d$is_model_study == 1 & d$has_abstract == 1, ]
n <- nrow(mod)
cat("Included studies:", n, "\n")

# ---------- Panel A: validation-practice distribution ----------
vpr <- data.frame(
  practice = c("External validation\n(independent cohort)", "Internal split/\ncross-validation", "Both mentioned",
               "Neither mentioned", "Explicit re-run of\nmodelling pipeline"),
  pct = c(100 * sum(mod$external_validation) / n,
          100 * sum(mod$internal_validation) / n,
          100 * sum(mod$external_validation & mod$internal_validation) / n,
          100 * sum(!mod$external_validation & !mod$internal_validation) / n,
          100 * sum(mod$rerun_pipeline) / n),
  n = c(sum(mod$external_validation), sum(mod$internal_validation),
        sum(mod$external_validation & mod$internal_validation),
        sum(!mod$external_validation & !mod$internal_validation),
        sum(mod$rerun_pipeline))
)
vpr$practice <- factor(vpr$practice, levels = rev(vpr$practice))
vpr$col <- c("#B2182B", "#2E6DA4", "#6F42C1", "#E64B35", "#C79100")[c(5, 4, 3, 2, 1)]

pA <- ggplot(vpr, aes(x = practice, y = pct, fill = practice)) +
  geom_bar(stat = "identity", width = 0.72) +
  geom_text(aes(label = paste0(round(pct, 1), "%\n(", n, " studies)")),
            hjust = -0.12, size = 3.4) +
  coord_flip() + ylim(0, 100) +
  scale_fill_manual(values = c("External validation\n(independent cohort)" = "#4C9F70",
                               "Internal split/\ncross-validation" = "#2E6DA4",
                               "Both mentioned" = "#6F42C1",
                               "Neither mentioned" = "#B2182B",
                               "Explicit re-run of\nmodelling pipeline" = "#0C7C94")) +
  labs(x = NULL, y = "Proportion of included studies (%)",
       title = paste0("A  Validation practices (n = ", n, " studies)")) +
  theme_bw(base_size = 12) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12, face = "bold"),
        panel.grid.major.y = element_blank())

# ---------- Panel B: trend by year ----------
yr <- as.data.frame(table(mod$year))
names(yr) <- c("year", "n")
yr$year <- as.integer(as.character(yr$year))
yr <- yr[yr$year >= 2017, ]
pB <- ggplot(yr, aes(x = year, y = n)) +
  geom_bar(stat = "identity", fill = "#2E6DA4", width = 0.7) +
  geom_text(aes(label = n), vjust = -0.4, size = 3.2) +
  labs(x = "Publication year", y = "Number of studies",
       title = "B  Publication trend of ESCC prognostic-model studies") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"))

# ---------- Panel C: distribution of reported AUCs ----------
aucs <- as.numeric(mod$max_auc[mod$max_auc != ""])
aucs <- aucs[!is.na(aucs)]
aucd <- data.frame(auc = aucs)
pC <- ggplot(aucd, aes(x = auc)) +
  geom_histogram(binwidth = 0.05, fill = "#C79100", color = "white") +
  geom_vline(xintercept = median(aucs), color = "#B2182B", linetype = 2, linewidth = 1) +
  geom_vline(xintercept = 0.8, color = "#2E7D46", linetype = 3) +
  annotate("text", x = median(aucs), y = Inf, vjust = 1.5, hjust = -0.15,
           label = paste0("Median ", round(median(aucs), 3)), size = 3.6, color = "#B2182B") +
  labs(x = "Reported AUC / C-index", y = "Number of studies",
       title = paste0("C  Reported discrimination (n = ", length(aucs), " studies extractable)")) +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 12, face = "bold"))

combined <- pA / (pB | pC)
ggsave("figures/fig12_literature_audit.png", combined,
       width = 11, height = 10, dpi = 300)
cat("Figure saved: figures/fig12_literature_audit.png\n")
cat("Median AUC:", round(median(aucs), 3), "| n =", length(aucs), "\n")
