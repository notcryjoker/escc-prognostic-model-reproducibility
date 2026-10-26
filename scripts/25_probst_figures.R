suppressMessages({ library(ggplot2); library(patchwork) })
d <- read.csv("results/literature_audit/probst_assessment.csv")
n <- nrow(d)
cat("n =", n, "\n")

# A: overall risk of bias
ov <- data.frame(
  risk = factor(c("Low", "Unclear", "High"), levels = c("Low", "Unclear", "High")),
  pct = c(100*sum(d$probast_overall=="Low")/n,
          100*sum(d$probast_overall=="Unclear")/n,
          100*sum(d$probast_overall=="High")/n),
  cnt = c(sum(d$probast_overall=="Low"), sum(d$probast_overall=="Unclear"), sum(d$probast_overall=="High")))
pA <- ggplot(ov, aes(x = risk, y = pct, fill = risk)) +
  geom_bar(stat = "identity", width = .65) +
  geom_text(aes(label = paste0(round(pct,1), "%\n(", cnt, ")")), vjust = -.4, size = 3.5) +
  ylim(0, 100) +
  scale_fill_manual(values = c("Low" = "#4C9F70", "Unclear" = "#C79100", "High" = "#B2182B")) +
  labs(x = NULL, y = "Proportion (%)", title = paste0("A  PROBAST overall risk of bias (n = ", n, ")")) +
  theme_bw(base_size = 11) + theme(legend.position = "none",
    plot.title = element_text(size = 12, face = "bold"))

# B: proportion of items rated at high risk of bias
it <- data.frame(
  item = c("4.8 Validation adequacy", "4.7 Performance reporting\n(discrimination + calibration)",
           "4.5 Predictor selection", "4.1 Adequacy of events", "4.4 Missing data handling"),
  high = c(100*sum(d$probast_4_8_validation=="High")/n,
           100*sum(d$probast_4_7_performance=="High")/n,
           100*sum(d$probast_4_5_selection=="High")/n,
           100*sum(d$probast_4_1_epv=="High")/n,
           100*sum(d$probast_4_4_missing=="High")/n),
  unclear = c(100*sum(d$probast_4_8_validation=="Unclear")/n,
              100*sum(d$probast_4_7_performance=="Unclear")/n,
              100*sum(d$probast_4_5_selection=="Unclear")/n,
              100*sum(d$probast_4_1_epv=="Unclear")/n,
              100*sum(d$probast_4_4_missing=="Unclear")/n))
it$item <- factor(it$item, levels = rev(it$item))
pB <- ggplot(it, aes(x = item, y = high)) +
  geom_bar(stat = "identity", width = .68, fill = "#B2182B") +
  geom_text(aes(label = paste0(round(high,1), "%")), hjust = -.12, size = 3.3) +
  coord_flip() + ylim(0, 100) +
  labs(x = NULL, y = "Proportion rated at high risk of bias (%)", title = "B  PROBAST items rated at high risk of bias") +
  theme_bw(base_size = 11) + theme(plot.title = element_text(size = 12, face = "bold"),
    panel.grid.major.y = element_blank())

ggsave("figures/fig13_probst.png", pA | pB, width = 11, height = 5.2, dpi = 300)
cat("Saved figures/fig13_probst.png\n")
cat("Overall High:", sum(d$probast_overall=="High"), "| Low:", sum(d$probast_overall=="Low"), "\n")
