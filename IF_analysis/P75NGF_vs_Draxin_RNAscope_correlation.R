# Example input vectors
draxin_dots_per_cell <- c(3.7, 5.9, 2.05, 1.99, 0.61, 1.72, 0.71, 0.30, 0.67, 2.70, 2.78, 1.97, 1.00, 9.90, 3.96, 2.80, 5.45)
p75_percent <- c(69.58, 96.61, 78.57, 17.68, 14.14, 0.57, 4.73, 16.55, 15.46, 1.79, 0.00, 2.22, 6.57, 99.1, 0.65, 2.7, 22.77)
pS100_percent <- c(76, 84.75, 50, 83.44, 58.58, 39.8, 90.2, 89.52, 85.22, 91.01, 75.38, 6.66, 87.78, 93.2, 76.97, 30.55, 92.07)

# Build data frame
df <- data.frame(
  P75_percent = p75_percent,
  DRAXIN_dots_per_cell = draxin_dots_per_cell,
  PS100_percent = pS100_percent
)

library(ggplot2)
library(patchwork)

# Spearman correlation for P75
ct <- cor.test(
  df$P75_percent,
  df$DRAXIN_dots_per_cell,
  method = "spearman",
  exact = FALSE
)

rho <- unname(ct$estimate)
pval <- ct$p.value
lab <- paste0(
  "Spearman \u03C1 = ", sprintf("%.2f", rho),
  "\n",
  "p = ", signif(pval, 2)
)

# Spearman correlation for S100B
ct_S100b <- cor.test(
  df$PS100_percent,
  df$DRAXIN_dots_per_cell,
  method = "spearman",
  exact = FALSE
)

rho_s100b <- unname(ct_S100b$estimate)
pval_s100b <- ct_S100b$p.value
lab_s100b <- paste0(
  "Spearman \u03C1 = ", sprintf("%.2f", rho_s100b),
  "\n",
  "p = ", signif(pval_s100b, 2)
)

# Plot for P75
p <- ggplot(df, aes(x = P75_percent, y = DRAXIN_dots_per_cell)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = lab,
    hjust = -0.05, vjust = 1.1,
    size = 3.5
  ) +
  theme_bw() +
  labs(
    x = "% P75NGFR+ cells per nerve",
    y = "Mean DRAXIN dots per cell",
    title = "P75NGFR vs DRAXIN"
  )

# Plot for S100B
pS100b <- ggplot(df, aes(x = PS100_percent, y = DRAXIN_dots_per_cell)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  annotate(
    "text",
    x = -Inf, y = Inf,
    label = lab_s100b,
    hjust = -0.05, vjust = 1.1,
    size = 3.5
  ) +
  theme_bw() +
  labs(
    x = "% PS100B+ cells per nerve",
    y = "Mean DRAXIN dots per cell",
    title = "PS100B vs DRAXIN"
  )

# Show plots
p
pS100b

# Combine plots side by side
combined_plot <- p + pS100b

# Save combined plot
ggsave(
  filename = "~/NAS/Aleksandra/Computational/R programming/paper_projects/GSEA/P75NGFvsDRAXINcorr.pdf",
  plot = combined_plot,
  width = 10,
  height = 5
)
