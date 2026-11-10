# ============================================================
# 04_P75NGF_vs_Draxin_RNAscope_correlation.R
#
# Spearman correlation between DRAXIN RNAscope signal (dots per cell)
# and two IF markers (P75NGFR%, PS100B%) per nerve.
#
# Input : input/IF_HCR_DRAXIN/P75NGF_vs_Draxin_RNAscope_data.csv
# Output: output/IF_HCR_DRAXIN/tables/correlation_stats.csv
#         output/IF_HCR_DRAXIN/figures/P75NGFvsDRAXINcorr.pdf
# ============================================================

suppressPackageStartupMessages({
  library(readr)
  library(ggplot2)
  library(patchwork)
})

# ------------------------------------------------------------
# paths
# ------------------------------------------------------------
data_file <- file.path("input", "IF_HCR_DRAXIN", "P75NGF_vs_Draxin_RNAscope_data.csv")

tables_dir  <- file.path("output", "IF_HCR_DRAXIN", "tables")
figures_dir <- file.path("output", "IF_HCR_DRAXIN", "figures")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(data_file)) {
  stop("Data file not found: ", data_file)
}

df <- readr::read_csv(data_file, show_col_types = FALSE)

# ------------------------------------------------------------
# Spearman correlation: DRAXIN vs P75NGFR
# ------------------------------------------------------------
ct <- cor.test(df$P75_percent, df$DRAXIN_dots_per_cell, method = "spearman", exact = FALSE)
rho <- unname(ct$estimate)
pval <- ct$p.value
lab <- paste0("Spearman \u03C1 = ", sprintf("%.2f", rho), "\n", "p = ", signif(pval, 2))

# ------------------------------------------------------------
# Spearman correlation: DRAXIN vs PS100B
# ------------------------------------------------------------
ct_S100b <- cor.test(df$PS100_percent, df$DRAXIN_dots_per_cell, method = "spearman", exact = FALSE)
rho_s100b <- unname(ct_S100b$estimate)
pval_s100b <- ct_S100b$p.value
lab_s100b <- paste0("Spearman \u03C1 = ", sprintf("%.2f", rho_s100b), "\n", "p = ", signif(pval_s100b, 2))

# ------------------------------------------------------------
# save correlation stats table
# ------------------------------------------------------------
stats_tbl <- tibble::tibble(
  comparison = c("P75NGFR_vs_DRAXIN", "PS100B_vs_DRAXIN"),
  spearman_rho = c(rho, rho_s100b),
  p_value = c(pval, pval_s100b),
  n = c(nrow(df), nrow(df))
)
readr::write_csv(stats_tbl, file.path(tables_dir, "correlation_stats.csv"))

# ------------------------------------------------------------
# plots
# ------------------------------------------------------------
p <- ggplot(df, aes(x = P75_percent, y = DRAXIN_dots_per_cell)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, label = lab, hjust = -0.05, vjust = 1.1, size = 3.5) +
  theme_bw() +
  labs(x = "% P75NGFR+ cells per nerve", y = "Mean DRAXIN dots per cell", title = "P75NGFR vs DRAXIN")

pS100b <- ggplot(df, aes(x = PS100_percent, y = DRAXIN_dots_per_cell)) +
  geom_point(size = 2, alpha = 0.85) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linewidth = 0.7) +
  annotate("text", x = -Inf, y = Inf, label = lab_s100b, hjust = -0.05, vjust = 1.1, size = 3.5) +
  theme_bw() +
  labs(x = "% PS100B+ cells per nerve", y = "Mean DRAXIN dots per cell", title = "PS100B vs DRAXIN")

print(p)
print(pS100b)

combined_plot <- p + pS100b

ggsave(
  filename = file.path(figures_dir, "P75NGFvsDRAXINcorr.pdf"),
  plot = combined_plot,
  width = 10,
  height = 5
)

message("Done.")
message("Table written to:  ", file.path(tables_dir, "correlation_stats.csv"))
message("Figure written to: ", file.path(figures_dir, "P75NGFvsDRAXINcorr.pdf"))
