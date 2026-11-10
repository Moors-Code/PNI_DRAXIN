# =============================================================================
# Figure 2E – SuperPlot
# Small dots = nerves (coloured by patient); large dots = patient medians,
# connected per patient; brackets = Holm-adjusted p values from fig2E_stats.R.
#
# Run from the project folder, after fig2E_stats.R:  Rscript scripts/fig2E_plot.R
#   in : output/fig2E_per_nerve.csv, output/fig2E_model_results.csv
#   out: output/Fig2E_cd45_superplot.pdf           (both panels)
#        output/Fig2E_cd45_superplot_full.pdf      (full perineural area)
#        output/Fig2E_cd45_superplot_tumorfree.pdf (non-tumour perineural area)
# Packages: install.packages("ggplot2")
# =============================================================================

# `%||%` is only in base R >= 4.4; some package versions use it, so define it for older R
if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x

library(ggplot2)

if (basename(getwd()) == "scripts") setwd("..")
set.seed(1)

GROUPS <- c("NI-Adj", "NI-TB", "InvN")
# Patient 3 = blue; no red/green (reserved for patients 1 and 2)
PATIENT_COLORS <- c("Patient 3" = "#1f77b4", "Patient 4" = "#ff7f0e",
                    "Patient 5" = "#9467bd", "Patient 6" = "#8c564b", "Patient 7" = "#17becf")
PANELS <- c(full = "Full perineural area", tumorfree = "Non-tumour perineural area")

d   <- read.csv("output/fig2E_per_nerve.csv")
res <- read.csv("output/fig2E_model_results.csv")

# nerves: one row per nerve per panel
nerves <- rbind(data.frame(patient = d$patient_label, group = d$group, panel = PANELS[["full"]],      density = d$density_full),
                data.frame(patient = d$patient_label, group = d$group, panel = PANELS[["tumorfree"]], density = d$density_tumorfree))
nerves$group <- factor(nerves$group, levels = GROUPS)
nerves$panel <- factor(nerves$panel, levels = PANELS)
nerves$x     <- as.numeric(nerves$group) + runif(nrow(nerves), -0.28, 0.28)

# patient medians, slightly offset per patient
medians <- aggregate(density ~ patient + group + panel, data = nerves, FUN = median)
offset  <- setNames(seq(-0.12, 0.12, length.out = 4), names(PATIENT_COLORS)[1:4])
medians$x <- as.numeric(medians$group) + offset[medians$patient]
group_med <- aggregate(density ~ group + panel, data = medians, FUN = median)
group_med$x <- as.numeric(group_med$group)

# brackets: Holm-adjusted p values from fig2E_stats.R
first  <- sub(" vs .*", "", res$comparison)
second <- sub(".* vs ", "", res$comparison)
pos    <- match(c(first), GROUPS); pos2 <- match(second, GROUPS)
res$panel <- factor(res$panel, levels = PANELS)
res$x0 <- pmin(pos, pos2); res$x1 <- pmax(pos, pos2)
wide   <- (res$x1 - res$x0) > 1                       # NI-Adj – InvN bracket goes on top
res$x0 <- ifelse(wide, res$x0, res$x0 + 0.05)         # small gap between neighbouring brackets
res$x1 <- ifelse(wide, res$x1, res$x1 - 0.05)
res$y  <- ifelse(wide, 1.2e7, 2.5e6)
res$label <- ifelse(res$p_holm < 0.001, "p < 0.001", sprintf("p = %.2g", res$p_holm))

superplot <- function(panels) {
  keep <- function(df) df[df$panel %in% panels, ]
  ggplot() +
    geom_point(data = keep(nerves), aes(x, density, colour = patient), size = 0.9, alpha = 0.35, stroke = 0) +
    geom_line(data = keep(medians), aes(x, density, colour = patient, group = patient), linewidth = 0.6) +
    geom_point(data = keep(medians), aes(x, density, fill = patient), shape = 21, size = 3, stroke = 0.4) +
    geom_segment(data = keep(group_med), aes(x = x - 0.3, xend = x + 0.3, y = density, yend = density), linewidth = 1) +
    geom_segment(data = keep(res), aes(x = x0, xend = x1, y = y, yend = y), linewidth = 0.3) +
    geom_segment(data = keep(res), aes(x = x0, xend = x0, y = y, yend = y / 1.6), linewidth = 0.3) +
    geom_segment(data = keep(res), aes(x = x1, xend = x1, y = y, yend = y / 1.6), linewidth = 0.3) +
    geom_text(data = keep(res), aes(x = (x0 + x1) / 2, y = y * 1.5, label = label), size = 3) +
    facet_wrap(~ panel, nrow = 1) +
    scale_colour_manual(values = PATIENT_COLORS, name = "Patient") +
    scale_fill_manual(values = PATIENT_COLORS, name = "Patient") +
    scale_x_continuous(breaks = 1:3, labels = GROUPS, limits = c(0.5, 3.5)) +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 200, base = 10),  # log scale that shows zeros
                       limits = c(0, 5e7), breaks = c(0, 1e3, 1e4, 1e5, 1e6),
                       labels = c("0", expression(10^3), expression(10^4), expression(10^5), expression(10^6))) +
    labs(x = NULL, y = expression(CD45^"+" ~ "density (cells / cm"^2*")")) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
}

ggsave("output/Fig2E_cd45_superplot.pdf",           superplot(PANELS),              width = 9, height = 4.3)
ggsave("output/Fig2E_cd45_superplot_full.pdf",      superplot(PANELS["full"]),      width = 5, height = 4.3)
ggsave("output/Fig2E_cd45_superplot_tumorfree.pdf", superplot(PANELS["tumorfree"]), width = 5, height = 4.3)
cat("Figures saved in output/\n")

