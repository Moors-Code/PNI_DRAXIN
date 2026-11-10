# =============================================================================
# Nerve markers – SuperPlots (same style as Figure 2E)
# Small dots = nerves (coloured by patient); large dots = patient medians,
# connected per patient; black bar = median of patient medians.
#
# Main figures:          NI-TB vs InvN (same section) with the mixed-model p value.
# Supplementary figures: NI-Adj, NI-TB, InvN, descriptive only (no statistics), because
#                        NI-Adj nerves lie on a separate section.
#
# Run from the project folder, after markers_stats.R:  Rscript scripts/markers_plot.R
#   in : output/markers_per_nerve.csv, output/markers_model_results.csv
#   out: output/markers_<analysis>.pdf          main, one per analysis
#        output/markers_all.pdf                 main, all panels
#        output/supp_markers_<analysis>.pdf     supplementary, three groups
#        output/supp_markers_all.pdf
# Packages: install.packages("ggplot2")
# =============================================================================

if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x
library(ggplot2)

if (basename(getwd()) == "scripts") setwd("..")
set.seed(1)

# Patient 3 = blue; no red/green (reserved for patients 1 and 2)
PATIENT_COLORS <- c("Patient 3" = "#1f77b4", "Patient 4" = "#ff7f0e",
                    "Patient 5" = "#9467bd", "Patient 6" = "#8c564b", "Patient 7" = "#17becf")
Y_LABEL <- "Normalised intensity\n(relative to section background, log scale)"

d_all <- read.csv("output/markers_per_nerve.csv")
res   <- read.csv("output/markers_model_results.csv")
labels <- unique(d_all$label)
ids    <- unique(d_all$analysis)

superplot <- function(groups, panels, with_stats) {
  d <- d_all[d_all$group %in% groups & d_all$label %in% panels, ]
  d$group <- factor(d$group, levels = groups)
  d$label <- factor(d$label, levels = labels)
  d$x <- as.numeric(d$group) + runif(nrow(d), -0.28, 0.28)

  medians <- aggregate(normalised ~ patient_label + group + label, data = d, FUN = median)
  pats    <- intersect(names(PATIENT_COLORS), unique(d$patient_label))
  offset  <- setNames(seq(-0.14, 0.14, length.out = length(pats)), pats)
  medians$x <- as.numeric(medians$group) + offset[medians$patient_label]
  group_med <- aggregate(normalised ~ group + label, data = medians, FUN = median)
  group_med$x <- as.numeric(group_med$group)

  p <- ggplot() +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_point(data = d, aes(x, normalised, colour = patient_label), size = 0.9, alpha = 0.35, stroke = 0) +
    geom_line(data = medians, aes(x, normalised, colour = patient_label, group = patient_label), linewidth = 0.6) +
    geom_point(data = medians, aes(x, normalised, fill = patient_label), shape = 21, size = 3, stroke = 0.4) +
    geom_segment(data = group_med, aes(x = x - 0.3, xend = x + 0.3, y = normalised, yend = normalised), linewidth = 1)

  if (with_stats) {
    top <- tapply(d$normalised, d$label, max)
    br  <- res[res$label %in% panels, ]
    br$label <- factor(br$label, levels = labels)
    br$y <- top[as.character(br$label)] * 1.6
    br$p_label <- ifelse(br$p_value < 0.001, "p < 0.001", sprintf("p = %.2g", br$p_value))
    p <- p +
      geom_segment(data = br, aes(x = 1, xend = 2, y = y, yend = y), linewidth = 0.3) +
      geom_segment(data = br, aes(x = 1, xend = 1, y = y, yend = y / 1.12), linewidth = 0.3) +
      geom_segment(data = br, aes(x = 2, xend = 2, y = y, yend = y / 1.12), linewidth = 0.3) +
      geom_text(data = br, aes(x = 1.5, y = y * 1.18, label = p_label), size = 3)
  }

  p + facet_wrap(~ label, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = PATIENT_COLORS, name = "Patient") +
    scale_fill_manual(values = PATIENT_COLORS, name = "Patient") +
    scale_x_continuous(breaks = seq_along(groups), labels = groups, limits = c(0.5, length(groups) + 0.5)) +
    scale_y_log10() +
    labs(x = NULL, y = Y_LABEL) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
}

MAIN <- c("NI-TB", "InvN")
SUPP <- c("NI-Adj", "NI-TB", "InvN")
for (i in seq_along(ids)) {
  ggsave(sprintf("output/markers_%s.pdf", ids[i]),      superplot(MAIN, labels[i], TRUE),  width = 4,   height = 4.3)
  ggsave(sprintf("output/supp_markers_%s.pdf", ids[i]), superplot(SUPP, labels[i], FALSE), width = 4.8, height = 4.3)
}
ggsave("output/markers_all.pdf",      superplot(MAIN, labels, TRUE),  width = 3.2 * length(labels) + 1.5, height = 4.3)
ggsave("output/supp_markers_all.pdf", superplot(SUPP, labels, FALSE), width = 4   * length(labels) + 1.5, height = 4.3)
cat("Figures saved in output/\n")

