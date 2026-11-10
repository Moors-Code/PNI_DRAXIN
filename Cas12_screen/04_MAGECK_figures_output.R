# ======================================================================
# MAGeCK_RRA_and_volcano_plots.R
#
# Re-draws the MAGeCK figures from the existing gene summary (MAGeCK is NOT
# re-run, so all values are exactly those of the original analysis):
#   1. RRA score distribution, negatively selected (depleted) genes
#   2. RRA score distribution, positively selected (enriched) genes
#      (genes ranked by RRA score, log scale, best on top; the top N_LABEL
#      genes of each direction -- MAGeCK's own rank -- as coloured points with
#      a legend, as in MAGeCK's default plots)
#   3. Volcano: -log10 FDR vs log2 fold change (Draxin vs control). Each gene's
#      FDR is taken from the direction of its fold change (neg|fdr for
#      LFC < 0, pos|fdr for LFC > 0), as MAGeCK tests the two tails separately.
#
# Figures are drawn at final print size (pt; 72 pt = 1 inch) with editable
# Helvetica text (base pdf device), for Illustrator.
#
# INPUT : <MAGECK_DIR>/T0_PLUS_CONTROL_PAIRED_TOTAL.gene_summary.txt
# OUTPUT: <MAGECK_DIR>/figures/MAGeCK_RRA_depleted.pdf, MAGeCK_RRA_enriched.pdf,
#         MAGeCK_volcano.pdf, MAGeCK_volcano_significant_and_labelled.csv
# ======================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- Settings --------------------------------------------------------------
MAGECK_DIR   <- "input_files/20260910/mageck_outputs/FILTERING_T0_PLUS_CONTROL_PAIRED_TOTAL/"
GENE_SUMMARY <- file.path(MAGECK_DIR, "T0_PLUS_CONTROL_PAIRED_TOTAL.gene_summary.txt")
OUT_DIR      <- file.path(MAGECK_DIR, "figures")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

CONTRAST_LABEL <- "DRAXIN vs Ctrl"
N_LABEL     <- 10        # RRA plots: top genes per direction shown in the legend (MAGeCK default 10)
N_VOLCANO   <- 4         # volcano: top genes per direction labelled
FDR_CUT     <- 0.05      # volcano: "Significant" = FDR below this (set to the threshold used in the paper)
FIG_W_PT    <- 200       # width of each figure (pt; 72 pt = 1 inch)
FIG_H_PT    <- 200       # height (pt)
FIG_FONT    <- 8         # pt: titles, axis titles
FIG_SMALL   <- 7         # pt: tick labels, legend, gene labels
# MAGeCK's own colour order for the top genes (as in its default plots)
TOP_COLS <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF",
              "#999999", "#66C2A5", "#FC8D62")

pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")
save_fig <- function(p, name) {
  ggsave(file.path(OUT_DIR, name), p, width = FIG_W_PT / 72, height = FIG_H_PT / 72, device = pdf_device)
  cat("Saved:", name, "\n")
}
theme_fig <- theme_classic(base_size = FIG_FONT, base_family = "Helvetica") +
  theme(axis.text = element_text(size = FIG_SMALL, colour = "black"),
        axis.line = element_line(linewidth = 0.3), axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = FIG_FONT, hjust = 0.5, lineheight = 0.95),
        legend.text = element_text(size = FIG_SMALL), legend.title = element_text(size = FIG_SMALL),
        plot.margin = margin(4, 6, 4, 4))
pt <- 1 / ggplot2::.pt
# legend inside the panel, top right (syntax differs before / from ggplot2 3.5)
legend_inside <- if (packageVersion("ggplot2") >= "3.5.0")
  theme(legend.position = "inside", legend.position.inside = c(0.98, 0.98)) else
    theme(legend.position = c(0.98, 0.98))
gene_name <- function(x) gsub("_", ";", x)          # MAGeCK multi-gene arrays: CDH1_CDH3 -> CDH1;CDH3

# ---- Data ------------------------------------------------------------------
# check.names = FALSE keeps MAGeCK's column names (neg|score, pos|fdr, ...)
gs <- read.delim(GENE_SUMMARY, check.names = FALSE)
need <- c("id", "neg|score", "neg|fdr", "neg|rank", "neg|lfc", "pos|score", "pos|fdr", "pos|rank", "pos|lfc")
if (!all(need %in% names(gs))) stop("Missing columns in the gene summary: ", paste(setdiff(need, names(gs)), collapse = ", "))
cat("Genes in the gene summary:", nrow(gs), "\n")

# ---- 1-2. RRA score distributions (as MAGeCK's plot: ranked RRA scores, log
# scale with the best scores on top, top genes as coloured points + legend) --
rra_plot <- function(dir) {
  d <- gs %>% transmute(id = gene_name(id), score = .data[[paste0(dir, "|score")]],
                        rank = .data[[paste0(dir, "|rank")]]) %>%
    arrange(score) %>% mutate(order = row_number())
  top <- d %>% filter(rank <= N_LABEL) %>% arrange(rank) %>% mutate(id = factor(id, id))
  ggplot(d, aes(order, score)) +
    geom_line(linewidth = 0.3, colour = "black") +
    geom_point(data = top, aes(colour = id), size = 1.8) +
    scale_colour_manual(values = setNames(TOP_COLS[seq_len(nrow(top))], levels(top$id)), name = NULL) +
    scale_y_continuous(trans = c("log10", "reverse"), labels = scales::label_log()) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 4), expand = expansion(mult = c(0.03, 0.05))) +
    labs(x = "Genes", y = "RRA score", title = paste0("RRA scores\n(", CONTRAST_LABEL, ")")) +
    guides(colour = guide_legend(override.aes = list(size = 1.8), keyheight = unit(3, "mm"))) +
    theme_fig +
    legend_inside +
    theme(legend.justification = c(1, 1),
          legend.background = element_rect(colour = "black", linewidth = 0.3),
          legend.margin = margin(2, 3, 2, 2), legend.text = element_text(size = FIG_SMALL - 1, face = "italic"))
}
save_fig(rra_plot("pos"), "MAGeCK_RRA_enriched.pdf")   # enriched arrays (resistance arm)
save_fig(rra_plot("neg"), "MAGeCK_RRA_depleted.pdf")   # depleted arrays (sensitising arm)

# ---- 3. Volcano: -log10 FDR vs log2 fold change ----------------------------------
# FDR from the direction of the fold change (MAGeCK tests both tails separately)
vol <- gs %>%
  transmute(id = gene_name(id), lfc = `neg|lfc`,
            fdr = ifelse(`neg|lfc` < 0, `neg|fdr`, `pos|fdr`),
            significant = factor(ifelse(fdr < FDR_CUT, "yes", "no"), c("no", "yes")),
            top = (`neg|rank` <= N_VOLCANO & lfc < 0) | (`pos|rank` <= N_VOLCANO & lfc > 0)) %>%
  mutate(neglog_fdr = -log10(pmax(fdr, .Machine$double.xmin)))
cat("Volcano: FDR <", FDR_CUT, "-- significant depleted:", sum(vol$significant == "yes" & vol$lfc < 0),
    "| enriched:", sum(vol$significant == "yes" & vol$lfc > 0), "\n")
x_max <- max(abs(vol$lfc), na.rm = TRUE) * 1.05
y_top <- max(vol$neglog_fdr, na.rm = TRUE)
p_vol <- ggplot(vol %>% arrange(significant), aes(lfc, neglog_fdr, colour = significant)) +
  geom_vline(xintercept = 0, colour = "black", linewidth = 0.4) +
  geom_point(size = 1.1, stroke = 0) +
  geom_text_repel(data = filter(vol, top), aes(label = id), colour = "black", size = (FIG_SMALL - 1) * pt,
                  segment.size = 0.2, min.segment.length = 0.2, max.overlaps = Inf, box.padding = 0.2) +
  # direction arrows above the panel, as in the figure
  annotate("segment", x = -0.15 * x_max, xend = -0.95 * x_max, y = y_top * 1.28, yend = y_top * 1.28,
           arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.3) +
  annotate("segment", x = 0.15 * x_max, xend = 0.95 * x_max, y = y_top * 1.28, yend = y_top * 1.28,
           arrow = arrow(length = unit(1.2, "mm")), linewidth = 0.3) +
  annotate("text", x = -0.55 * x_max, y = y_top * 1.31, label = "depleted\n(sensitising arm)",
           vjust = 0, size = (FIG_SMALL - 1) * pt, lineheight = 0.9) +
  annotate("text", x = 0.55 * x_max, y = y_top * 1.31, label = "enriched\n(resistance arm)",
           vjust = 0, size = (FIG_SMALL - 1) * pt, lineheight = 0.9) +
  scale_colour_manual(values = c(no = "grey70", yes = "#E64B35"), name = "Significant", drop = FALSE) +
  scale_x_continuous(limits = c(-x_max, x_max), breaks = scales::breaks_pretty(n = 5)) +
  scale_y_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = c(0.02, 0.02))) +
  coord_cartesian(ylim = c(0, y_top * 1.5), clip = "off") +
  labs(x = paste0("log2 fold change\n(", CONTRAST_LABEL, ")"), y = "-log10(FDR)") +
  guides(colour = guide_legend(override.aes = list(size = 2), keyheight = unit(3, "mm"))) +
  theme_fig +
  theme(axis.line.x = element_blank(), legend.position = "right", legend.margin = margin(0, 0, 0, 0))
save_fig(p_vol, "MAGeCK_volcano.pdf")

write.csv(vol %>% filter(top | significant == "yes") %>% arrange(lfc),
          file.path(OUT_DIR, "MAGeCK_volcano_significant_and_labelled.csv"), row.names = FALSE)
cat("Done.\n")
