# ======================================================================
# 02_spatial_association_with_nerve_step3_prioritisation_plot.R
#
# Figure 3C: how the spatial test prioritises the DE genes.
# Every DE gene (InvN vs NI, FDR < 0.1, log2FC > 0) is placed by its
# association with the invaded-nerve mask (y) and with the interface ring
# (x), both as delta L = Lee's L of the real mask minus the mean of its 200
# relocated shams (mean over patients). Colour = result of step1:
#   core signature (significant and above the shams in every patient) or
#   not in the core signature.
# Reference genes (NOT used for selection, visual check only):
#   tumour programmes (top genes of the tumour NMF factors) and the nerve
#   programme (top genes of the nerve factor that was used to define the
#   nerve masks -- their nerve association is partly expected by
#   construction, so they are shown, not used as a positive control).
#
# INPUTS: tables/lee_sham_DE_genes_scores.csv (step1),
#         input/TableS2_NMF_top_genes.csv
# OUTPUTS: figures/fig3C_prioritisation.pdf, tables/fig3C_reference_and_labelled_genes.csv
# ======================================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
TABLES_DIR   <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
FIGURES_DIR  <- file.path(PROJECT_DIR, "output", "02_spatial_association", "figures")
dir.create(FIGURES_DIR, recursive = TRUE, showWarnings = FALSE)
# Publication PDFs: base pdf device, standard PostScript font (Helvetica) kept as
# editable text, no dingbat glyphs for points -- opens cleanly in Illustrator /
# Acrobat (cairo_pdf embeds Linux fonts such as DejaVu that Adobe does not know).
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

TUMOUR_FACTORS <- c("NMF_4", "NMF_13", "NMF_14", "NMF_19")   # epithelial / tumour programmes
NERVE_FACTORS  <- c("NMF_27")                       # nerve programme (defines the nerve masks)
TOP_N_FACTOR   <- 50                                # top genes per factor used as reference
PANEL_A_GENES  <- c("NES", "KRT19")   # the example genes of panel A (step2 GENES): always marked,
# so the reader can find them in this plot
N_LABEL        <- 5      # labels are chosen from the data, not by hand: the top N DE genes by nerve delta L

# ---- 1. Data ---------------------------------------------------------------
res <- read.csv(file.path(TABLES_DIR, "lee_sham_DE_genes_scores.csv"))   # DE_pass family only (step1)
if (!"signature_core" %in% names(res))                       # older step1 output
  res$signature_core <- res$family == "DE_pass" & res$nerve_status == "pass" & res$nerve_n_patients_exceed_q95 == 3
de <- res %>% filter(family == "DE_pass", nerve_status != "not_tested")

nmf <- read.csv(file.path(PROJECT_DIR, "input", "TableS2_NMF_top_genes.csv"))
top_of <- function(f) unique(unlist(nmf[nmf$component %in% f, 2:(TOP_N_FACTOR + 1)]))
tumour_ref <- setdiff(top_of(TUMOUR_FACTORS), "")
nerve_ref  <- setdiff(top_of(NERVE_FACTORS), "")
both <- intersect(tumour_ref, nerve_ref)
if (length(both)) { tumour_ref <- setdiff(tumour_ref, both); nerve_ref <- setdiff(nerve_ref, both) }

de <- de %>% mutate(
  tier = factor(case_when(signature_core ~ "Core signature (all patients)",
                          TRUE ~ "Not prioritised"),
                levels = c("Core signature (all patients)", "Not prioritised")),
  reference = case_when(Gene %in% tumour_ref ~ "Tumour programme (NMF)",
                        Gene %in% nerve_ref  ~ "Nerve programme (NMF, defines nerve masks)",
                        TRUE ~ NA_character_))

# Reference genes are top genes of whole-tissue NMF factors; only those that are
# also DE (InvN vs NI) can appear in the plot. Report how many that is.
ref_cov <- data.frame(reference = c("Tumour programme (NMF)", "Nerve programme (NMF, defines nerve masks)"),
                      n_top_genes = c(length(tumour_ref), length(nerve_ref)),
                      n_in_DE_plot = c(sum(tumour_ref %in% de$Gene), sum(nerve_ref %in% de$Gene)))
cat("NMF reference genes present among the plotted DE genes:\n"); print(ref_cov)

ref_tab <- de %>% filter(!is.na(reference)) %>% count(reference, tier)
cat("Reference genes among DE genes, by tier:\n"); print(ref_tab)

# ---- 1b. Genes to label (data-driven) ------------------------------------------
top_by <- function(d, col, why) d %>% filter(is.finite(.data[[col]])) %>%
  slice_max(.data[[col]], n = N_LABEL, with_ties = FALSE) %>% transmute(Gene, label_reason = why)
label_tab <- top_by(de, "nerve_delta_L_mean", "Top nerve delta L (all DE genes)")
cat("\nLabelled genes:", nrow(label_tab), "\n")
write.csv(de %>% filter(!is.na(reference) | Gene %in% label_tab$Gene) %>%
            left_join(label_tab, by = "Gene") %>%
            select(Gene, reference, label_reason, tier, nerve_delta_L_mean, interface_delta_L_mean,
                   nerve_padj, nerve_consistency),
          file.path(TABLES_DIR, "fig3C_reference_and_labelled_genes.csv"), row.names = FALSE)

# ---- 2. Plot ---------------------------------------------------------------
tier_col <- c("Core signature (all patients)" = "#2166AC",
              "Not prioritised" = "grey85")
ref_col  <- c("Tumour programme (NMF)" = "#D6604D", "Nerve programme (NMF, defines nerve masks)" = "#1A9850")
n_lab <- de %>% count(tier) %>% mutate(lab = paste0(tier, " (n = ", n, ")"))
lims <- range(c(de$nerve_delta_L_mean, de$interface_delta_L_mean), na.rm = TRUE)

p <- ggplot(de, aes(interface_delta_L_mean, nerve_delta_L_mean)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
  geom_point(data = filter(de, tier == "Not prioritised"), aes(fill = tier), shape = 21, colour = "transparent", size = 0.9) +
  geom_point(data = filter(de, tier != "Not prioritised"), aes(fill = tier), shape = 21, colour = "transparent", size = 1.5) +
  geom_point(data = filter(de, !is.na(reference)), aes(colour = reference), shape = 21, fill = NA,
             size = 2.4, stroke = 0.7) +
  geom_point(data = filter(de, Gene %in% PANEL_A_GENES), shape = 22, fill = NA, colour = "black", size = 3.2, stroke = 0.8) +
  geom_text_repel(data = filter(de, Gene %in% PANEL_A_GENES), aes(label = Gene), size = 3, fontface = "bold.italic",
                  nudge_x = 0.01, min.segment.length = 0, segment.size = 0.3) +
  geom_text_repel(data = filter(de, Gene %in% setdiff(label_tab$Gene, PANEL_A_GENES)), aes(label = Gene), size = 2.6, fontface = "italic",
                  min.segment.length = 0, segment.size = 0.2, box.padding = 0.35, max.overlaps = 40) +
  scale_fill_manual(values = tier_col, labels = setNames(n_lab$lab, n_lab$tier), name = "Spatial test (step1)") +
  scale_colour_manual(values = ref_col, name = "Reference (not used for selection)") +
  labs(x = expression("Association with the interface ring ("*Delta*"L vs relocated shams)"),
       y = expression("Association with the invaded nerve ("*Delta*"L vs relocated shams)"),
       title = "Spatial prioritisation of genes up in invaded nerves",
       caption = paste0("Each point: one DE gene (InvN vs NI). Delta L = Lee's L of the real mask minus the mean of 200 ",
                        "geometry-matched relocated masks, mean of 3 patients.\nDashed line: equal association with ",
                        "nerve and interface. Reference: top ", TOP_N_FACTOR, " genes of tumour (",
                        paste(TUMOUR_FACTORS, collapse = ", "), ") and nerve (", paste(NERVE_FACTORS, collapse = ", "),
                        ") NMF factors\n(",
                        ref_cov$n_in_DE_plot[1], " and ", ref_cov$n_in_DE_plot[2], " of them are DE and plotted).\n",
                        "Squares: example genes of panel A. Labels: top ", N_LABEL, " genes by nerve delta L.")) +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 3)),
         colour = guide_legend(ncol = 1, override.aes = list(size = 3))) +
  theme_classic(base_size = 10) +
  theme(legend.position = "right", plot.caption = element_text(hjust = 0, size = 7),
        plot.title = element_text(face = "bold"))
out <- file.path(FIGURES_DIR, "fig3C_prioritisation.pdf")
ggsave(out, p, width = 8.5, height = 5.8, device = pdf_device)
cat("Saved:", out, "\n")
