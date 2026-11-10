# ============================================================
# 03d_Figure2_ssGSEA_plots.R
#
# READ-AND-PLOT ONLY: loads the nerve-level pathway score bundle
# written by 03c_Figure2_ssGSEA_compute.R and produces:
#   - the full pathway score table (all pathways, all nerves)
#   - Panel 1: complement/fibroblast/leukocyte pathways
#   - Panel 2: interferon/IL6-JAK-STAT3/prostaglandin pathways
#
# For the repeated-measures ANOVA across nerve types (not used in the
# paper), see Figure2_ssGSEA_stats_supplementary.R.
#
# Input : processed/Figure2_GSEA/ssGSEA/nerve_scores_data.rds
# Output: output/Figure2_GSEA/tables/nerve_pathway_scores.csv
#         output/Figure2_GSEA/figures/complememnt_rolling.pdf
#         output/Figure2_GSEA/figures/pathways_figure2.pdf
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(ggplot2)
})

# ------------------------------------------------------------
# paths + bundle
# ------------------------------------------------------------
processed_dir <- file.path("processed", "Figure2_GSEA", "ssGSEA")
tables_dir    <- file.path("output", "Figure2_GSEA", "tables")
figures_dir   <- file.path("output", "Figure2_GSEA", "figures")
dir.create(tables_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

bundle_path <- file.path(processed_dir, "nerve_scores_data.rds")
stopifnot(file.exists(bundle_path))
bundle <- readRDS(bundle_path)

nerve_scores <- bundle$nerve_scores

# ------------------------------------------------------------
# full, unfiltered pathway score table
# ------------------------------------------------------------
readr::write_csv(nerve_scores, file.path(tables_dir, "nerve_pathway_scores.csv"))

# ------------------------------------------------------------
# shared plot settings
# ------------------------------------------------------------
facet_scales <- "free_y"
base_size    <- 11
panel_aspect <- 0.5

make_pathway_panel <- function(nerve_scores, gene_sets_to_plot, facet_ncol) {
  df <- nerve_scores %>%
    filter(pathway %in% gene_sets_to_plot) %>%
    mutate(
      pathway = factor(pathway, levels = gene_sets_to_plot),
      nerve_types = factor(nerve_types)
    )

  ggplot(df, aes(x = nerve_types, y = z_score)) +
    geom_violin(fill = "grey85", color = NA, alpha = 0.6, trim = FALSE) +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white", color = "grey30", linewidth = 0.4) +
    geom_jitter(aes(color = patient_id), width = 0.10, size = 0.9, alpha = 0.95) +
    facet_wrap(~ pathway, scales = facet_scales, ncol = facet_ncol) +
    theme_classic(base_size = base_size) +
    theme(
      aspect.ratio = panel_aspect,
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(size = base_size * 0.85),
      legend.position = "right",
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
      panel.spacing = unit(6, "pt")
    ) +
    ylab("ssGSEA score (SD units, per nerve)") +
    xlab("Nerve type") +
    labs(color = "Patient ID")
}

# ------------------------------------------------------------
# Panel 1: complement / stromal / leukocyte pathways
# ------------------------------------------------------------
panel1_pathways <- c(
  "C5: GOBP_LEUKOCYTE_TETHERING_OR_ROLLING",
  "C2: REACTOME_INITIAL_TRIGGERING_OF_COMPLEMENT",
  "C4: GAVISH_FIBROBLASTS_COMPLEMENT"
)

p1 <- make_pathway_panel(nerve_scores, panel1_pathways, facet_ncol = 1)

ggsave(
  filename = file.path(figures_dir, "complememnt_rolling.pdf"),
  plot = p1, width = 5, height = 10, units = "in", dpi = 600
)

# ------------------------------------------------------------
# Panel 2: interferon / IL6-JAK-STAT3 / prostaglandin pathways
# ------------------------------------------------------------
panel2_pathways <- c(
  "H: IFNA_RESPONSE",
  "H: IFNG_RESPONSE",
  "H: IL6_JAK_STAT3",
  "C5: GOBP_RESPONSE_TO_PROSTAGLANDIN_E"
)

p2 <- make_pathway_panel(nerve_scores, panel2_pathways, facet_ncol = 2)

ggsave(
  filename = file.path(figures_dir, "pathways_figure2.pdf"),
  plot = p2, width = 10, height = 10, units = "in", dpi = 600
)

message("Done.")
message("Table written to:   ", tables_dir)
message("Figures written to: ", figures_dir)
