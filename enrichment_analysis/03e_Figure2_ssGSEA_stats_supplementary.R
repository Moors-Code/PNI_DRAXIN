# ============================================================
# Figure2_ssGSEA_stats_supplementary.R
#
# SUPPLEMENTARY / OPTIONAL -- not part of the paper's reported figures.
# Tests whether pathway activity (ssGSEA z-scores) differs across the
# three nerve categories using repeated-measures ANOVA, accounting for
# patient-level pairing. Kept for reference/reviewer questions.
#
# Input : processed/Figure2_GSEA/ssGSEA/nerve_scores_data.rds
# Output: output/Figure2_GSEA/tables/supplementary/ANOVA_pathway_differences.csv
# ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

processed_dir <- file.path("processed", "Figure2_GSEA", "ssGSEA")
tables_dir    <- file.path("output", "Figure2_GSEA", "tables", "supplementary")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

bundle_path <- file.path(processed_dir, "nerve_scores_data.rds")
stopifnot(file.exists(bundle_path))
bundle <- readRDS(bundle_path)
nerve_scores <- bundle$nerve_scores

# ------------------------------------------------------------
# patient-level aggregation + repeated-measures ANOVA per pathway
# ------------------------------------------------------------
patient_scores <- nerve_scores %>%
  group_by(patient_id, nerve_types, pathway) %>%
  summarise(mean_z = mean(z_score), .groups = "drop")

run_rm_anova <- function(df_one_pathway) {
  df_one_pathway <- droplevels(df_one_pathway)
  fit <- aov(mean_z ~ nerve_types + Error(patient_id / nerve_types), data = df_one_pathway)
  sm <- summary(fit)
  p <- sm[["Error: patient_id:nerve_types"]][[1]][["Pr(>F)"]][1]
  tibble::tibble(p_value = p)
}

anova_tbl <- patient_scores %>%
  group_by(pathway) %>%
  group_modify(~ run_rm_anova(.x)) %>%
  ungroup() %>%
  mutate(FDR = p.adjust(p_value, method = "BH"))

readr::write_csv(anova_tbl, file.path(tables_dir, "ANOVA_pathway_differences.csv"))

message("Done (supplementary). Table written to: ", tables_dir)
