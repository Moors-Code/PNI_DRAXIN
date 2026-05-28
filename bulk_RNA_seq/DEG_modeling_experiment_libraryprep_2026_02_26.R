library(edgeR)
library(stringr)
library(dplyr)
library(ggplot2)
library(tibble)
library(limma)
library(readr)

outdir <- "analysis_experiement_library_prepmodeled_2026_02_26"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# =========================================================
# 1. Read raw deduplicated count matrix
# =========================================================
count_data <- read.csv("data/counts_dedup.csv", row.names = 1, check.names = FALSE)

# =========================================================
# 2. Keep only ctrl and DRAXIN-related samples
# =========================================================
dataset_sub <- count_data[, startsWith(colnames(count_data), "ctrl") |
                            startsWith(colnames(count_data), "DRAXIN"),
                          drop = FALSE]

samples <- colnames(dataset_sub)

# =========================================================
# 3. Build metadata from sample names
# =========================================================
meta <- data.frame(
  sample = samples,
  condition = case_when(
    str_detect(samples, "DRAXIN") & str_detect(samples, "diluted|1_2000") ~ "diluted_DRAXIN",
    str_detect(samples, "DRAXIN") ~ "DRAXIN",
    str_detect(samples, "^ctrl") ~ "ctrl",
    TRUE ~ NA_character_
  ),
  experiment = case_when(
    str_detect(samples, "18_04_2024") ~ "exp1",
    str_detect(samples, "6_06_2024")  ~ "exp2",
    str_detect(samples, "2_05_2024")  ~ "exp3",
    TRUE ~ NA_character_
  ),
  stringsAsFactors = FALSE
)

# =========================================================
# 4. Restrict to experiments 1 and 2
# =========================================================
meta2 <- meta %>%
  filter(experiment %in% c("exp1", "exp2"),
         condition %in% c("ctrl", "DRAXIN", "diluted_DRAXIN")) %>%
  filter(!is.na(condition), !is.na(experiment))

# =========================================================
# 5. Define batch variable
# =========================================================
meta2$batch3 <- case_when(
  meta2$sample %in% c("ctrl_rep1_18_04_2024",
                      "DRAXIN_rep1_18_04_2024") ~ "exp1_altprep",
  meta2$experiment == "exp1" ~ "exp1_main",
  meta2$experiment == "exp2" ~ "exp2_main"
)

meta2$condition <- factor(meta2$condition,
                          levels = c("ctrl", "DRAXIN", "diluted_DRAXIN"))

meta2$batch3 <- factor(meta2$batch3,
                       levels = c("exp1_main", "exp1_altprep", "exp2_main"))

# =========================================================
# 6. Align counts to metadata
# =========================================================
counts2 <- dataset_sub[, meta2$sample, drop = FALSE]
stopifnot(identical(colnames(counts2), meta2$sample))

# =========================================================
# 7. Build design matrix
# =========================================================
design <- model.matrix(~ batch3 + condition, data = meta2)
print(colnames(design))

# =========================================================
# 8. Create edgeR object and filter lowly expressed genes
# =========================================================
#filtered genes
y <- DGEList(counts = counts2)
keep <- filterByExpr(y, group = meta2$condition, min.count = 5)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# all genes
y_all <- DGEList(counts = counts2)
y_all <- calcNormFactors(y_all)
expr_logCPM_TMM_all <- cpm(y_all, log = TRUE, prior.count = 1)


# =========================================================
# 8b. Export heatmap input files
# =========================================================
heatmap_input_dir <- file.path(outdir, "EDGE_RNAseq_GSEA", "Heatmaps_input")
dir.create(heatmap_input_dir, recursive = TRUE, showWarnings = FALSE)

expr_logCPM_TMM <- cpm(y, log = TRUE, prior.count = 1)

meta2_export <- meta2 %>%
  as.data.frame()

stopifnot(identical(colnames(expr_logCPM_TMM), meta2_export$sample))

saveRDS(
  expr_logCPM_TMM,
  file.path(heatmap_input_dir, "expr_logCPM_TMM.rds")
)

saveRDS(
  expr_logCPM_TMM_all,
  file.path(heatmap_input_dir, "expr_logCPM_TMM_all.rds")
)

readr::write_csv(
  meta2_export,
  file.path(heatmap_input_dir, "meta2.csv")
)

# =========================================================
# 9. Estimate dispersion and fit quasi-likelihood model
# =========================================================
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)

# =========================================================
# 10. Differential expression testing
# =========================================================

# DRAXIN vs ctrl
res_draxin <- glmQLFTest(fit, coef = "conditionDRAXIN")
tab_draxin <- topTags(res_draxin, n = Inf)$table

# diluted DRAXIN vs ctrl
res_dil <- glmQLFTest(fit, coef = "conditiondiluted_DRAXIN")
tab_dil <- topTags(res_dil, n = Inf)$table

# DRAXIN vs diluted DRAXIN
res_high_vs_low <- glmQLFTest(
  fit,
  contrast = makeContrasts(conditionDRAXIN - conditiondiluted_DRAXIN, levels = design)
)
tab_high_vs_low <- topTags(res_high_vs_low, n = Inf)$table

# =========================================================
# 11. Save result tables
# =========================================================
write.csv(tab_draxin, file.path(outdir, "edgeR_expBlocked_DRAXIN_vs_ctrl_ALL.csv"))
write.csv(tab_dil, file.path(outdir, "edgeR_expBlocked_diluted_vs_ctrl_ALL.csv"))
write.csv(tab_high_vs_low, file.path(outdir, "edgeR_expBlocked_DRAXIN_vs_diluted_ALL.csv"))

# =========================
# Volcano plots
# =========================
make_volcano <- function(tab, title,
                         fdr_cutoff = 0.05,
                         lfc_cutoff = 1,
                         label_top_n = 15) {
  
  df <- tab %>%
    tibble::rownames_to_column("gene") %>%
    mutate(
      FDR = ifelse(is.na(FDR), 1, FDR),
      logFC = ifelse(is.na(logFC), 0, logFC),
      neglog10FDR = -log10(pmax(FDR, 1e-300)),
      sig = (FDR < fdr_cutoff) & (abs(logFC) >= lfc_cutoff),
      direction = case_when(
        sig & logFC > 0 ~ "Up",
        sig & logFC < 0 ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  # label top significant genes only
  lab <- df %>%
    filter(sig) %>%
    arrange(FDR) %>%
    head(label_top_n)
  
  n_up <- sum(df$direction == "Up", na.rm = TRUE)
  n_down <- sum(df$direction == "Down", na.rm = TRUE)
  xmax <- max(abs(df$logFC), na.rm = TRUE)
  
  ggplot(df, aes(x = logFC, y = neglog10FDR)) +
    geom_point(aes(color = direction), size = 1.8, alpha = 0.8) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    geom_hline(yintercept = -log10(fdr_cutoff),
               linetype = "dashed", linewidth = 0.4) +
    geom_text(data = lab,
              aes(label = gene),
              vjust = -0.5,
              size = 3,
              check_overlap = TRUE) +
    coord_cartesian(xlim = c(-xmax, xmax)) +
    scale_color_manual(values = c("Up" = "firebrick", "Down" = "steelblue", "NS" = "grey70")) +
    theme_classic(base_size = 12) +
    labs(
      title = title,
      subtitle = paste0(
        "Up: ", n_up, " | Down: ", n_down,
        "  (FDR < ", fdr_cutoff, ", |log2FC| >= ", lfc_cutoff, ")"
      ),
      x = "log2 fold-change",
      y = "-log10(FDR)",
      color = NULL
    ) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
}

p1 <- make_volcano(
  tab_draxin,
  "High DRAXIN vs ctrl",
  lfc_cutoff = 0.3,
  fdr_cutoff = 0.2,
  label_top_n = 15
)

p2 <- make_volcano(
  tab_dil,
  "Diluted DRAXIN vs ctrl",
  lfc_cutoff = 0.3,
  fdr_cutoff = 0.2,
  label_top_n = 15
)

p3 <- make_volcano(
  tab_high_vs_low,
  "High DRAXIN vs Diluted DRAXIN",
  lfc_cutoff = 0.3,
  fdr_cutoff = 0.2,
  label_top_n = 15
)

print(p1)
print(p2)
print(p3)

ggsave(file.path(outdir, "volcano_DRAXIN_vs_ctrl.pdf"),
       plot = p1, width = 7, height = 6)

ggsave(file.path(outdir, "volcano_diluted_vs_ctrl.pdf"),
       plot = p2, width = 7, height = 6)

ggsave(file.path(outdir, "volcano_DRAXIN_vs_diluted.pdf"),
       plot = p3, width = 7, height = 6)

# =========================
# PCA of samples
# =========================

logcpm <- cpm(y, log = TRUE, prior.count =1)

pca <- prcomp(t(logcpm), scale. = FALSE)

pca_df <- data.frame(
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  condition = meta2$condition,
  experiment = meta2$experiment,
  sample = colnames(logcpm)
)

p_pca <- ggplot(pca_df, aes(PC1, PC2,
                            color = condition,
                            shape = experiment,
                            label = sample)) +
  geom_point(size = 4) +
  geom_text(vjust = -0.7, size = 3, check_overlap = TRUE) +
  theme_classic(base_size = 12) +
  labs(
    title = "PCA of samples",
    x = paste0("PC1 (", round(100 * summary(pca)$importance[2, 1], 1), "%)"),
    y = paste0("PC2 (", round(100 * summary(pca)$importance[2, 2], 1), "%)")
  )

print(p_pca)

ggsave(file.path(outdir, "PCA_samples.pdf"),
       plot = p_pca, width = 7, height = 6)
