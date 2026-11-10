# ======================================================================
# 02_spatial_association_with_nerve_step5_ORA.R
#
# Figure 3F: over-representation analysis (ORA) of the core signature.
#
# Query: the core signature genes (step1: significant combined test and real
# Lee's L above the 95th percentile of the relocated masks in all 3 patients).
#
# Backgrounds
#   MAIN figure  -- genes upregulated in InvN vs NI (DE_pass, FDR < 0.1,
#       log2FC > 0): every gene that entered the spatial test, including those
#       that failed it, i.e. the set the signature was selected from. Question:
#       which functions distinguish the reproducibly nerve-associated genes
#       within the invasion-upregulated response? (Signals that are up with
#       invasion in general, e.g. from adjacent tumour tissue, are in both sets
#       and cancel out.)
#   SUPPLEMENT (tables) -- all genes tested in the InvN vs NI DE analysis
#       (TableS3): functions of the signature relative to all nerve-expressed
#       genes.
#
# Figures (main background, all core genes):
#   1. g:Profiler Manhattan plot (gostplot / publish_gostplot): all significant
#      terms of GO:BP, GO:MF, Reactome and KEGG; labelled: the N_LABEL most
#      significant GO terms of size MIN_TERM_SIZE-MAX_TERM_SIZE
#      (LABEL_DRIVERS_ONLY = TRUE: only g:Profiler driver terms).
#   2. GO:BP dot plot: top N_DOT GO:BP terms by adjusted p (no redundancy
#      removal: top terms sharing genes reflect one programme, visible from
#      the gene lines).
#   3. Optional (MAKE_SIMPLIFIED = TRUE): redundancy-reduced dot plot --
#      clusterProfiler simplify() if installed, else g:Profiler driver terms.
#   In the dot plots, under each term: its N_LABEL_GENES genes with the
#   largest nerve delta L (step1), i.e. the most nerve-associated genes behind
#   the term (a fixed rule; full gene lists in the tables).
#
# Tables: all terms for both backgrounds, for the whole query list.
#
# INPUTS: tables/lee_sham_signature_core.csv, tables/lee_sham_DE_genes_scores.csv
#         (step1), input/TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv
# OUTPUTS: figures/fig3F_ORA_manhattan.pdf, fig3F_ORA_dotplot_GO_BP.pdf
#          (+ fig3F_ORA_dotplot_GO_BP_simplified.pdf if MAKE_SIMPLIFIED);
#          tables/fig3F_ORA_all_terms_bg_upregulated.csv (main),
#          tables/TableS_ORA_all_terms_bg_DE_tested.csv (supplement),
#          tables/fig3F_ORA_dotplot_terms_shown.csv,
#          tables/fig3F_enrichGO_BP.csv (only if clusterProfiler is used)
# ======================================================================

library(dplyr)
library(ggplot2)
library(gprofiler2)

select <- dplyr::select; filter <- dplyr::filter; mutate <- dplyr::mutate
arrange <- dplyr::arrange; count <- dplyr::count; rename <- dplyr::rename

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
TABLES_DIR    <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
FIGURES_DIR   <- file.path(PROJECT_DIR, "output", "02_spatial_association", "figures")
for (d in c(TABLES_DIR, FIGURES_DIR)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
# Publication PDFs: base pdf device, standard PostScript font (Helvetica) kept as
# editable text, no dingbat glyphs for points -- opens cleanly in Illustrator.
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

SOURCES         <- c("GO:BP", "GO:MF", "GO:CC", "REAC", "KEGG")   # queried (all in the tables)
PLOT_SOURCES    <- c("GO:BP", "GO:MF", "REAC", "KEGG")            # Manhattan plot (GO:CC: tables only)
MIN_TERM_SIZE   <- 10        # term size limits for labels / dot plots (clusterProfiler defaults)
MAX_TERM_SIZE   <- 500
N_LABEL         <- 10        # Manhattan plot: terms labelled
LABEL_DRIVERS_ONLY <- FALSE  # Manhattan labels: FALSE = N_LABEL most significant GO terms (default);
# TRUE = only g:Profiler driver terms
MAKE_SIMPLIFIED <- FALSE     # also draw the redundancy-reduced dot plot (driver terms / simplify())
N_DOT           <- 10        # dot plots: terms shown
N_LABEL_GENES   <- 3         # dot plots: genes printed under each term
# Dot plot layout, drawn at final print size (text sizes in pt as printed)
DOT_WIDTH       <- 4.2       # inches (about half a page width)
DOT_FONT        <- 8         # pt, term labels
DOT_WRAP        <- 30        # characters per line for term names
DOT_SHOW_GENES  <- TRUE      # second line under each term with its top-delta L genes
DOT_CAPTION     <- FALSE     # method details go in the figure legend instead
SIMPLIFY_CUTOFF <- 0.7       # clusterProfiler::simplify (Wang); package default
ALL_LABEL       <- "All core genes"

# ---- 1. Gene sets and backgrounds ------------------------------------------
core <- unique(trimws(read.csv(file.path(TABLES_DIR, "lee_sham_signature_core.csv"))$Gene))
lee  <- read.csv(file.path(TABLES_DIR, "lee_sham_DE_genes_scores.csv"))        # DE_pass family (step1)
dL   <- lee %>% select(Gene, nerve_delta_L_mean)
de_tested <- read.csv(file.path(INPUT_DIR, "TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv"), check.names = FALSE)
bg_tested <- unique(trimws(de_tested$Gene))
bg_up     <- intersect(unique(trimws(lee$Gene)), bg_tested)
stopifnot(all(core %in% bg_up))
BACKGROUNDS <- list(upregulated = bg_up, DE_tested = bg_tested)
BG_TEXT <- c(upregulated = "genes upregulated in InvN vs NI (all genes entering the spatial test)",
             DE_tested   = "genes tested in the InvN vs NI DE analysis")
cat("Core signature genes: ", length(core), "\nBackgrounds: upregulated = ", length(bg_up),
    ", DE-tested = ", length(bg_tested), "\n", sep = "")

# One query: the whole prioritised list, not split into subgroups.
queries <- setNames(list(core), ALL_LABEL)

# ---- 2. g:Profiler, all queries x both backgrounds -------------------------
# Queried fresh on every run (online database; version and date stored).
HAS_HIGHLIGHT <- "highlight" %in% names(formals(gost))
if (!HAS_HIGHLIGHT) warning("gprofiler2 < 0.2.2: no driver-term highlighting. Update with install.packages('gprofiler2').")
run_gost <- function(genes, bg) {
  args <- list(query = genes, organism = "hsapiens", custom_bg = bg, correction_method = "fdr",
               evcodes = TRUE, sources = SOURCES)
  if (HAS_HIGHLIGHT) args$highlight <- TRUE
  do.call(gost, args)
}
gost_main <- NULL
ora <- bind_rows(lapply(names(BACKGROUNDS), function(b) bind_rows(lapply(names(queries), function(q) {
  r <- run_gost(queries[[q]], BACKGROUNDS[[b]])
  if (b == "upregulated" && q == ALL_LABEL) gost_main <<- r
  if (is.null(r)) return(NULL)
  r$result %>% mutate(query_set = q, n_genes_in_query_set = length(queries[[q]]), background = b,
                      background_size = length(BACKGROUNDS[[b]]))
}))))
gp_version <- tryCatch(get_version_info()$gprofiler_version, error = function(e) NA_character_)
if (nrow(ora) == 0) stop("No significant terms for any query.")
if (!"highlighted" %in% names(ora)) ora$highlighted <- NA
ora <- ora %>% mutate(fold_enrichment = (intersection_size / query_size) / (term_size / effective_domain_size))
cat("\nSignificant terms per background x query:\n")
print(count(ora, background, query_set, name = "significant_terms") %>% as.data.frame())
if (is.null(gost_main)) stop("No significant terms for all core genes against the upregulated background.")

# ---- 3. Manhattan plot (main background, all core genes) -------------------
g_plot <- gost_main; g_plot$result <- g_plot$result %>% filter(source %in% PLOT_SOURCES)
lab_pool <- g_plot$result %>% filter(grepl("^GO:", source), term_size >= MIN_TERM_SIZE, term_size <= MAX_TERM_SIZE)
if (LABEL_DRIVERS_ONLY && HAS_HIGHLIGHT) lab_pool <- lab_pool %>% filter(highlighted %in% TRUE)
label_ids <- lab_pool %>% arrange(p_value) %>% slice_head(n = N_LABEL) %>% pull(term_id)
man <- publish_gostplot(gostplot(g_plot, capped = FALSE, interactive = FALSE), highlight_terms = label_ids)
out_man <- file.path(FIGURES_DIR, "fig3F_ORA_manhattan.pdf")
h_man <- 4 + 0.25 * length(label_ids)
if (inherits(man, "ggplot")) ggsave(out_man, man, width = 9, height = h_man, device = pdf_device) else {
  pdf_device(out_man, width = 9, height = h_man); grid::grid.draw(man); dev.off() }
cat("Saved: fig3F_ORA_manhattan.pdf (labelled: ", length(label_ids),
    if (LABEL_DRIVERS_ONLY) " driver terms" else " most significant GO terms", ")\n", sep = "")

# ---- 4. GO:BP dot plots (main background, all core genes) ------------------
top_genes <- function(g) {
  g <- data.frame(Gene = strsplit(g, ",")[[1]]) %>% left_join(dL, by = "Gene") %>% arrange(desc(nerve_delta_L_mean))
  paste0(paste(head(g$Gene, N_LABEL_GENES), collapse = ", "), if (nrow(g) > N_LABEL_GENES) ", ..." else "")
}
bp_main <- gost_main$result %>% filter(source == "GO:BP", term_size >= MIN_TERM_SIZE, term_size <= MAX_TERM_SIZE)
as_dot <- function(d) d %>% transmute(term = term_name, term_id, p_adj = p_value, count = intersection_size,
                                      gene_ratio = intersection_size / query_size, genes = intersection)

# 4a. not simplified: top N_DOT GO:BP terms by adjusted p
dot_full <- bp_main %>% arrange(p_value) %>% slice_head(n = N_DOT) %>% as_dot()
full_method <- "g:Profiler (GO:BP), no redundancy removal"

# 4b. simplified (optional, MAKE_SIMPLIFIED)
USE_CP <- MAKE_SIMPLIFIED && all(vapply(c("clusterProfiler", "org.Hs.eg.db"), requireNamespace, logical(1), quietly = TRUE))
if (USE_CP) {
  suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
  ego <- enrichGO(gene = core, universe = bg_up, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 1,
                  minGSSize = MIN_TERM_SIZE, maxGSSize = MAX_TERM_SIZE)
  ego_df <- if (is.null(ego)) data.frame() else as.data.frame(ego)
  if (nrow(ego_df) > 0) {
    ego_s <- clusterProfiler::simplify(ego, cutoff = SIMPLIFY_CUTOFF, by = "p.adjust", select_fun = min, measure = "Wang")
    write.csv(ego_df %>% mutate(kept_after_simplify = ID %in% as.data.frame(ego_s)$ID),
              file.path(TABLES_DIR, "fig3F_enrichGO_BP.csv"), row.names = FALSE)
    dot_simpl <- as.data.frame(ego_s) %>% arrange(p.adjust) %>% slice_head(n = N_DOT) %>%
      transmute(term = Description, term_id = ID, p_adj = p.adjust, count = Count,
                gene_ratio = Count / as.numeric(sub(".*/", "", GeneRatio)), genes = gsub("/", ",", geneID))
  } else dot_simpl <- data.frame()
  simpl_method <- paste0("clusterProfiler enrichGO (GO:BP), redundant terms removed with simplify() (Wang, cutoff ",
                         SIMPLIFY_CUTOFF, ")")
} else if (MAKE_SIMPLIFIED) {
  message("clusterProfiler / org.Hs.eg.db not installed: simplified dot plot uses g:Profiler driver terms.")
  dot_simpl <- (if (HAS_HIGHLIGHT) filter(bp_main, highlighted %in% TRUE) else bp_main) %>%
    arrange(p_value) %>% slice_head(n = N_DOT) %>% as_dot()
  simpl_method <- paste0("g:Profiler (GO:BP), ", if (HAS_HIGHLIGHT) "redundancy removed by driver-term highlighting"
                         else "no redundancy removal available (update gprofiler2)")
} else {
  dot_simpl <- data.frame(); simpl_method <- NA
}

dot_plot <- function(d, method, file) {
  if (nrow(d) == 0) { cat("No GO:BP terms for", file, "\n"); return(invisible(d)) }
  d <- d %>% mutate(top_genes = vapply(genes, top_genes, character(1)),
                    label = if (DOT_SHOW_GENES) paste0(stringr::str_wrap(term, DOT_WRAP), "\n", top_genes)
                    else stringr::str_wrap(term, DOT_WRAP))
  p <- ggplot(d %>% mutate(label = reorder(label, gene_ratio)), aes(gene_ratio, label, size = count, colour = p_adj)) +
    geom_point() +
    scale_colour_gradient(low = "#e31a1c", high = "#2c7fb8", name = "Adj. p",
                          breaks = scales::breaks_pretty(n = 3), labels = function(x) format(x, digits = 1)) +
    scale_size_continuous(range = c(1.5, 5), name = "Genes", breaks = scales::breaks_pretty(n = 3)) +
    scale_x_continuous(breaks = scales::breaks_pretty(n = 3), expand = expansion(mult = 0.12)) +
    labs(x = "Gene ratio", y = NULL,
         caption = if (DOT_CAPTION) paste0(length(core), " core genes; background = ", length(bg_up),
                                           " DE-up genes; BH; ", method, "; top ", N_DOT, " by adjusted p.\n",
                                           "Genes: the ", N_LABEL_GENES, " term genes with the largest nerve delta L.")) +
    theme_bw(base_size = DOT_FONT) +
    theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_line(linewidth = 0.25),
          axis.text.y = element_text(size = DOT_FONT, colour = "black", lineheight = 0.85),
          axis.text.x = element_text(size = DOT_FONT - 1, colour = "black"),
          legend.position = "right", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0),
          legend.title = element_text(size = DOT_FONT - 1), legend.text = element_text(size = DOT_FONT - 2),
          legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.35, "cm"),
          plot.caption = element_text(hjust = 0, size = DOT_FONT - 2), plot.margin = margin(4, 6, 4, 4))
  h <- 0.8 + nrow(d) * (if (DOT_SHOW_GENES) 0.33 else 0.24)
  ggsave(file.path(FIGURES_DIR, file), p, width = DOT_WIDTH, height = h, device = pdf_device)
  cat("Saved:", file, "(", DOT_WIDTH, "x", round(h, 1), "in )\n")
  invisible(d)
}
shown_full  <- dot_plot(dot_full,  full_method,  "fig3F_ORA_dotplot_GO_BP.pdf")
shown_simpl <- if (MAKE_SIMPLIFIED) dot_plot(dot_simpl, simpl_method, "fig3F_ORA_dotplot_GO_BP_simplified.pdf") else data.frame()

# ---- 5. Tables ---------------------------------------------------------------
write_terms <- function(b, file) {
  write.csv(ora %>% filter(background == b) %>%
              mutate(parents = vapply(parents, paste, character(1), collapse = ";"),
                     labelled_in_manhattan = b == "upregulated" & query_set == ALL_LABEL & term_id %in% label_ids,
                     in_dotplot = b == "upregulated" & query_set == ALL_LABEL & term_id %in% shown_full$term_id,
                     in_dotplot_simplified = b == "upregulated" & query_set == ALL_LABEL & term_id %in% shown_simpl$term_id,
                     gprofiler_version = gp_version, query_date = as.character(Sys.Date())) %>%
              rename(driving_genes = intersection) %>% arrange(query_set, p_value),
            file.path(TABLES_DIR, file), row.names = FALSE)
  cat("Saved:", file, "\n")
}
write_terms("upregulated", "fig3F_ORA_all_terms_bg_upregulated.csv")
write_terms("DE_tested",   "TableS_ORA_all_terms_bg_DE_tested.csv")
write.csv(bind_rows(mutate(shown_full, plot = "not simplified"),
                    if (nrow(shown_simpl)) mutate(shown_simpl, plot = "simplified")),
          file.path(TABLES_DIR, "fig3F_ORA_dotplot_terms_shown.csv"), row.names = FALSE)
cat("\nDone.\n")
