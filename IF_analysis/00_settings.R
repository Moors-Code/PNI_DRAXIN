###############################################################################
# 00_settings.R - shared by all scripts (paths, which input folder, patient colours)
#
# Folder structure (run every script from the project folder, or from scripts/):
#   input/Cancer_analysis/                    raw HALO/Horizon CSV exports
#   input/Cancer_analysis_with_background/    raw exports incl. CD3/CD45 background
#   input/Nerve_analysis/                     raw nerve exports
#   input/Nerve_analysis_with_control/        raw nerve exports incl. control channels
#   input/thresholds.csv                      per-block thresholds (Image_Name + *_threshold)
#   output/<input folder>/cleaned/...         written by 01_clean_cancer.R / 02_clean_nerve.R
#   output/Cancer_final/, output/Nerve_final/, output/Nerve_NESTIN_ITGB4/,
#   output/Representative_examples/           written by the analysis scripts
###############################################################################

if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input")) stop("Run the scripts from the project folder (the one containing input/).")
if (!exists("%||%", mode = "function")) `%||%` <- function(x, y) if (is.null(x)) y else x

# ---- Which cleaned data the analysis scripts use --------------------------------
CANCER_INPUT <- "Cancer_analysis_with_background"   # or "Cancer_analysis"
NERVE_INPUT  <- "Nerve_analysis"                    # or "Nerve_analysis_with_control"

THRESHOLD_FILE <- file.path("input", "thresholds.csv")
cleaned_file   <- function(folder) file.path("output", folder, "cleaned", "combined",
                                             "All_images_combined.csv")

# ---- Blocks whose name does not start with the patient number ------------------
# (the cleaning scripts take the text before "bl" as Patient_ID; "2777_block6" -> "2777_")
PATIENT_OF_BLOCK <- c("2777_block6" = "27777")

# ---- Patient colours (same as the QuPath/COMET figures) -------------------------
# Patients 1 and 2 (red, green) are not in these analyses; their colours are not used.
PATIENTS <- data.frame(
  Patient_ID = c("70865",     "56945",     "80830",     "23885",     "27777",     "60456"),
  Label      = c("Patient 3", "Patient 4", "Patient 5", "Patient 6", "Patient 7", "Patient 8"),
  Colour     = c("#1f77b4",   "#ff7f0e",   "#9467bd",   "#8c564b",   "#17becf",   "#e6ab02"),
  stringsAsFactors = FALSE)

# Main block per patient gets the patient colour; any further block of the same
# patient a darker shade of it (e.g. 80830 bl3, 2777_block6)
MAIN_BLOCK <- c("70865" = "70865 bl15", "56945" = "56945 bl20", "80830" = "80830 bl4",
                "23885" = "23885 bl5",  "27777" = "27777 bl9",  "60456" = "60456 bl5")

darker <- function(col, f = 0.6) {
  rgb(t(col2rgb(col) * f), maxColorValue = 255)
}

# Colours and legend labels for a set of blocks (Image_Name) with their Patient_ID
block_colours <- function(blocks, patients) {
  blocks <- as.character(blocks); patients <- as.character(patients)
  unknown <- setdiff(unique(patients), PATIENTS$Patient_ID)
  if (length(unknown)) stop("Patient(s) without a colour in 00_settings.R (PATIENTS): ",
                            paste(unknown, collapse = ", "))
  base <- PATIENTS$Colour[match(patients, PATIENTS$Patient_ID)]
  main <- blocks == MAIN_BLOCK[patients] | is.na(MAIN_BLOCK[patients])
  setNames(ifelse(main, base, darker(base)), blocks)
}
block_labels <- function(blocks, patients) {
  blocks <- as.character(blocks); patients <- as.character(patients)
  lab <- PATIENTS$Label[match(patients, PATIENTS$Patient_ID)]
  setNames(paste0(lab, " (", blocks, ")"), blocks)
}

# ---- Shared SuperPlot (same format as the QuPath tenascin / p75NTR / GAP43 figures) ----
# df    : one row per ROI / nerve with columns panel, group, value, block (Image_Name)
# stats : one row per panel with columns panel, p_label
# Small dots = ROIs / nerves (coloured by patient/block); large dots = block medians,
# connected between groups; black bar = median of block medians; bracket = p value.
# Export size: 4 x 4.3 inch for one panel, (3.2 x panels + 1.5) x 4.3 inch for several.
superplot_2group <- function(df, group_levels, group_labels, colours, labels, stats,
                             ylab, log_y = TRUE, ref_line = NULL, seed = 1) {
  library(ggplot2)
  set.seed(seed)
  df <- df[!is.na(df$value) & (!log_y | df$value > 0), ]
  df$group <- factor(df$group, levels = group_levels)
  df$panel <- factor(df$panel, levels = unique(df$panel))
  df$x <- as.numeric(df$group) + runif(nrow(df), -0.28, 0.28)
  blocks <- unique(df$block); blocks <- blocks[order(labels[blocks])]
  # legend: "Patient N" only; the block is added only if a patient has several blocks here
  short <- sub(" \\(.*$", "", labels[blocks])
  labels[blocks] <- ifelse(short %in% short[duplicated(short)], labels[blocks], short)

  med <- aggregate(value ~ block + group + panel, data = df, FUN = median)
  off <- setNames(seq(-0.14, 0.14, length.out = max(length(blocks), 2))[seq_along(blocks)], blocks)
  med$x <- as.numeric(med$group) + off[med$block]
  gmed <- aggregate(value ~ group + panel, data = med, FUN = median)
  gmed$x <- as.numeric(gmed$group)

  top <- tapply(df$value, df$panel, max); rng <- tapply(df$value, df$panel, function(v) diff(range(c(0, v))))
  br <- stats[stats$panel %in% levels(df$panel), ]
  br$panel <- factor(br$panel, levels = levels(df$panel))
  pk <- as.character(br$panel)
  br$y    <- if (log_y) top[pk] * 1.6  else top[pk] + 0.08 * rng[pk]
  br$tick <- if (log_y) br$y / 1.12   else br$y - 0.03 * rng[pk]
  br$ylab <- if (log_y) br$y * 1.18   else br$y + 0.05 * rng[pk]

  p <- ggplot()
  if (!is.null(ref_line)) p <- p + geom_hline(yintercept = ref_line, linetype = "dashed",
                                              colour = "grey60", linewidth = 0.3)
  p <- p +
    geom_point(data = df, aes(x, value, colour = block), size = 0.9, alpha = 0.35, stroke = 0) +
    geom_line(data = med, aes(x, value, colour = block, group = block), linewidth = 0.6) +
    geom_point(data = med, aes(x, value, fill = block), shape = 21, size = 3, stroke = 0.4) +
    geom_segment(data = gmed, aes(x = x - 0.3, xend = x + 0.3, y = value, yend = value), linewidth = 1) +
    geom_segment(data = br, aes(x = 1, xend = 2, y = y, yend = y), linewidth = 0.3) +
    geom_segment(data = br, aes(x = 1, xend = 1, y = y, yend = tick), linewidth = 0.3) +
    geom_segment(data = br, aes(x = 2, xend = 2, y = y, yend = tick), linewidth = 0.3) +
    geom_text(data = br, aes(x = 1.5, y = ylab, label = p_label), size = 3) +
    facet_wrap(~ panel, nrow = 1, scales = "free_y") +
    scale_colour_manual(values = colours[blocks], limits = blocks, labels = labels[blocks], name = "Patient") +
    scale_fill_manual(values = colours[blocks], limits = blocks, labels = labels[blocks], name = "Patient") +
    scale_x_continuous(breaks = 1:2, labels = group_labels[group_levels], limits = c(0.5, 2.5)) +
    labs(x = NULL, y = ylab) +
    theme_classic(base_size = 11) +
    theme(strip.background = element_blank(), strip.text = element_text(size = 11))
  if (log_y) p + scale_y_log10()
  else p + scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05)))
}

save_superplot <- function(p, file_base, n_panels, formats = c("pdf", "png")) {
  width <- if (n_panels == 1) 4 else 3.2 * n_panels + 1.5
  for (f in formats) ggplot2::ggsave(paste0(file_base, ".", f), p, width = width, height = 4.3,
                                     dpi = 300, bg = "white")
}

