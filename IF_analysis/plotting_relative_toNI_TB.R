###############################################################################
# IF nerve markers -- one figure, one scale: intensity relative to NI-TB
#
# Combines the per-nerve outputs of the two IF analyses and plots every marker
# as intensity RELATIVE TO THE NI-TB NERVES OF THE SAME SECTION:
#   relative = nerve value / geometric mean of the NI-TB nerves in that section
#   -> NI-TB centres on 1, InvN shows the within-section fold difference.
# A ratio is always > 0 (no negative values); it is drawn on a log2 axis so that
# 0.5x and 2x are equally far from 1, with ratio tick labels (0.5, 1, 2).
# The geometric mean matches the statistics, which model LOG intensity: the plot
# shows the same within-section fold difference the mixed models estimate.
# Within one section the background normalisation of markers_stats.R is a single
# factor, so it cancels: raw and normalised values give the same ratio.
#
# Nothing is re-analysed. P values are read from the two statistics outputs
# and Holm-adjusted together over all panels in this figure.
#
# Inputs (set NESTIN_DIR below):
#   NESTIN_DIR/Per_nerve_values.csv   and  NESTIN_DIR/Statistics.csv
#       (05_nerve_NESTIN_ITGB4_analysis.R; measure = median cell intensity)
#   MARKERS_DIR/markers_per_nerve.csv and  MARKERS_DIR/markers_model_results.csv
#       (markers_stats.R; GAP43, p75NTR, tenascin)
# Outputs (OUT_DIR):
#   IF_markers_relative_to_NITB_per_nerve.csv   one row per nerve and marker
#   IF_markers_relative_to_NITB_pvalues.csv     raw and Holm-adjusted p per panel
#   Fig_IF_markers_relative_to_NITB.pdf/.png    all panels side by side
###############################################################################
library(ggplot2)

# ---- Settings -----------------------------------------------------------------
# Outputs of 05_nerve_NESTIN_ITGB4_analysis.R, which lives in the Horizon
# analysis project: this script reads Per_nerve_values.csv (per-nerve median
# NESTIN and ITGB4 intensities) and Statistics.csv (fold changes and p values
# for the "median cell intensity" measure). Set this path before running.
NESTIN_DIR  <- "<path to Horizon_analysis/output/Nerve_NESTIN_ITGB4>"

MARKERS_DIR <- "output/"                      # markers_per_nerve.csv, markers_model_results.csv
OUT_DIR     <- "output/IF_markers_relative_to_NITB"          # where the combined figure goes
MIN_CELLS   <- 10                            # as in 05_nerve_NESTIN_ITGB4_analysis.R

# Panels, in plotting order. source: "nestin" (05 script) or "markers" (markers_stats.R);
# key: marker name in Statistics.csv, or analysis id in markers_model_results.csv.
PANELS <- data.frame(
  label  = c("ITGB4",  "Nestin", "GAP43",       "p75NTR, nerve", "p75NTR, perineural ring", "Tenascin-C, 0-50 um ring"),
  source = c("nestin", "nestin", "markers",     "markers",       "markers",                 "markers"),
  key    = c("ITGB4",  "NESTIN", "gap43_nerve", "p75ntr_nerve",  "p75ntr_perineural", "tenascin_immune_ring"),
  stringsAsFactors = FALSE)

# Patient labels and colours (same as the other IF figures). Patients 1 and 2
# (red, green) are not in these analyses; their colours are not used.
PATIENTS <- data.frame(
  Patient_ID = c("70865",     "56945",     "80830",     "23885",     "27777",     "60456"),
  Label      = c("Patient 3", "Patient 4", "Patient 5", "Patient 6", "Patient 7", "Patient 8"),
  Colour     = c("#1f77b4",   "#ff7f0e",   "#9467bd",   "#8c564b",   "#17becf",   "#e6ab02"),
  stringsAsFactors = FALSE)
# Main block per patient gets the patient colour; any further block of the same
# patient a darker shade of it (e.g. 80830 bl3, 2777_block6)
MAIN_BLOCK <- c("70865" = "70865 bl15", "56945" = "56945 bl20", "80830" = "80830 bl4",
                "23885" = "23885 bl5",  "27777" = "27777 bl9",  "60456" = "60456 bl5")
darker <- function(col, f = 0.6) rgb(t(col2rgb(col) * f), maxColorValue = 255)
GROUPS <- c("NI-TB", "InvN")

# Layout of the combined figure -- edit freely:
#   PANEL_ORDER: which panels to show and in which order (labels as in PANELS$label;
#                panels left out here appear only as single files)
#   GRID_NCOL / GRID_NROW: grid of the combined figure (e.g. 3 x 2, 6 x 1); set one to
#                NULL to let it follow from the other
PANEL_ORDER <- c("ITGB4", "Nestin", "p75NTR, nerve", "p75NTR, perineural ring", "GAP43", "Tenascin-C, 0-50 um ring")
GRID_NCOL   <- 2
GRID_NROW   <- 3
FIG_PANEL_W_MM <- 38      # width per panel (single and combined figures)
FIG_PANEL_H_MM <- 55      # height per panel row
FONT_PT        <- 7
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
pdf_device <- function(filename, width, height, ...)
  grDevices::pdf(file = filename, width = width, height = height, useDingbats = FALSE,
                 family = "Helvetica", encoding = "ISOLatin1.enc")

# ---- 1. Per-nerve values from both analyses -----------------------------------
read_req <- function(f) { if (!file.exists(f)) stop("Missing input: ", f); read.csv(f, check.names = FALSE) }

nest <- read_req(file.path(NESTIN_DIR, "Per_nerve_values.csv"))
nest <- nest[nest$N_cells >= MIN_CELLS, ]
nest$group <- ifelse(nest$Invasion == "Invaded", "InvN", ifelse(nest$Invasion == "Uninvaded", "NI-TB", NA))
nest_long <- do.call(rbind, lapply(PANELS$key[PANELS$source == "nestin"], function(m) {
  col <- paste0("Median_", m, "_intensity")
  if (!col %in% names(nest)) stop("Column ", col, " not in Per_nerve_values.csv")
  data.frame(source = "nestin", key = m, section = nest$Image_Name, nerve = paste(nest$Image_Name, nest$ROI),
             patient = as.character(nest$Patient_ID), group = nest$group, value = nest[[col]])
}))

mk <- read_req(file.path(MARKERS_DIR, "markers_per_nerve.csv"))
mk_long <- data.frame(source = "markers", key = mk$analysis, section = mk$image_name, nerve = paste(mk$image_name, mk$nerve_id),
                      patient = as.character(mk$patient), group = mk$group, value = mk$normalised)
mk_long <- mk_long[mk_long$key %in% PANELS$key[PANELS$source == "markers"], ]

d <- rbind(nest_long, mk_long)
d <- d[d$group %in% GROUPS & is.finite(d$value) & d$value > 0, ]

# Only sections with both nerve types (the comparison is within sections, as in the statistics)
both <- with(d, tapply(group, paste(key, section), function(g) all(GROUPS %in% g)))
d <- d[both[paste(d$key, d$section)], ]

# ---- 2. Relative to the NI-TB nerves of the same section ----------------------
ref <- with(d[d$group == "NI-TB", ], tapply(log(value), paste(key, section), mean))   # log geometric mean
d$relative <- exp(log(d$value) - ref[paste(d$key, d$section)])
d$panel    <- factor(PANELS$label[match(d$key, PANELS$key)], levels = PANELS$label)
d$group    <- factor(d$group, levels = GROUPS)
unknown <- setdiff(unique(d$patient), PATIENTS$Patient_ID)
if (length(unknown)) stop("Patient ID(s) not in PATIENTS: ", paste(unknown, collapse = ", "))
d$patient_label <- PATIENTS$Label[match(d$patient, PATIENTS$Patient_ID)]

# Block colours and labels, decided WITHIN each analysis (the two analyses may spell
# block names differently, so a block is identified by analysis + name). The main
# block of a patient (MAIN_BLOCK; names compared ignoring case, spaces, "_" and "-")
# gets the patient colour and the label "Patient N"; any further block of the same
# patient in that analysis a darker shade (2nd x0.6, 3rd x0.6^2) and "Patient N (<block>)".
# A patient with a single block in an analysis always gets the patient colour.
norm_name <- function(x) tolower(gsub("[ _-]", "", x))
blocks <- unique(d[, c("source", "patient", "patient_label", "section")])
blocks <- blocks[order(blocks$source, blocks$patient, blocks$section), ]
blocks$colour <- NA_character_; blocks$block_label <- NA_character_
for (sc in unique(blocks$source)) for (p in unique(blocks$patient[blocks$source == sc])) {
  i <- which(blocks$source == sc & blocks$patient == p)
  main <- which(norm_name(blocks$section[i]) == norm_name(MAIN_BLOCK[p]))
  if (!length(main)) {
    if (length(i) > 1) message(sc, ", patient ", p, ": main block '", MAIN_BLOCK[p], "' not found among ",
                               paste(blocks$section[i], collapse = ", "), " -> first block used as main")
    main <- 1
  }
  ord  <- c(i[main[1]], i[-main[1]])
  base <- PATIENTS$Colour[PATIENTS$Patient_ID == p]
  for (k in seq_along(ord)) {
    blocks$colour[ord[k]]      <- if (k == 1) base else darker(base, 0.6^(k - 1))
    blocks$block_label[ord[k]] <- if (k == 1) blocks$patient_label[ord[k]] else
      paste0(blocks$patient_label[ord[k]], " (", sub(paste0("^", p, "[ _-]*"), "", blocks$section[ord[k]]), ")")
  }
}
lv <- unique(blocks[order(match(blocks$patient, PATIENTS$Patient_ID), blocks$block_label), c("block_label", "colour")])
block_cols <- setNames(lv$colour, lv$block_label)
d$block_label <- factor(blocks$block_label[match(paste(d$source, d$section), paste(blocks$source, blocks$section))],
                        levels = names(block_cols))
write.csv(d[, c("panel", "patient_label", "section", "block_label", "nerve", "group", "value", "relative")],
          file.path(OUT_DIR, "IF_markers_relative_to_NITB_per_nerve.csv"), row.names = FALSE)

cat("Nerves per panel and group:\n"); print(table(d$panel, d$group))
missing <- setdiff(PANELS$label, unique(as.character(d$panel)))
if (length(missing)) message("No data for: ", paste(missing, collapse = ", "))

# Section (block) medians (large dots) and the median of block medians (black bar)
sec <- aggregate(relative ~ panel + source + section + block_label + group, data = d, FUN = median)
bar <- aggregate(relative ~ panel + group, data = sec, FUN = median)

# ---- 3. P values from the two statistics outputs, Holm over all panels --------
st <- read_req(file.path(NESTIN_DIR, "Statistics.csv"))
st <- st[st$Measure == "median cell intensity", ]
ms <- read_req(file.path(MARKERS_DIR, "markers_model_results.csv"))
PANELS$p <- ifelse(PANELS$source == "nestin",
                   st$p[match(PANELS$key, st$Marker)],
                   ms$p_value[match(PANELS$key, ms$analysis)])
PANELS$fold <- ifelse(PANELS$source == "nestin",
                      st$Estimate[match(PANELS$key, st$Marker)],
                      ms$fold[match(PANELS$key, ms$analysis)])
ok <- PANELS$label %in% unique(as.character(d$panel)) & !is.na(PANELS$p)
PANELS$p_holm <- NA
PANELS$p_holm[ok] <- p.adjust(PANELS$p[ok], "holm")
write.csv(PANELS[, c("label", "source", "key", "fold", "p", "p_holm")],
          file.path(OUT_DIR, "IF_markers_relative_to_NITB_pvalues.csv"), row.names = FALSE)
cat("\nP values (Holm over", sum(ok), "panels):\n"); print(PANELS[, c("label", "fold", "p", "p_holm")], row.names = FALSE)

fmt_p <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "p < 0.001", paste0("p = ", signif(p, 2))))
lab <- data.frame(panel = factor(PANELS$label, levels = PANELS$label), label = fmt_p(PANELS$p_holm))
lab <- lab[lab$panel %in% d$panel, ]

# ---- 4. Figures -----------------------------------------------------------------
# (a) one file per marker; (b) one combined figure with the panels in PANEL_ORDER on a
# GRID_NCOL x GRID_NROW grid. Every panel has its own y-axis: markers differ in range
# (e.g. ITGB4 vs tenascin), so a common axis would flatten some of them.
pt <- 1 / ggplot2::.pt
bad <- setdiff(PANEL_ORDER, PANELS$label)
if (length(bad)) stop("PANEL_ORDER has label(s) not in PANELS: ", paste(bad, collapse = ", "))

make_plot <- function(dd, nrow = NULL, ncol = NULL) {
  dd  <- droplevels(dd)
  ss  <- droplevels(sec[sec$panel %in% dd$panel, ])
  bb  <- bar[bar$panel %in% dd$panel, ]
  ll  <- lab[lab$panel %in% dd$panel, ]
  # p label just above the data of each panel (log scale)
  top <- do.call(rbind, lapply(split(dd$relative, dd$panel, drop = TRUE), function(v) {
    r <- log(range(v)); data.frame(y = exp(r[2] + 0.12 * diff(r) + 0.05)) }))
  ll$y <- top$y[match(as.character(ll$panel), rownames(top))]
  cols <- block_cols[levels(dd$block_label)]
  legend_ncol <- min(length(cols), 2 * if (is.null(ncol)) 1 else ncol)   # wrap the legend to the figure width
  ggplot(dd, aes(group, relative)) +
    geom_hline(yintercept = 1, linetype = "dashed", colour = "grey60", linewidth = 0.25) +
    geom_point(aes(colour = block_label), position = position_jitter(width = 0.12, height = 0, seed = 1),
               size = 0.6, alpha = 0.45, stroke = 0) +
    geom_line(data = ss, aes(group = paste(source, section), colour = block_label), linewidth = 0.4) +
    geom_point(data = ss, aes(fill = block_label), shape = 21, colour = "black", size = 2, stroke = 0.3) +
    geom_errorbar(data = bb, aes(ymin = relative, ymax = relative), width = 0.45, linewidth = 0.6, colour = "black") +
    geom_text(data = ll, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, size = FONT_PT * pt, vjust = 0) +
    scale_y_continuous(trans = "log2", breaks = scales::breaks_log(n = 5, base = 2),
                       labels = function(x) format(signif(x, 2), drop0trailing = TRUE, trim = TRUE),
                       expand = expansion(mult = c(0.05, 0.12))) +
    scale_colour_manual(values = cols, name = NULL) +
    scale_fill_manual(values = cols, name = NULL) +
    guides(colour = guide_legend(ncol = legend_ncol), fill = guide_legend(ncol = legend_ncol)) +
    facet_wrap(~ panel, nrow = nrow, ncol = ncol, scales = "free_y") +
    labs(x = NULL, y = "Intensity relative to NI-TB\n(same section, log2 scale)") +
    theme_classic(base_size = FONT_PT, base_family = "Helvetica") +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold", size = FONT_PT),
          axis.text = element_text(colour = "black", size = FONT_PT - 1),
          axis.line = element_line(linewidth = 0.25), axis.ticks = element_line(linewidth = 0.25),
          panel.spacing = unit(3, "mm"), legend.position = "bottom", legend.key.size = unit(3, "mm"))
}
save_fig <- function(p, file, ncol, nrow) {
  w <- FIG_PANEL_W_MM * ncol / 25.4 + 0.45                     # + y-axis title
  h <- FIG_PANEL_H_MM * nrow / 25.4 + 0.35                     # + legend
  ggsave(paste0(file, ".pdf"), p, width = w, height = h, device = pdf_device)
  ggsave(paste0(file, ".png"), p, width = w, height = h, dpi = 300, bg = "white")
}

# (a) one file per marker
for (pn in levels(droplevels(d$panel))) {
  p <- make_plot(d[d$panel == pn, ])
  save_fig(p, file.path(OUT_DIR, paste0("Fig_IF_relative_to_NITB_", gsub("[^A-Za-z0-9]+", "_", pn))), 1, 1)
}

# (b) combined figure, panels in PANEL_ORDER on the chosen grid
dc <- d[d$panel %in% PANEL_ORDER, ]
dc$panel <- factor(as.character(dc$panel), levels = intersect(PANEL_ORDER, unique(as.character(dc$panel))))
n  <- nlevels(dc$panel)
nc <- if (!is.null(GRID_NCOL)) GRID_NCOL else ceiling(n / GRID_NROW)
nr <- if (!is.null(GRID_NROW)) GRID_NROW else ceiling(n / nc)
if (nc * nr < n) stop("GRID_NCOL x GRID_NROW (", nc, " x ", nr, ") is smaller than the number of panels (", n, ")")
sec_c <- sec; bar_c <- bar; lab_c <- lab
sec$panel <- factor(as.character(sec$panel), levels = levels(dc$panel))
bar$panel <- factor(as.character(bar$panel), levels = levels(dc$panel))
lab$panel <- factor(as.character(lab$panel), levels = levels(dc$panel))
p_all <- make_plot(dc, nrow = nr, ncol = nc)
sec <- sec_c; bar <- bar_c; lab <- lab_c
grid_tag <- paste0(nc, "x", nr)
save_fig(p_all, file.path(OUT_DIR, paste0("Fig_IF_markers_relative_to_NITB_", grid_tag)), nc, nr)
cat("\nSaved: one figure per marker + combined figure (", grid_tag, ", order: ",
    paste(levels(dc$panel), collapse = " | "), ") in ", OUT_DIR, "\n", sep = "")
