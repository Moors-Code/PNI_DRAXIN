# ============================================================
# Batch extract ROI + stable markers + RB columns (detected by CONTENT)
# Column NAMES may vary; RB VALUES are stable across files.
# ============================================================

library(readr)
library(dplyr)
library(stringr)
library(purrr)
library(tidyr)

# ---- set paths ----
# Folder holding the raw per-ROI exports of the COMET immune panel (one CSV per
# ROI, with the ROI number in the file name). Set this before running.
RAW_EXPORT_DIR <- "<path to the raw COMET immune-panel ROI exports>"

in_dir   <- RAW_EXPORT_DIR
out_file <- file.path("input", "combined_selected_columns.csv")

if (!dir.exists(in_dir))
  stop("Set RAW_EXPORT_DIR to the folder with the raw ROI exports (currently: ", in_dir, ")")

# ---- ROI from filename (ROI1, ROI12, ...) ----
extract_roi <- function(path) {
  fn <- basename(path)
  m <- str_match(fn, "(?i)ROI\\s*\\d+")[, 1]
  if (is.na(m)) return(NA_character_)
  str_replace_all(toupper(m), "\\s+", "")
}

# ---- RB detection by CONTENT (values) ----
detect_rb_columns <- function(dat) {
  
  # any RB columns (names can be "Rule.Based.class.group", "Rule-Based class group", with dots, etc.)
  nms <- names(dat)
  rb_cols <- nms[str_detect(tolower(nms), "rule") &
                   str_detect(tolower(nms), "based") &
                   str_detect(tolower(nms), "class") &
                   str_detect(tolower(nms), "group")]
  
  RB_CD45   <- NA_character_
  RB_Immune <- NA_character_
  RB_PDL1   <- NA_character_
  
  # helper to get unique non-NA strings
  uniq_vals <- function(x) unique(trimws(as.character(x[!is.na(x) & x != ""])))
  
  for (cn in rb_cols) {
    u <- uniq_vals(dat[[cn]])
    
    # --- CD45 column: contains BOTH CD45+ and CD45- ---
    if (is.na(RB_CD45) && ("CD45+" %in% u) && ("CD45-" %in% u)) {
      RB_CD45 <- cn
      next
    }
    
    # --- PDL1 combo column: contains CD45+PDL1- / CD45+PDL1+ / CD45-PDL1- ---
    if (is.na(RB_PDL1) &&
        any(c("CD45+PDL1-", "CD45+PDL1+", "CD45-PDL1-") %in% u)) {
      RB_PDL1 <- cn
      next
    }
    
    # --- Immune types column: contains any of these ---
    immune_markers <- c(
      "CD4 T cells", "CD8 T cells", "CD8- CD4- T cells", "CD20+ cells",
      "Other immune cells", "CD11c DCs", "CD68+ macro", "NK cells", "non immune cells"
    )
    if (is.na(RB_Immune) && any(immune_markers %in% u)) {
      RB_Immune <- cn
      next
    }
  }
  
  list(
    rb_cols = rb_cols,
    RB_CD45 = RB_CD45,
    RB_Immune = RB_Immune,
    RB_PDL1 = RB_PDL1
  )
}

# ---- stable marker columns by flexible pattern matching (your files use dots) ----
# This matches the "Nuclei.Mean.Intensity..T4.IFNG." style etc.
pick_first <- function(nms, pattern) {
  hit <- grep(pattern, nms, value = TRUE, ignore.case = TRUE)
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

get_stable_cols <- function(dat) {
  nms <- names(dat)
  
  c(
    Nuc_TGFB1 = pick_first(nms, "Nuclei.*Mean.*Intensity.*T2.*TGFB1"),
    Nuc_GZMB  = pick_first(nms, "Nuclei.*Mean.*Intensity.*T3.*GZMB"),
    Nuc_IFNG  = pick_first(nms, "Nuclei.*Mean.*Intensity.*T4.*IFNG"),
    Nuc_IL12A = pick_first(nms, "Nuclei.*Mean.*Intensity.*T6.*IL12A"),
    Nuc_HIF1A = pick_first(nms, "Nuclei.*Mean.*Intensity.*T7.*HIF1A"),
    Nuc_IL1B  = pick_first(nms, "Nuclei.*Mean.*Intensity.*T11.*IL1B"),
    Nuc_CXCL9 = pick_first(nms, "Nuclei.*Mean.*Intensity.*T12.*CXCL9"),
    
    Nuc_KI67  = pick_first(nms, "Nuclei.*Mean.*Intensity.*KI67"),
    
    Nuc_IFNG_dots  = pick_first(nms, "Nuclei.*Dot.*count.*Dots.*T4.*IFNG"),
    Nuc_CXCL9_dots = pick_first(nms, "Nuclei.*Dot.*count.*Dots.*T12.*CXCL9"),
    Nuc_IL1B_dots  = pick_first(nms, "Nuclei.*Dot.*count.*Dots.*T11.*IL1B")
  )
}

# ---- list files ----
files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No CSV files found in: ", in_dir)

# ---- read + combine ----
combined <- map_dfr(files, function(f) {
  
  dat <- suppressMessages(read_csv(f, show_col_types = FALSE))
  roi <- extract_roi(f)
  
  rb <- detect_rb_columns(dat)
  st <- get_stable_cols(dat)
  
  message("\nFILE: ", basename(f))
  message("  ROI_ID: ", roi)
  message("  RB_CD45:   ", rb$RB_CD45)
  message("  RB_Immune: ", rb$RB_Immune)
  message("  RB_PDL1:   ", rb$RB_PDL1)
  
  # warn if missing
  if (is.na(rb$RB_CD45))   warning("Missing RB_CD45 in: ", basename(f))
  if (is.na(rb$RB_Immune)) warning("Missing RB_Immune in: ", basename(f))
  if (is.na(rb$RB_PDL1))   warning("Missing RB_PDL1 in: ", basename(f))
  
  missing_st <- names(st)[is.na(st)]
  if (length(missing_st) > 0) {
    warning("Missing stable cols in ", basename(f), ": ", paste(missing_st, collapse = " | "))
  }
  
  dat %>%
    transmute(
      ROI_ID = roi,
      
      RB_CD45   = if (!is.na(rb$RB_CD45))   trimws(as.character(.data[[rb$RB_CD45]]))   else NA_character_,
      RB_Immune = if (!is.na(rb$RB_Immune)) trimws(as.character(.data[[rb$RB_Immune]])) else NA_character_,
      RB_PDL1   = if (!is.na(rb$RB_PDL1))   trimws(as.character(.data[[rb$RB_PDL1]]))   else NA_character_,
      
      Nuc_TGFB1 = if (!is.na(st[["Nuc_TGFB1"]])) .data[[st[["Nuc_TGFB1"]]]] else NA_real_,
      Nuc_GZMB  = if (!is.na(st[["Nuc_GZMB"]]))  .data[[st[["Nuc_GZMB"]]]]  else NA_real_,
      Nuc_IFNG  = if (!is.na(st[["Nuc_IFNG"]]))  .data[[st[["Nuc_IFNG"]]]]  else NA_real_,
      Nuc_IL12A = if (!is.na(st[["Nuc_IL12A"]])) .data[[st[["Nuc_IL12A"]]]] else NA_real_,
      Nuc_HIF1A = if (!is.na(st[["Nuc_HIF1A"]])) .data[[st[["Nuc_HIF1A"]]]] else NA_real_,
      Nuc_IL1B  = if (!is.na(st[["Nuc_IL1B"]]))  .data[[st[["Nuc_IL1B"]]]]  else NA_real_,
      Nuc_CXCL9 = if (!is.na(st[["Nuc_CXCL9"]])) .data[[st[["Nuc_CXCL9"]]]] else NA_real_,
      
      Nuc_KI67  = if (!is.na(st[["Nuc_KI67"]]))  .data[[st[["Nuc_KI67"]]]]  else NA_real_,
      
      Nuc_IFNG_dots  = if (!is.na(st[["Nuc_IFNG_dots"]]))  .data[[st[["Nuc_IFNG_dots"]]]]  else NA_real_,
      Nuc_CXCL9_dots = if (!is.na(st[["Nuc_CXCL9_dots"]])) .data[[st[["Nuc_CXCL9_dots"]]]] else NA_real_,
      Nuc_IL1B_dots  = if (!is.na(st[["Nuc_IL1B_dots"]]))  .data[[st[["Nuc_IL1B_dots"]]]]  else NA_real_
    )
})

# ---- write ----
dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
write_csv(combined, out_file)

message("\nDone. Wrote: ", out_file)
message("Rows combined: ", nrow(combined))
message("Unique ROI_IDs: ", paste(sort(unique(na.omit(combined$ROI_ID))), collapse = ", "))
