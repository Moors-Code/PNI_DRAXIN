###############################################################################
# Clean HALO / Horizon feature tables for NERVE analysis
#
# For every .csv in input/Nerve_analysis and input/Nerve_analysis_with_control
# (each folder separately) this script:
#   1. Extracts Image_Name, Patient_ID, Block, ROI, ROI_number and Type
#      from "Annotation Group" (date prefixes such as "20260912_" removed;
#      combined patches such as C_12_13 get ROI_number "12;13")
#        .../Sequential IF_60456 bl5/UN_3/Cell Segmentation
#        -> Image_Name = "60456 bl5", ROI = "UN_3", Type = "Uninvaded"
#        ROI starting with "INV" (e.g. INV_1) -> Type = "Invaded"
#        Any other ROI name -> Type = "Unknown" (and a warning is printed)
#   2. Merges duplicated columns (same name, one per ROI) into one column,
#      keeping the non-empty value for each row.
#   3. Removes all "Rule-Based class group" columns (not needed for nerves).
#   4. Saves <file>_cleaned.csv into output/<folder>/cleaned/  (one per export)
#   5. COMBINES all cleaned files of the same image (e.g. batch exports
#      "part 1" + "part 2") into output/<folder>/cleaned/combined/:
#        <Image_Name>_combined.csv  and  All_images_combined.csv
#      A Source_File column records which export each cell came from.
#   6. Reports are made per IMAGE on the combined data, so ROIs split over
#      batches are not falsely reported as missing:
#        Cleaning_report.csv          one row per image
#        Cleaning_report_per_ROI.csv  one row per image x ROI (Yes/No)
#        Issues_report.csv            ROIs with any "No", missing
#                                     consecutive ROIs, and ROIs that appear
#                                     in more than one export (double export)
#      If an ROI is in several exports and one is more complete (more
#      analyses present), only that one is kept and this is reported.
#
# Base R only - no packages needed.
###############################################################################

# ---- Settings -------------------------------------------------------------
source(file.path(if (basename(getwd()) == "scripts") "." else "scripts", "00_settings.R"))

# Input folders cleaned by this script (a folder that does not exist is skipped)
nerve_folders <- c("Nerve_analysis", "Nerve_analysis_with_control")

analyses_base <- c(
  "Cells/Mean Intensity (S100b - Cy5 (1))",
  "Cells/Mean Intensity (P75NGF - TRITC (1))",
  "Cells/Mean Intensity (ITGB4 - TRITC (1))",
  "Cells/Mean Intensity (NESTIN - Cy5 (1))",
  "Cells/Dot count (Dots T3 DRAXIN - Cy5)",
  "Cells/Dot count (Dots T4 GAP43 - Cy7)"
)
# Control channels in the "_with_control" folder: any further "Cells/Mean Intensity (...)"
# or "Cells/Dot count (...)" column found there is added to the checks automatically.
# To check a fixed list instead, give the column names here.
analyses_control <- c()

# ---- Helpers --------------------------------------------------------------

# Take the first non-empty value across several columns, row by row
coalesce_cols <- function(df_part) {
  out <- rep(NA_character_, nrow(df_part))
  for (j in seq_len(ncol(df_part))) {
    v    <- df_part[[j]]
    fill <- is.na(out) & !is.na(v) & v != ""
    out[fill] <- v[fill]
  }
  out
}

# Count rows where duplicate columns hold DIFFERENT non-empty values
count_conflicts <- function(df_part) {
  m <- as.matrix(df_part)
  m[m == ""] <- NA
  sum(apply(m, 1, function(r) length(unique(r[!is.na(r)])) > 1))
}

# Parse "Annotation Group" into Image_Name / ROI / Type
parse_annotation <- function(x) {
  parts <- strsplit(x, "/", fixed = TRUE)
  get_part <- function(p, i) if (length(p) >= i) p[i] else NA_character_

  slide_folder <- vapply(parts, get_part, character(1), i = 2)
  roi          <- vapply(parts, get_part, character(1), i = 3)

  # Image name = text after the last "IF_" in the slide folder name
  # (e.g. "...HiPlex + Sequential IF_60456 bl5" -> "60456 bl5").
  # If "IF_" is missing, fall back to the text after the last "_".
  image_name <- ifelse(grepl("IF_", slide_folder),
                       sub(".*IF_", "", slide_folder),
                       sub(".*_", "", slide_folder))
  # Newer slides add a date ("IF_20260912_80830 bl4") -> drop it ("80830 bl4")
  image_name <- sub("^\\d{8}_", "", image_name)

  type <- ifelse(is.na(roi),                                   NA_character_,
          ifelse(grepl("^INV", roi, ignore.case = TRUE),    "Invaded",
          ifelse(grepl("^UN",  roi, ignore.case = TRUE),    "Uninvaded",
                                                             "Unknown")))

  image_name <- trimws(image_name)
  roi        <- trimws(roi)
  # Patient and block: "80830 bl4" -> Patient_ID "80830", Block "bl4"
  patient <- trimws(sub("\\s*bl.*$", "", image_name, ignore.case = TRUE))
  block   <- ifelse(grepl("bl\\s*\\d+", image_name, ignore.case = TRUE),
                    sub(".*?(bl\\s*\\d+).*", "\\1", image_name, ignore.case = TRUE),
                    NA_character_)
  # ROI numbers: "C_12_13" -> "12;13" (patch next to two nerves),
  #              "C_FAR_3" -> "3", "INV_5" -> "5"
  roi_num <- ifelse(grepl("(_\\d+)+$", roi),
                    gsub("_", ";", sub("^_", "", regmatches(roi, regexpr("(_\\d+)+$", roi)))),
                    NA_character_)
  # blocks whose name does not start with the patient number (00_settings.R)
  fix <- image_name %in% names(PATIENT_OF_BLOCK)
  patient[fix] <- PATIENT_OF_BLOCK[image_name[fix]]
  data.frame(Image_Name = image_name,
             Patient_ID = patient,
             Block      = block,
             ROI        = roi,
             ROI_number = roi_num,
             Type       = type,
             stringsAsFactors = FALSE)
}

# Natural sort so UN_2 comes before UN_10
# (combined patches such as C_12_13 are sorted by their first number)
nat_sort <- function(x) x[order(sub("(_\\d+)+$", "", x),
                                suppressWarnings(as.numeric(
                                  sub("^_?(\\d+).*$", "\\1",
                                      sub("^.*?((_\\d+)+)$", "\\1", x)))))]

# Stack data frames whose columns may differ (missing columns -> NA)
rbind_fill <- function(dfs) {
  all_cols <- unique(unlist(lapply(dfs, names)))
  do.call(rbind, lapply(dfs, function(d) {
    for (m in setdiff(all_cols, names(d))) d[[m]] <- NA
    d[, all_cols, drop = FALSE]
  }))
}

safe_name <- function(x) gsub("[^A-Za-z0-9._-]+", "_", x)

yes_no <- function(df, col) {
  if (col %in% names(df) && any(!is.na(df[[col]]) & df[[col]] != "")) "Yes" else "No"
}

# ---- Clean one input folder ---------------------------------------------
clean_folder <- function(input_dir, output_dir, analyses) {
  combined_dir <- file.path(output_dir, "combined")
  dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

  # ---- Main loop ------------------------------------------------------------
  files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE,
                      ignore.case = TRUE)
  files <- files[!grepl("_cleaned\\.csv$|Cleaning_report", basename(files))]

  if (length(files) == 0) stop("No CSV files found in ", input_dir)

  cleaned_list <- list()

  for (f in files) {
    message("Processing: ", basename(f))

    # Read everything as text so empty cells stay empty and names are kept as-is
    raw <- read.csv(f, check.names = FALSE, colClasses = "character",
                    na.strings = c("", "NA"), encoding = "UTF-8")
    nm  <- names(raw)
    # Byte-level clean-up so it works in any system locale:
    #   strip byte-order mark, and make units plain ASCII (μm² -> um2)
    nm  <- enc2utf8(nm)
    nm  <- gsub("\xef\xbb\xbf", "", nm, useBytes = TRUE)        # BOM
    nm  <- gsub("\xce\xbc|\xc2\xb5", "u", nm, useBytes = TRUE)  # μ / µ
    nm  <- gsub("\xc2\xb2", "2", nm, useBytes = TRUE)           # ²
    Encoding(nm) <- "UTF-8"

    # Unnamed first column = row index from the export
    nm[nm == "" | is.na(nm)] <- "Row_ID"

    # Drop Rule-Based class group columns (not needed for nerve analysis)
    keep <- !grepl("^Rule-Based class group", nm)
    raw  <- raw[, keep, drop = FALSE]
    nm   <- nm[keep]
    base_nm <- nm

    # Merge columns sharing the same (base) name, keeping first-appearance order
    merged <- list()
    conflicts <- c()
    for (u in unique(base_nm)) {
      idx  <- which(base_nm == u)
      part <- raw[, idx, drop = FALSE]
      if (length(idx) > 1) {
        n_conf <- count_conflicts(part)
        if (n_conf > 0) conflicts[u] <- n_conf
        merged[[u]] <- coalesce_cols(part)
      } else {
        merged[[u]] <- part[[1]]
      }
    }
    clean <- as.data.frame(merged, check.names = FALSE, stringsAsFactors = FALSE)

    if (length(conflicts) > 0) {
      warning(basename(f), ": rows with conflicting values in merged columns -> ",
              paste0(names(conflicts), " (", conflicts, ")", collapse = "; "),
              ". The first non-empty value was kept.")
    }

    # Convert numeric / logical columns back to proper types
    clean[] <- lapply(clean, function(v) type.convert(v, as.is = TRUE))

    # Add Image_Name, ROI, Type right after Annotation Group
    ann   <- parse_annotation(clean[["Annotation Group"]])
    pos   <- match("Annotation Group", names(clean))
    clean <- cbind(clean[, seq_len(pos), drop = FALSE], ann,
                   clean[, -seq_len(pos), drop = FALSE])

    if (any(ann$Type == "Unknown", na.rm = TRUE)) {
      warning(basename(f), ": ROI(s) not starting with INV or UN -> Type 'Unknown': ",
              paste(unique(ann$ROI[ann$Type == "Unknown"]), collapse = ", "))
    }

    # Save cleaned file (one per export)
    out_name <- sub("\\.csv$", "_cleaned.csv", basename(f), ignore.case = TRUE)
    write.csv(clean, file.path(output_dir, out_name),
              row.names = FALSE, na = "", fileEncoding = "UTF-8")

    # Keep it for combining; remember which export each cell came from
    clean <- cbind(Source_File = basename(f), clean, stringsAsFactors = FALSE)
    cleaned_list[[f]] <- clean
  }

  # ---- Combine cleaned files per image --------------------------------------
  all_data <- rbind_fill(cleaned_list)
  rownames(all_data) <- NULL
  images   <- sort(unique(all_data$Image_Name))

  # ---- ROIs exported more than once (e.g. a missing ROI re-exported) --------
  # For every image x ROI found in more than one export file:
  #   * if one export has MORE analyses present than all others, keep only that
  #     one (e.g. today's empty C_5 is replaced by tomorrow's re-export of C_5);
  #   * otherwise keep all copies and flag them as a possible double export.
  replaced_log <- list()
  n_present <- function(df) sum(vapply(analyses, function(a) yes_no(df, a) == "Yes",
                                       logical(1)))
  drop_rows <- rep(FALSE, nrow(all_data))
  key       <- paste(all_data$Image_Name, all_data$ROI, sep = "\r")
  for (k in unique(key)) {
    idx  <- which(key == k)
    srcs <- unique(all_data$Source_File[idx])
    if (length(srcs) < 2) next
    score <- vapply(srcs, function(s0)
      n_present(all_data[idx[all_data$Source_File[idx] == s0], , drop = FALSE]),
      numeric(1))
    best <- srcs[score == max(score)]
    if (length(best) == 1) {
      worse <- setdiff(srcs, best)
      drop_rows[idx[all_data$Source_File[idx] %in% worse]] <- TRUE
      replaced_log[[length(replaced_log) + 1]] <- data.frame(
        Image_Name = all_data$Image_Name[idx[1]], ROI = all_data$ROI[idx[1]],
        Type = all_data$Type[idx[1]],
        Issue = "ROI replaced by more complete export",
        N_cells = sum(all_data$Source_File[idx] == best),
        Details = paste0("Kept: ", best, " (", max(score), "/", length(analyses),
                         " analyses). Dropped: ", paste(worse, collapse = " | "),
                         " (", paste(score[worse], collapse = ", "), "/",
                         length(analyses), ")"),
        stringsAsFactors = FALSE)
    }
  }
  if (any(drop_rows)) {
    message(sum(drop_rows), " cell(s) dropped from less complete duplicate exports")
    all_data <- all_data[!drop_rows, , drop = FALSE]
    rownames(all_data) <- NULL
  }

  for (img in images) {
    d <- all_data[all_data$Image_Name == img, , drop = FALSE]
    write.csv(d, file.path(combined_dir, paste0(safe_name(img), "_combined.csv")),
              row.names = FALSE, na = "", fileEncoding = "UTF-8")
    n_src <- length(unique(d$Source_File))
    message("Image '", img, "': ", n_src, " export file(s) combined, ",
            nrow(d), " cells")
  }
  write.csv(all_data, file.path(combined_dir, "All_images_combined.csv"),
            row.names = FALSE, na = "", fileEncoding = "UTF-8")

  # ---- Reports per image (on combined data) ---------------------------------
  report_file <- list()
  report_roi  <- list()

  for (img in images) {
    d        <- all_data[all_data$Image_Name == img, , drop = FALSE]
    roi_type <- unique(d[, c("ROI", "Type")])
    rois     <- nat_sort(unique(d$ROI))
    inv      <- nat_sort(roi_type$ROI[roi_type$Type == "Invaded"])
    uninv    <- nat_sort(roi_type$ROI[roi_type$Type == "Uninvaded"])
    unk      <- nat_sort(roi_type$ROI[roi_type$Type == "Unknown"])
    srcs     <- unique(d$Source_File)

    rep_row <- data.frame(
      Image_Name        = img,
      N_source_files    = length(srcs),
      Source_files      = paste(srcs, collapse = " | "),
      N_cells           = nrow(d),
      N_unique_ROIs     = length(rois),
      N_Invaded_ROIs    = length(inv),
      Invaded_ROIs      = paste(inv, collapse = ", "),
      N_Uninvaded_ROIs  = length(uninv),
      Uninvaded_ROIs    = paste(uninv, collapse = ", "),
      N_Unknown_ROIs    = length(unk),
      stringsAsFactors  = FALSE, check.names = FALSE
    )
    for (a in analyses) rep_row[[a]] <- yes_no(d, a)
    report_file[[img]] <- rep_row

    for (r in rois) {
      sub_df <- d[d$ROI == r, , drop = FALSE]
      rr <- data.frame(Image_Name = img, ROI = r, Type = sub_df$Type[1],
                       N_cells = nrow(sub_df),
                       Source_files = paste(unique(sub_df$Source_File),
                                            collapse = " | "),
                       N_source_files = length(unique(sub_df$Source_File)),
                       stringsAsFactors = FALSE, check.names = FALSE)
      for (a in analyses) rr[[a]] <- yes_no(sub_df, a)
      report_roi[[paste(img, r)]] <- rr
    }
  }

  report_file <- do.call(rbind, report_file)
  report_roi  <- do.call(rbind, report_roi)
  rownames(report_file) <- NULL
  rownames(report_roi)  <- NULL

  write.csv(report_file, file.path(output_dir, "Cleaning_report.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(report_roi,  file.path(output_dir, "Cleaning_report_per_ROI.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")

  # ---- Final issues report (per image, on combined data) --------------------
  #   * ROIs with "No" for at least one analysis
  #   * missing consecutive ROI numbers in each series (INV, UN, ...),
  #     checked from 1 up to the highest number present
  #   * ROIs found in more than one export file (possible double export)
  issues <- replaced_log
  add_issue <- function(img, roi, type, issue, n, details) {
    issues[[length(issues) + 1]] <<- data.frame(
      Image_Name = img, ROI = roi, Type = type, Issue = issue,
      N_cells = n, Details = details, stringsAsFactors = FALSE)
  }

  for (img in images) {
    rr <- report_roi[report_roi$Image_Name == img, , drop = FALSE]

    for (i in seq_len(nrow(rr))) {
      missing_a <- analyses[unlist(rr[i, analyses]) == "No"]
      if (length(missing_a) > 0)
        add_issue(img, rr$ROI[i], rr$Type[i],
                  if (length(missing_a) == length(analyses))
                    "All analyses missing" else "Analysis missing",
                  rr$N_cells[i], paste(missing_a, collapse = "; "))
      if (rr$N_source_files[i] > 1)
        add_issue(img, rr$ROI[i], rr$Type[i], "ROI in more than one file",
                  rr$N_cells[i],
                  paste("Possible double export - check:", rr$Source_files[i]))
    }

    # Combined patches (e.g. C_8_9, C_12_13) count for every number they contain
    prefix0 <- sub("(_\\d+)+$", "", rr$ROI)
    nums    <- lapply(regmatches(rr$ROI, regexpr("(_\\d+)+$", rr$ROI)),
                      function(x) as.integer(strsplit(sub("^_", "", x), "_")[[1]]))
    has_num <- grepl("(_\\d+)+$", rr$ROI)
    prefix  <- rep(prefix0[has_num], lengths(nums))
    num     <- unlist(nums)
    for (p in unique(prefix[!is.na(num)])) {
      present <- sort(num[prefix == p & !is.na(num)])
      gaps    <- setdiff(seq_len(max(present)), present)
      for (g in gaps)
        add_issue(img, paste0(p, "_", g), rr$Type[prefix0 == p][1],
                  "ROI missing", 0,
                  paste0("Not in any export of this image (series ", p,
                         "_1 to ", p, "_", max(present), ")"))
    }
  }

  if (length(issues) > 0) {
    issues <- do.call(rbind, issues)
    issues <- issues[order(issues$Image_Name, issues$Issue), ]
    rownames(issues) <- NULL
    message("\n", nrow(issues), " issue(s) found - see Issues_report.csv:")
    print(issues[, c("Image_Name", "ROI", "Type", "Issue")], row.names = FALSE)
  } else {
    issues <- data.frame(Image_Name = character(), ROI = character(),
                         Type = character(), Issue = character(),
                         N_cells = integer(), Details = character())
    message("\nNo issues found: all ROIs complete and consecutive.")
  }

  # Images without any issue are listed too, so every image appears once
  ok_img <- setdiff(images, issues$Image_Name)
  if (length(ok_img) > 0)
    issues <- rbind(issues, data.frame(Image_Name = ok_img, ROI = "", Type = "",
                                       Issue = "No issues", N_cells = NA,
                                       Details = "", stringsAsFactors = FALSE))

  write.csv(issues, file.path(output_dir, "Issues_report.csv"),
            row.names = FALSE, na = "", fileEncoding = "UTF-8")

  message("Done. ", length(files), " file(s) cleaned, ", length(images),
          " image(s) combined. Output in: ", output_dir)
  invisible(NULL)
}

# ---- Run for each input folder -------------------------------------------------
for (folder in nerve_folders) {
  input_dir <- file.path("input", folder)
  if (!dir.exists(input_dir)) { message("Skipped (folder not found): ", input_dir); next }
  if (!length(list.files(input_dir, "\\.csv$", ignore.case = TRUE))) {
    message("Skipped (no CSV files): ", input_dir); next }
  analyses <- analyses_base
  if (grepl("control", folder)) {
    extra <- analyses_control
    if (!length(extra)) {                       # detect control columns from the headers
      hdr <- unique(unlist(lapply(list.files(input_dir, "\\.csv$", full.names = TRUE, ignore.case = TRUE),
                                  function(f) names(read.csv(f, nrows = 1, check.names = FALSE)))))
      extra <- setdiff(grep("^Cells/(Mean Intensity|Dot count) \\(", hdr, value = TRUE), analyses_base)
    }
    if (length(extra)) message("Control columns checked as well: ", paste(extra, collapse = "; "))
    analyses <- c(analyses, extra)
  }
  message("\n==== ", folder, " ====")
  clean_folder(input_dir, file.path("output", folder, "cleaned"), analyses)
}

