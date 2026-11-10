# ======================================================================
# 02_spatial_association_with_nerve_step1_lee_sham.R
#
# Spatial association of candidate genes with INVADED nerves:
# Lee's L, computed per patient, and calibrated against relocated sham
# masks of the same nerve shapes in the same section.
#
# QUESTION
#   Is a gene's expression more strongly spatially associated with the
#   invaded-nerve mask than with masks of exactly the same shapes placed
#   elsewhere in the same diseased section?
#
# DESIGN
#   * One Diseased section per patient; everything is computed per section
#     (= per patient). No statistic or permutation mixes sections.
#   * Mask: all invaded-nerve spots of the section.
#   * Population (spots entering the statistic): the section's spots minus
#       - non-invaded nerve spots,
#       - the interface of every real nerve (contact rule of the tumour-
#         proximal definition: nneighbors = 2, i.e. one of the 12 nearest
#         spots within 46 coordinate units is a nerve spot; 98-205 um),
#       - for a sham mask, its own interface (same contact rule), so real
#         and sham masks are both surrounded by an excluded ring.
#     Every spot contacting any nerve is excluded (no assignment of a spot
#     to a single nerve).
#   * Graph: k = 6 nearest neighbours within the section, restricted to
#     the population and row-standardised (as spdep::subset.listw).
#   * Statistic: Lee's L between expression (logcounts) and the 0/1 mask
#     (formula of spdep::lee; checked against spdep::lee below).
#   * Null: N_SHAM full sham masks per section. In each, every invaded
#     nerve's shape is relocated (rules as for the plasma shams: centroid
#     >= 300 um from the original, snap tolerance 0.6 x spot spacing,
#     >= 90 % of spots snapped; the shape is moved by a spot-to-spot shift so
#     it keeps its exact form on the spot grid) onto eligible spots only (no nerve, no
#     interface, no overlap with other shapes of the same mask). Replicates
#     with a failed placement are discarded and redrawn; all attempts are
#     recorded.
#   * Per gene and patient, the real L is RANKED among its N_SHAM sham L
#     values: u = (#{sham < real} + 0.5) / (N_SHAM + 1). If real and sham
#     masks are exchangeable (no association with the nerve location), the
#     rank is uniform. Descriptive per patient: z = (L_real - mean sham) /
#     SD sham, delta_L = L_real - mean sham, and whether L_real exceeds the
#     sham 95th percentile.
#   * PRIMARY test (prespecified, nerve mask, DE_pass family): rank scores
#     qnorm(u) are combined over patients with equal weight,
#     S = sum qnorm(u) / sqrt(n_patients), and p = P(S_null >= S), where the
#     null distribution of S is enumerated EXACTLY from uniform ranks in each
#     patient (all (N_SHAM + 1)^n_patients combinations). The test is exact
#     under exchangeability, needs no normality or SD estimate, and a single
#     extreme patient cannot dominate (each rank score is bounded), so strong
#     evidence requires consistency across patients. Smallest attainable p =
#     (1 / (N_SHAM + 1))^n_patients. BH across genes within each family;
#     selected = padj < ALPHA.
#   * SECONDARY (description, not selection): the same test for the
#     INTERFACE -- mask = interface ring of the invaded nerves, population =
#     background + that ring (nerve spots removed); null = the rings around
#     the same sham masks (restricted to background spots). Genes are then
#     described as nerve only / nerve + interface / interface only / neither.
#   * NULL AUDIT (code check): N_AUDIT held-out sham triplets are put through
#     the identical pipeline as if they were the real nerves (ranked among
#     the remaining shams). ~5 % of genes should have p < 0.05 and BH should
#     find ~0 genes. This verifies the procedure; exchangeability of real and
#     sham masks is supported by the matched mask QC, not by the audit.
#   * Real and sham masks are built the same way: the mask is an island
#     (its own contact ring removed) inside the same background; population
#     size, mask size, ring size, graph edges and isolated mask spots are
#     recorded for every real and sham mask (lee_sham_mask_QC.csv).
#   * CORE SIGNATURE (prespecified before the 200-sham results): DE_pass
#     genes that pass the combined test AND whose L_real exceeds the 95th
#     percentile of their own shams in EVERY patient (signature_core).
#     All DE_pass genes with their scores, core flag included, are in
#     lee_sham_DE_genes_scores.csv.
#   * No spot-level permutation (lee.mc): with ~10^4 spatially
#     autocorrelated spots it rejects for almost any patterned gene.
#   * Gene families (BH separately): DE_pass (FDR < 0.1, log2FC > 0 in InvN
#     vs NI, TableS3) = PRIMARY family, the only source of the signature;
#     LR_only (candidate pool genes via the LR analysis) = EXPLORATORY.
#
# INPUTS (read only): input/Integrated_NMF.rds, input/DE_SPE_annotaion.csv,
#   input/TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv,
#   input/gene_pool_collapsed_with_provenance.csv,
#   input/um_per_pixel_CALIBRATED.rds
# OUTPUTS: output/02_spatial_association/tables/lee_sham_all_candidates.csv,
#   lee_sham_DE_genes_scores.csv (read by steps 2-3), lee_sham_LR_only_scores.csv,
#   lee_sham_signature_core.csv (read by steps 4-5), lee_sham_null_audit.csv,
#   lee_sham_placement_QC.csv,
#   lee_sham_mask_QC.csv, processed/lee_sham_per_patient.rds
# ======================================================================

library(SpatialExperiment)
library(Matrix)
library(RANN)
library(dplyr)
library(tidyr)

select <- dplyr::select; filter <- dplyr::filter; rename <- dplyr::rename
mutate <- dplyr::mutate; summarise <- dplyr::summarise; distinct <- dplyr::distinct
group_by <- dplyr::group_by; arrange <- dplyr::arrange

# ---- 0. Directory structure -------------------------------------------
# Paths are relative to the project folder (the one containing input/).
# Run from there, or from scripts/ -- this moves up one level if needed.
if (basename(getwd()) == "scripts") setwd("..")
if (!dir.exists("input"))
  stop("Run this script from the project folder (the one containing input/).")

PROJECT_DIR   <- "."
INPUT_DIR     <- file.path(PROJECT_DIR, "input")
PROCESSED_DIR <- file.path(PROJECT_DIR, "processed")
TABLES_DIR    <- file.path(PROJECT_DIR, "output", "02_spatial_association", "tables")
dir.create(TABLES_DIR,    recursive = TRUE, showWarnings = FALSE)
dir.create(PROCESSED_DIR, recursive = TRUE, showWarnings = FALSE)

N_SHAM            <- 200    # full sham masks per section (set e.g. 20 for a quick test run)
N_AUDIT           <- 50     # null audit: held-out sham triplets analysed as if they were real
MAX_TRIES_SHAPE   <- 200    # placement attempts per nerve shape
MAX_REDRAWS       <- 5 * N_SHAM
MIN_RELOC_UM      <- 300
SNAP_TOL_FRAC     <- 0.6
MIN_SNAP_FRAC     <- 0.9
K_GRAPH           <- 6
IF_NNEIGHBORS     <- 2      # contact rule of the tumour-proximal definition
IF_RADIUS         <- IF_NNEIGHBORS * 23
IF_K              <- IF_NNEIGHBORS * 6 + 1
ALPHA             <- 0.05
INCLUDE_LR_ONLY   <- TRUE
EXPR_ASSAY        <- "logcounts"
SEED              <- 42
PATIENT_LABELS    <- c("54703" = "P1", "6249" = "P2", "70865" = "P3")

# ---- 1. Data -------------------------------------------------------------
# Reading the object is the slow step: an `spe` already in memory is reused
# (clear everything else first with rm(list = setdiff(ls(), "spe"))). The spot
# annotation is always re-read from the csv, so it is never stale.
if (!exists("spe") || !inherits(spe, "SpatialExperiment")) {
  spe <- readRDS(file.path(INPUT_DIR, "Integrated_NMF.rds"))
} else cat("Using the spe already in memory (", ncol(spe), " spots)\n", sep = "")
# Attach the spot annotation BY BARCODE, checked against the object's own
# identifiers before anything is overwritten (spatial coordinates belong to
# the object's columns, so a row-order mismatch would silently move labels).
attach_annotation <- function(spe, path) {
  ann <- read.csv(path)
  if (nrow(ann) != ncol(spe)) stop("Annotation has ", nrow(ann), " rows but the object has ", ncol(spe), " spots.")
  if (anyDuplicated(ann$barcode)) stop("Duplicated barcodes in the annotation file.")
  keys <- list(colData_barcode = if ("barcode" %in% names(colData(spe))) as.character(spe$barcode),
               colnames = colnames(spe))
  keys <- Filter(function(k) !is.null(k) && !anyNA(match(k, ann$barcode)), keys)
  if (!length(keys)) stop("Annotation barcodes do not match the object's barcodes or column names.")
  ann <- ann[match(keys[[1]], ann$barcode), , drop = FALSE]
  if ("sample_id" %in% names(colData(spe)) && !all(as.character(spe$sample_id) == as.character(ann$sample_id)))
    stop("Annotation sample_id does not match the object's sample_id for the same barcodes.")
  colData(spe) <- DataFrame(ann); colnames(spe) <- ann$barcode
  spe
}
spe <- attach_annotation(spe, file.path(INPUT_DIR, "DE_SPE_annotaion.csv"))
stopifnot(all(c("NerveID", "Invaded", "context", "sample_id", "patient_id") %in% colnames(colData(spe))),
          EXPR_ASSAY %in% assayNames(spe))
um_per_pixel <- readRDS(file.path(INPUT_DIR, "um_per_pixel_CALIBRATED.rds"))
coords   <- spatialCoords(spe)
nerve_id <- as.character(spe$NerveID)

# ---- 2. Gene families (straight from the DE table) ----------------------
full_de <- read.csv(file.path(INPUT_DIR, "TableS3_Full_DE_invaded_vs_non_invaded_nerve.csv"), check.names = FALSE)
names(full_de) <- trimws(names(full_de))
full_de <- full_de %>%
  transmute(Gene = trimws(Gene), DE_log2FC = `log2FC (logarithm of fold change)`, DE_padj = `adj.P.Val (FDR)`,
            DE_pass = !is.na(DE_padj) & !is.na(DE_log2FC) & DE_padj < 0.1 & DE_log2FC > 0) %>%
  distinct(Gene, .keep_all = TRUE)
gene_pool <- read.csv(file.path(INPUT_DIR, "gene_pool_collapsed_with_provenance.csv"))
# Symbol correction: the LR resource lists placental growth factor's pair as
# "PIGF^FLT1", but the official symbol of placental growth factor is PGF
# (NCBI 5228); PIGF (NCBI 5281) is an unrelated GPI-anchor gene. The pool
# entry is therefore tested as PGF.
SYMBOL_FIXES <- c(PIGF = "PGF")
fix <- trimws(gene_pool$Gene) %in% names(SYMBOL_FIXES)
if (any(fix)) {
  cat("Gene pool symbol corrections:", paste(names(SYMBOL_FIXES), SYMBOL_FIXES, sep = " -> ", collapse = ", "), "\n")
  gene_pool$Gene[fix] <- SYMBOL_FIXES[trimws(gene_pool$Gene[fix])]
}
de_pass_genes <- full_de$Gene[full_de$DE_pass]
lr_only_genes <- setdiff(trimws(gene_pool$Gene), de_pass_genes)
candidates <- data.frame(Gene = c(de_pass_genes, if (INCLUDE_LR_ONLY) lr_only_genes)) %>%
  left_join(full_de, by = "Gene") %>%
  mutate(family = ifelse(Gene %in% de_pass_genes, "DE_pass", "LR_only"))
genes <- intersect(candidates$Gene, rownames(spe))
cat("Candidates:", nrow(candidates), "| in expression matrix:", length(genes), "\n")
if (length(genes) < 0.9 * nrow(candidates)) stop("Gene identifiers do not match rownames(spe).")

# ---- 3. Helpers ------------------------------------------------------------
# Contacts (border-script rule) within one section; returns list of global indices
contacts_in_section <- function(idx) {
  nn <- nn2(coords[idx, , drop = FALSE], coords[idx, , drop = FALSE],
            k = min(IF_K, length(idx)), searchtype = "radius", radius = IF_RADIUS)
  lapply(seq_along(idx), function(i) { nb <- nn$nn.idx[i, -1]; idx[nb[nb != 0]] })
}
# spots (global idx) whose contacts include any spot of `mask` (exact rule)
ring_of <- function(mask, idx, neigh) {
  inm <- rep(FALSE, ncol(spe)); inm[mask] <- TRUE
  idx[vapply(idx, function(i) !inm[i] && any(inm[neigh[[i]]]), logical(1))]
}
# Row-standardised kNN graph of the section, restricted to `pop` (global idx)
graph_on <- function(knn_idx, pop) {
  keep <- rep(FALSE, ncol(spe)); keep[pop] <- TRUE
  pos  <- integer(ncol(spe)); pos[pop] <- seq_along(pop)
  i <- rep(seq_along(pop), each = K_GRAPH); j <- as.vector(t(knn_idx[pop, , drop = FALSE]))
  ok <- keep[j]
  A <- sparseMatrix(i = i[ok], j = pos[j[ok]], x = 1, dims = c(length(pop), length(pop)))
  rs <- rowSums(A)
  Diagonal(x = ifelse(rs > 0, 1 / rs, 0)) %*% A
}
# Lee's L of every gene with the 0/1 mask, over population `pop` (global idx).
# Same formula as spdep::lee, L = n/S2 * sum_i lag(x)_i lag(y)_i / (|x - xbar| |y - ybar|),
# rearranged so the section's expression matrix Xs (spots x genes, dense,
# built once per section) is never copied or centred per mask:
#   sum_i lag(x)_i lag(y)_i = (W' lag(y))' x - xbar * sum_i lag(y)_i rowsum_i
lee_L <- function(W, pop, mask, Xs, X2s, sec_pos) {
  pl <- sec_pos[pop]; n <- length(pop)
  ind <- numeric(nrow(Xs)); ind[pl] <- 1
  mu  <- as.vector(crossprod(Xs, ind)) / n
  sq  <- as.vector(crossprod(X2s, ind))
  ssx <- sq - n * mu^2
  y  <- as.numeric(pop %in% mask); yc <- y - mean(y)
  rs <- rowSums(W); S2 <- sum(rs^2)
  ly <- as.vector(W %*% yc)
  v  <- numeric(nrow(Xs)); v[pl] <- as.vector(crossprod(W, ly))
  num <- as.vector(crossprod(Xs, v)) - mu * sum(ly * rs)
  L <- (n / S2) * num / (sqrt(pmax(ssx, 0)) * sqrt(sum(yc^2)))
  L[!(ssx > 1e-9 * pmax(sq, 1))] <- NA                                  # no variance in this population
  setNames(L, colnames(Xs))
}

# Graph/population description of one comparison (for the mask QC table)
mask_qc <- function(W, pop, mask, ring_n, s, mask_name, test) {
  pos <- match(mask, pop); rs <- rowSums(W)
  data.frame(sample_id = s, mask = mask_name, test = test, population = length(pop),
             mask_spots = length(mask), ring_spots = ring_n, edges = sum(W != 0),
             isolated_mask_spots = sum(rs[pos] == 0), isolated_all = sum(rs == 0))
}
# Per gene and patient: real L vs its shams (descriptive summary)
calibrate <- function(L_real, L_sham, test) {
  mu <- colMeans(L_sham, na.rm = TRUE); sdv <- apply(L_sham, 2, sd, na.rm = TRUE)
  q95 <- apply(L_sham, 2, quantile, 0.95, na.rm = TRUE)
  z <- (L_real - mu) / sdv; z[!is.finite(z)] <- NA                    # zero sham variance -> NA, logged later
  M <- nrow(L_sham); bad <- is.na(L_real) | colSums(is.na(L_sham)) > 0
  below <- colSums(sweep(L_sham, 2, L_real, "<"), na.rm = TRUE)
  ties  <- colSums(sweep(L_sham, 2, L_real, "=="), na.rm = TRUE)
  u <- (below + 0.5 * ties + 0.5) / (M + 1); u[bad] <- NA
  data.frame(Gene = colnames(L_sham), test = test,
             L_real = L_real, L_sham_mean = mu, L_sham_sd = sdv, L_sham_q95 = q95,
             delta_L = L_real - mu, z = z, rank_u = u, rank_score = qnorm(u), n_sham = M,
             emp_p = (1 + colSums(sweep(L_sham, 2, L_real, ">="), na.rm = TRUE)) / (nrow(L_sham) + 1),
             exceeds_q95 = L_real > q95, row.names = NULL)
}
# Exact null distribution of S = sum qnorm(u_p) / sqrt(P) under uniform ranks;
# returns a function giving P(S_null >= S)
exact_null <- function(Ms) {
  sc <- lapply(Ms, function(M) qnorm(((0:M) + 0.5) / (M + 1)))
  tot <- Reduce(function(a, b) as.vector(outer(a, b, "+")), sc) / sqrt(length(Ms))
  tot <- sort(tot); N <- length(tot)
  function(S) ifelse(is.na(S), NA, (N - findInterval(S - 1e-9, tot)) / N)
}

# ---- 4. Per section: masks, shams, Lee's L ------------------------------
dis_samples <- unique(spe$sample_id[spe$context == "Diseased"])
pts_of <- sapply(dis_samples, function(s) unique(as.character(spe$patient_id[spe$sample_id == s])))
if (any(duplicated(pts_of))) stop("More than one Diseased section per patient -- per-patient design assumes one.")
expr_all <- assay(spe, EXPR_ASSAY)[genes, , drop = FALSE]

set.seed(SEED)
per_patient <- list(); placement_qc <- list(); lee_check <- NULL; sham_store <- list(); mask_qc_all <- list()
for (s in dis_samples) {
  pid <- pts_of[[s]]
  idx <- which(spe$sample_id == s)
  if (!all(spe$context[idx] %in% "Diseased")) stop(s, ": section contains spots not annotated as Diseased.")
  inv_nerves <- unique(nerve_id[idx][!is.na(nerve_id[idx]) & spe$Invaded[idx] %in% "Invaded"])
  n_na_inv <- sum(!is.na(nerve_id[idx]) & is.na(spe$Invaded[idx]))
  if (n_na_inv > 0) warning(s, ": ", n_na_inv, " nerve spots with missing 'Invaded' label -- treated as not invaded.")
  if (length(inv_nerves) == 0) { warning("No invaded nerves in ", s); next }
  cat("\n== Section", s, "(patient", pid, "):", length(idx), "spots,", length(inv_nerves), "invaded nerves ==\n")
  
  neigh_list <- contacts_in_section(idx); neigh <- vector("list", ncol(spe)); neigh[idx] <- neigh_list
  is_nerve <- !is.na(nerve_id)
  real_mask <- idx[nerve_id[idx] %in% inv_nerves]
  ni_nerve  <- idx[is_nerve[idx] & !nerve_id[idx] %in% inv_nerves]
  iface_all <- ring_of(idx[is_nerve[idx]], idx, neigh)              # every spot touching any nerve
  eligible  <- setdiff(idx, c(idx[is_nerve[idx]], iface_all))       # background spots
  cat("  invaded-nerve spots:", length(real_mask), "| non-invaded nerve:", length(ni_nerve),
      "| interface (excluded):", length(iface_all), "| background:", length(eligible), "\n")
  nerve_sp <- idx[is_nerve[idx]]
  if_dist <- nn2(coords[nerve_sp, , drop = FALSE], coords[iface_all, , drop = FALSE], k = 1)$nn.dists[, 1] * um_per_pixel
  cat("  interface distance to nearest nerve spot (um): min", round(min(if_dist)), "| median",
      round(median(if_dist)), "| max", round(max(if_dist)), "(expected about 98-205)\n")
  
  # kNN graph on the whole section (global indices), restricted per mask later
  kn <- nn2(coords[idx, , drop = FALSE], coords[idx, , drop = FALSE], k = K_GRAPH + 1)$nn.idx[, -1, drop = FALSE]
  knn_idx <- matrix(NA_integer_, ncol(spe), K_GRAPH); knn_idx[idx, ] <- matrix(idx[kn], ncol = K_GRAPH)
  
  # --- sham masks ---
  spacing <- median(nn2(coords[idx, ], coords[idx, ], k = 2)$nn.dists[, 2])
  snap_tol <- SNAP_TOL_FRAC * spacing; min_reloc <- MIN_RELOC_UM / um_per_pixel
  elig_lgl <- rep(FALSE, ncol(spe)); elig_lgl[eligible] <- TRUE
  place_one <- function(nid, used) {
    sp <- idx[nerve_id[idx] == nid & !is.na(nerve_id[idx])]
    cen <- colMeans(coords[sp, , drop = FALSE])
    # Anchor the shape on one of its own spots (the one nearest its centroid) and
    # move that spot onto another spot: the shift is then spot-to-spot, so the
    # shape lands on the spot grid without distortion or fragmentation.
    anchor <- coords[sp[which.min(colSums((t(coords[sp, , drop = FALSE]) - cen)^2))], ]
    rel <- sweep(coords[sp, , drop = FALSE], 2, anchor)
    for (t in seq_len(MAX_TRIES_SHAPE)) {
      c0 <- coords[sample(eligible, 1), ]
      if (sqrt(sum((c0 - anchor)^2)) < min_reloc) next              # centroid moves by the same vector
      nn <- nn2(coords[idx, , drop = FALSE], sweep(rel, 2, c0, "+"), k = 1)
      ok <- nn$nn.dists[, 1] <= snap_tol
      if (mean(ok) < MIN_SNAP_FRAC) next
      g <- idx[nn$nn.idx[ok, 1]]
      if (anyDuplicated(g) || !all(elig_lgl[g]) || any(used[g])) next
      return(list(spots = g, tries = t, snap = mean(ok)))
    }
    list(spots = NULL, tries = MAX_TRIES_SHAPE, snap = NA)
  }
  shams <- list(); sham_snap <- numeric(0); draws <- 0; failed <- 0; tries_total <- 0
  while (length(shams) < N_SHAM && draws < MAX_REDRAWS) {
    draws <- draws + 1; used <- rep(FALSE, ncol(spe)); spots <- integer(0); ok_all <- TRUE; snaps <- numeric(0)
    for (nid in sample(inv_nerves)) {
      r <- place_one(nid, used); tries_total <- tries_total + r$tries
      if (is.null(r$spots)) { ok_all <- FALSE; break }
      used[r$spots] <- TRUE; spots <- c(spots, r$spots); snaps <- c(snaps, r$snap)
    }
    if (ok_all) { shams[[length(shams) + 1]] <- spots; sham_snap <- c(sham_snap, mean(snaps)) } else failed <- failed + 1
  }
  if (length(shams) < N_SHAM) warning(s, ": only ", length(shams), " sham masks could be placed.")
  placement_qc[[s]] <- data.frame(sample_id = s, patient_id = pid, n_invaded_nerves = length(inv_nerves),
                                  real_mask_spots = length(real_mask), sham_masks = length(shams),
                                  sham_draws = draws, sham_draws_failed = failed, shape_attempts = tries_total,
                                  background_spots = length(eligible), interface_excluded = length(iface_all))
  cat("  sham masks:", length(shams), "of", draws, "draws (", failed, "failed )\n")
  
  iface_inv <- intersect(ring_of(real_mask, idx, neigh), iface_all)   # ring of the invaded nerves
  cat("  interface of invaded nerves:", length(iface_inv), "spots\n")
  
  # Expression of the section, dense, built once (spots x genes)
  Xs  <- t(as.matrix(expr_all[, idx, drop = FALSE])); X2s <- Xs^2
  sec_pos <- integer(ncol(spe)); sec_pos[idx] <- seq_along(idx)
  
  # --- Lee's L: real mask ---
  pop_real <- c(eligible, real_mask)
  W  <- graph_on(knn_idx, pop_real)
  qc_rows <- list(mask_qc(W, pop_real, real_mask, length(intersect(ring_of(real_mask, idx, neigh), iface_all)),
                          s, "real", "nerve"))
  L_real <- lee_L(W, pop_real, real_mask, Xs, X2s, sec_pos)
  
  # one-off check of the formula against spdep::lee (first section, 3 genes)
  if (is.null(lee_check) && requireNamespace("spdep", quietly = TRUE)) {
    lw <- spdep::mat2listw(as(W, "CsparseMatrix"), style = "W", zero.policy = TRUE)
    gg <- which(!is.na(L_real))[c(1:3, length(which(!is.na(L_real))))]
    gg <- unique(c(gg, which.max(L_real)))
    ref <- sapply(gg, function(g) spdep::lee(Xs[sec_pos[pop_real], g], as.numeric(pop_real %in% real_mask), lw,
                                             length(pop_real), zero.policy = TRUE)$L)
    lee_check <- max(abs(ref - L_real[gg]))
    cat("  check vs spdep::lee, max |difference|:", signif(lee_check, 3), "\n")
    if (lee_check > 1e-6) stop("Lee's L formula does not match spdep::lee.")
  }
  
  # --- Lee's L: real interface (nerve spots removed, ring = mask) ---
  pop_if <- c(eligible, iface_inv); W_if <- graph_on(knn_idx, pop_if)
  L_real_if <- lee_L(W_if, pop_if, iface_inv, Xs, X2s, sec_pos)
  qc_rows[[2]] <- mask_qc(W_if, pop_if, iface_inv, length(iface_inv), s, "real", "interface")
  
  # --- Lee's L: sham masks (nerve test: shape = mask, own ring excluded;
  #     interface test: own ring = mask, shape excluded) ---
  L_sham <- L_sham_if <- matrix(NA_real_, length(shams), length(genes), dimnames = list(NULL, genes))
  for (m in seq_along(shams)) {
    sm   <- shams[[m]]
    ring <- intersect(ring_of(sm, idx, neigh), eligible)
    pop  <- c(setdiff(eligible, c(sm, ring)), sm); Wm <- graph_on(knn_idx, pop)
    L_sham[m, ] <- lee_L(Wm, pop, sm, Xs, X2s, sec_pos)
    pop2 <- setdiff(eligible, sm); Wm2 <- graph_on(knn_idx, pop2)   # ring is part of eligible, so it is in pop2
    stopifnot(all(ring %in% pop2), !any(sm %in% pop2))
    L_sham_if[m, ] <- lee_L(Wm2, pop2, ring, Xs, X2s, sec_pos)
    qc_rows[[length(qc_rows) + 1]] <- cbind(mask_qc(Wm, pop, sm, length(ring), s, paste0("sham_", m), "nerve"),
                                            mean_snapped_fraction = sham_snap[m])
    qc_rows[[length(qc_rows) + 1]] <- mask_qc(Wm2, pop2, ring, length(ring), s, paste0("sham_", m), "interface")
    if (m %% 20 == 0) cat("  sham", m, "/", length(shams), "\n")
  }
  per_patient[[s]] <- bind_rows(calibrate(L_real, L_sham, "nerve"), calibrate(L_real_if, L_sham_if, "interface")) %>%
    mutate(sample_id = s, patient_id = pid)
  sham_store[[s]] <- list(nerve = L_sham, interface = L_sham_if)
  mask_qc_all[[s]] <- bind_rows(qc_rows)
  rm(Xs, X2s); gc(verbose = FALSE)
}
per_patient <- bind_rows(per_patient); placement_qc <- bind_rows(placement_qc); mask_qc_all <- bind_rows(mask_qc_all)
saveRDS(list(per_patient = per_patient, sham_store = sham_store), file.path(PROCESSED_DIR, "lee_sham_per_patient.rds"))
write.csv(placement_qc, file.path(TABLES_DIR, "lee_sham_placement_QC.csv"), row.names = FALSE)
write.csv(mask_qc_all,  file.path(TABLES_DIR, "lee_sham_mask_QC.csv"), row.names = FALSE)
cat("\nMask QC (median over shams vs real, nerve test):\n")
print(mask_qc_all %>% filter(test == "nerve") %>% mutate(kind = ifelse(mask == "real", "real", "sham")) %>%
        group_by(sample_id, kind) %>%
        summarise(population = median(population), mask_spots = median(mask_spots), ring_spots = median(ring_spots),
                  isolated_mask_spots = median(isolated_mask_spots), .groups = "drop"))
zero_var <- per_patient %>% filter(is.na(rank_u)) %>% count(test, patient_id, name = "genes_without_variance")
if (nrow(zero_var) > 0) { cat("\nGenes with no variance across shams (not tested in that patient):\n"); print(zero_var) }

# ---- 5. Combine patients, BH per family and test --------------------------
lab <- function(p) ifelse(p %in% names(PATIENT_LABELS), PATIENT_LABELS[p], p)
Ms <- sapply(names(sham_store), function(s) nrow(sham_store[[s]]$nerve))
p_exact <- exact_null(Ms)
cat("\nExact null: ", prod(Ms + 1), " rank combinations (", paste(Ms, collapse = "/"),
    " shams per patient); smallest attainable p = ", signif(1 / prod(Ms + 1), 3), "\n", sep = "")

# Null audit: held-out sham triplets, each ranked among the remaining shams
null_audit <- function(tst) {
  secs <- names(sham_store); P <- length(secs)
  p_audit <- exact_null(Ms - 1)
  fam_de <- genes %in% de_pass_genes
  set.seed(SEED + 2)
  cat("Null audit (", tst, "): ", N_AUDIT, " held-out sham triplets ...\n", sep = "")
  bind_rows(lapply(seq_len(N_AUDIT), function(k) {
    S <- Reduce(`+`, lapply(secs, function(s) {
      L <- sham_store[[s]][[tst]]; h <- sample.int(nrow(L), 1); rest <- L[-h, , drop = FALSE]
      below <- colSums(sweep(rest, 2, L[h, ], "<")); ties <- colSums(sweep(rest, 2, L[h, ], "=="))
      qnorm((below + 0.5 * ties + 0.5) / (nrow(rest) + 1)) })) / sqrt(P)
    pf <- p_audit(S)
    data.frame(test = tst, audit = k, frac_p_below_005 = mean(pf < 0.05, na.rm = TRUE),
               frac_p_below_001 = mean(pf < 0.01, na.rm = TRUE),
               BH_hits_DE_pass = sum(p.adjust(pf[fam_de], "BH") < ALPHA, na.rm = TRUE))
  }))
}

summarise_test <- function(tst) {
  pp <- per_patient %>% filter(test == tst)
  comb <- pp %>% group_by(Gene) %>%
    summarise(n_patients = sum(!is.na(rank_u)), L_real_mean = mean(L_real), delta_L_mean = mean(delta_L),
              rank_score = sum(rank_score) / sqrt(n()),
              n_patients_above_sham_median = sum(rank_u > 0.5),
              n_patients_exceed_q95 = sum(exceeds_q95), .groups = "drop") %>%
    mutate(p_comb = ifelse(n_patients == length(Ms), p_exact(rank_score), NA))
  wide <- pp %>% mutate(P = lab(patient_id)) %>% select(Gene, P, L_real, delta_L, rank_u, z) %>%
    pivot_wider(names_from = P, values_from = c(L_real, delta_L, rank_u, z))
  out <- candidates %>% select(Gene, family) %>%
    left_join(comb, by = "Gene") %>% left_join(wide, by = "Gene") %>%
    group_by(family) %>% mutate(padj = p.adjust(p_comb, "BH")) %>% ungroup() %>%
    mutate(status = case_when(is.na(p_comb) ~ "not_tested",
                              padj < ALPHA & rank_score > 0 ~ "pass", TRUE ~ "fail"),
           consistency = paste0(n_patients_exceed_q95, "/", n_patients)) %>%
    select(-family)
  names(out)[-1] <- paste0(tst, "_", names(out)[-1])
  out
}
result <- candidates %>%
  left_join(summarise_test("nerve"), by = "Gene") %>%
  left_join(summarise_test("interface"), by = "Gene") %>%
  mutate(
    not_tested_reason = case_when(
      !Gene %in% rownames(spe)       ~ "gene not in expression matrix",
      nerve_status == "not_tested"   ~ "no variance in at least one patient's masks",
      TRUE                           ~ NA_character_),
    family_role = ifelse(family == "DE_pass", "primary", "exploratory"),
    category = case_when(
      nerve_status == "not_tested"                           ~ "Not tested",
      nerve_status == "pass" & interface_status == "pass"    ~ "Nerve and interface",
      nerve_status == "pass"                                 ~ "Nerve only",
      interface_status == "pass"                             ~ "Interface only",
      TRUE                                                   ~ "Neither")) %>%
  select(Gene, family, family_role, category, nerve_status, nerve_padj, nerve_L_real_mean, nerve_delta_L_mean,
         nerve_consistency, interface_status, interface_padj, interface_consistency,
         DE_log2FC, DE_padj, everything()) %>%
  arrange(family, factor(nerve_status, c("pass", "fail", "not_tested")), nerve_padj, desc(nerve_delta_L_mean))

audit <- bind_rows(null_audit("nerve"), null_audit("interface"))
write.csv(audit, file.path(TABLES_DIR, "lee_sham_null_audit.csv"), row.names = FALSE)
cat("\n=== Null audit: held-out sham triplets analysed as if real (", N_AUDIT, " per test) ===\n", sep = "")
print(audit %>% group_by(test) %>%
        summarise(mean_frac_p_below_005 = round(mean(frac_p_below_005), 3),
                  mean_frac_p_below_001 = round(mean(frac_p_below_001), 4),
                  audits_with_any_BH_hit = sum(BH_hits_DE_pass > 0),
                  max_BH_hits = max(BH_hits_DE_pass), .groups = "drop"))
cat("Calibrated if fractions are ~0.05 / ~0.01 and BH hits are ~0 (a few audits with 1-2 hits is expected).\n")
result <- result %>%
  mutate(signature_core = family == "DE_pass" & nerve_status == "pass" &
           nerve_n_patients_exceed_q95 == length(Ms) & nerve_n_patients == length(Ms)) %>%
  relocate(signature_core, .after = category)
# Outputs, one file per gene family (downstream scripts read ONLY the DE file
# or the core signature; the LR-only family is exploratory and kept apart):
#   lee_sham_all_candidates.csv     both families together (complete record)
#   lee_sham_DE_genes_scores.csv    DE_pass family: every upregulated gene with its scores
#   lee_sham_LR_only_scores.csv     LR_only family (exploratory, not used for selection)
#   lee_sham_signature_core.csv     core signature (subset of the DE file)
write.csv(result, file.path(TABLES_DIR, "lee_sham_all_candidates.csv"), row.names = FALSE)
de_scores <- result %>% filter(family == "DE_pass")
write.csv(de_scores, file.path(TABLES_DIR, "lee_sham_DE_genes_scores.csv"), row.names = FALSE)
write.csv(result %>% filter(family == "LR_only"), file.path(TABLES_DIR, "lee_sham_LR_only_scores.csv"), row.names = FALSE)
signature <- de_scores %>% filter(nerve_status == "pass")
write.csv(signature %>% filter(signature_core), file.path(TABLES_DIR, "lee_sham_signature_core.csv"), row.names = FALSE)

cat("\n=== Lee's L vs matched shams, exact rank test (BH < ", ALPHA, " within family) ===\n", sep = "")
print(placement_qc)
print(table(result$family, result$category))
cat("\nSignature (DE_pass, nerve test):", nrow(signature), "genes; patients above sham 95th percentile:\n")
print(table(signature$nerve_n_patients_exceed_q95))
cat("Core signature (above sham 95th percentile in every patient):", sum(signature$signature_core), "genes\n")
cat("\nSaved: lee_sham_all_candidates.csv, lee_sham_DE_genes_scores.csv, lee_sham_LR_only_scores.csv,\n",
    " lee_sham_signature_core.csv, lee_sham_placement_QC.csv in\n ", TABLES_DIR, "\n")

