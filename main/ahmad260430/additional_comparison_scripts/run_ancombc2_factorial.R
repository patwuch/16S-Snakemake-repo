#!/usr/bin/env Rscript
# Fit ANCOM-BC2 as a Treatment x Radiation factorial for ahmad260430, with an
# interaction term.
#
# Why this is a standalone script and not part of the Snakemake pipeline: the
# pipeline runs ancombc2 with `group = <column>`, which switches on the
# package's multi-group machinery (global test, mdFDR pairwise contrasts,
# structural zeros). With an interaction in the formula those tests cover the
# group's main effects only, i.e. one radiation state, and they compare levels
# as differences of separately bias-corrected coefficients.
#
# So every fit here is `group = NULL`: per-coefficient inference (LFC/SE/W/p/q
# for each model term, interaction terms included). `--pairwise` covers every
# treatment pair at the radiation reference level, each as a direct coefficient
# of a fit whose reference is one of the pair (factor_pairs() in
# ancombc2_common.R), with one BH pass across the pairs.
#
# Needs ANCOMBC >= 2.14.0: run it in containers/ancombc2-env.sif.
#
# `Group`'s six levels are a 3x2 factorial (Treatment x Radiation) flattened into
# one column; metadata.tsv carries the decomposed `Treatment` and `Radiation`
# columns alongside it.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_ancombc2_factorial.R [--db Silva138|Greengenes2|both] [--pairwise]
#
# The BIOM/taxonomy readers below are copied from scripts/run_ancombc2.R rather
# than sourced: that file is a Snakemake `script:` target and expects `snakemake@`
# bindings at load time. ancombc2_common.R supplies require_ancombc() and
# factor_pairs(); it is sourced first, so this file's own readers and collapse
# policy take precedence over its copies.

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
# Results in tables/, plots in figures/, prose at the root. Routing by what the
# file is keeps the folder navigable however many parameterizations run.
TAB_DIR   <- file.path(OUT_DIR, "tables")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

# canonical_lineage()/write_res(): the shared .res interface the pipeline's
# effect-size barplot and cladogram renderer consumes, so ANCOM-BC2 coefficients
# can be drawn in the same chart format as every other method.
source(file.path(REPO_ROOT, "scripts", "write_res.R"))
LOG_PREFIX <- "[factorial]"
source(file.path(THIS_DIR, "ancombc2_common.R"))

suppressPackageStartupMessages({
  library(biomformat)
  library(ANCOMBC)
})
require_ancombc("2.14.0")

# --- Settings. Mirror config/config_ahmad260430.yaml's ANCOMBC2_* block. ---
MIN_SAMPLE_PRESENCE <- 2
PRV_CUT             <- 0.05
LIB_CUT             <- 1000
P_ADJ_METHOD        <- "BH"
ALPHA               <- 0.05
PSEUDO_SENS         <- TRUE

# Taxonomic resolution is a property of the DATABASE, not of the method: the
# same policy run_lefse_huttenhower.py already applies. Verified on this study's
# own taxonomy - Silva138 assigns no species at all here (0/2325 features, so
# levels 6 and 7 are identical), while Greengenes2 resolves species for 1369/2325
# and gains 41 extra groups at level 7. Collapsing also makes every feature a
# taxonomy path, so the .res tree maps 1:1 instead of hundreds of identically
# annotated ASVs piling onto one node.
COLLAPSE_LEVEL <- c(Silva138 = 6, Greengenes2 = 7)   # 6 = genus, 7 = species

collapse_to_level <- function(counts, tax_map, level) {
  if (level <= 0) return(counts)
  paths <- vapply(rownames(counts), function(t) {
    s <- unname(tax_map[trimws(t)])
    if (length(s) == 0 || is.na(s)) return(NA_character_)
    parts <- trimws(strsplit(s, ";")[[1]])
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) return(NA_character_)
    paste(parts[seq_len(min(level, length(parts)))], collapse = ";")
  }, character(1))
  drop <- is.na(paths)
  if (any(drop))
    message("[factorial] ", sum(drop), " feature(s) had no taxonomy and were dropped")
  rowsum(counts[!drop, , drop = FALSE], group = paths[!drop], reorder = FALSE)
}

# Reference level for each factor. Alphabetical ordering would make "IR" the
# Radiation reference ("IR" < "None"), which silently flips the sign of every
# radiation coefficient, so both are pinned explicitly.
TREATMENT_LEVELS <- c("Control", "IFA", "Inosine")
RADIATION_LEVELS <- c("None", "IR")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
db_arg <- arg_val("--db", "both")
dbs <- if (db_arg == "both") c("Silva138", "Greengenes2") else db_arg

# Day17 is the default, not an option: it is the only timepoint present in all
# three treatment arms. Day25 exists for IFA alone, so including it aliases the
# Day17/Day25 difference into the Treatment coefficients (the balance check below
# reports this). Pass "--timepoint all" to fit every sample anyway, or another
# value to select a different timepoint.
timepoint_filter <- arg_val("--timepoint", "Day17")
if (identical(timepoint_filter, "all")) timepoint_filter <- NULL
out_suffix <- if (is.null(timepoint_filter)) "_all-timepoints" else paste0("_", timepoint_filter)

# The interaction model spends 6 parameters on 24 samples (4 per cell), so a
# taxon seen in only a couple of samples cannot support an interaction
# estimate — its standard error comes back NaN. Raising this trades taxa
# tested for estimates that actually exist. Runs at a non-default value get
# their own output file so prevalence settings can be compared side by side.
MIN_SAMPLE_PRESENCE_DEFAULT <- MIN_SAMPLE_PRESENCE
MIN_SAMPLE_PRESENCE <- as.integer(arg_val("--min-presence", MIN_SAMPLE_PRESENCE))
if (is.na(MIN_SAMPLE_PRESENCE)) stop("--min-presence must be an integer")
if (MIN_SAMPLE_PRESENCE != MIN_SAMPLE_PRESENCE_DEFAULT)
  out_suffix <- paste0(out_suffix, "_minpres", MIN_SAMPLE_PRESENCE)

# Every coefficient is conditional on the other factor's reference level, so
# `RadiationIR` is the radiation effect *within the reference treatment arm* only.
# The radiation effect inside a non-reference arm is beta_Radiation +
# beta_interaction, whose standard error needs the covariance between them and
# so cannot be read off the table. Refitting with that arm as the reference
# makes ancombc2 estimate it directly instead.
treatment_ref <- arg_val("--treatment-reference", TREATMENT_LEVELS[1])
if (!(treatment_ref %in% TREATMENT_LEVELS))
  stop("--treatment-reference must be one of: ", paste(TREATMENT_LEVELS, collapse = ", "))
if (!identical(treatment_ref, TREATMENT_LEVELS[1])) {
  TREATMENT_LEVELS <- c(treatment_ref, setdiff(TREATMENT_LEVELS, treatment_ref))
  out_suffix  <- paste0(out_suffix, "_ref", treatment_ref)
}

# Treatment coefficients are reported at the Radiation reference level, so reading a
# treatment contrast among irradiated mice means making IR the reference. Same
# reparameterization trick as --treatment-reference, on the other factor.
rad_ref <- arg_val("--radiation-reference", RADIATION_LEVELS[1])
if (!(rad_ref %in% RADIATION_LEVELS))
  stop("--radiation-reference must be one of: ", paste(RADIATION_LEVELS, collapse = ", "))
if (!identical(rad_ref, RADIATION_LEVELS[1])) {
  RADIATION_LEVELS <- c(rad_ref, setdiff(RADIATION_LEVELS, rad_ref))
  out_suffix <- paste0(out_suffix, "_radref", rad_ref)
}

# "--treatment IFA" restricts to a single arm and drops Treatment from the model, giving
# a plain within-arm radiation test. Deliberately NOT the same analysis as
# reading RadiationIR off the factorial at that reference level: it pools no
# variance across arms, and re-estimates the sampling-fraction bias from only
# these samples. Useful precisely because it corroborates (or fails to) under
# weaker assumptions.
# "--collapse-level N" overrides the per-database policy for every database;
# 0 keeps ASV-level features (what earlier runs of this script did).
collapse_override <- arg_val("--collapse-level", NULL)
if (!is.null(collapse_override)) {
  collapse_override <- as.integer(collapse_override)
  if (is.na(collapse_override)) stop("--collapse-level must be an integer")
  out_suffix <- paste0(out_suffix, "_lvl", collapse_override)
}

# "--pairwise": every treatment pair at the radiation reference level. A pair is
# the direct coefficient of a fit whose reference is one of its two arms (IFA
# and Inosine against Control from the Control-reference fit, Inosine against
# IFA from the IFA-reference fit), and one BH pass covers all pairs' p-values.
pairwise_mode <- "--pairwise" %in% args
if (pairwise_mode && "--treatment-reference" %in% args)
  stop("--pairwise fits each treatment reference it needs; drop --treatment-reference")
if (pairwise_mode) out_suffix <- paste0(out_suffix, "_pairwise")

treatment_subset <- arg_val("--treatment", NULL)
if (!is.null(treatment_subset)) {
  if (!(treatment_subset %in% TREATMENT_LEVELS))
    stop("--treatment must be one of: ", paste(TREATMENT_LEVELS, collapse = ", "))
  out_suffix <- paste0(out_suffix, "_treatment", treatment_subset)
  if (pairwise_mode) stop("--pairwise needs every treatment arm; drop --treatment")
}

# ---------- BIOM format helpers ----------
.is_hdf5_biom <- function(path) {
  magic <- tryCatch(readBin(path, raw(), n = 4), error = function(e) raw(0))
  length(magic) >= 4 && identical(magic[1:4], as.raw(c(0x89, 0x48, 0x44, 0x46)))
}

.read_hdf5_biom <- function(path) {
  if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("Package 'rhdf5' is required for HDF5 BIOM files.")
  obs_ids  <- as.character(rhdf5::h5read(path, "observation/ids"))
  samp_ids <- as.character(rhdf5::h5read(path, "sample/ids"))
  data_v   <- rhdf5::h5read(path, "observation/matrix/data")
  indices  <- as.integer(rhdf5::h5read(path, "observation/matrix/indices")) + 1L
  indptr   <- as.integer(rhdf5::h5read(path, "observation/matrix/indptr"))  + 1L
  mat <- matrix(0L, nrow = length(obs_ids), ncol = length(samp_ids),
                dimnames = list(obs_ids, samp_ids))
  for (i in seq_along(obs_ids)) {
    col_idx <- indices[indptr[i]:(indptr[i + 1L] - 1L)]
    mat[i, col_idx] <- data_v[indptr[i]:(indptr[i + 1L] - 1L)]
  }
  mat
}

.read_biom_auto <- function(path) {
  if (.is_hdf5_biom(path)) .read_hdf5_biom(path)
  else as.matrix(biomformat::biom_data(biomformat::read_biom(path)))
}

# ---------- Taxonomy ----------
read_taxonomy_map <- function(path) {
  tax_df <- read.delim(path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
  col <- intersect(c("Taxon", "taxonomy"), colnames(tax_df))[1]
  if (is.na(col)) return(setNames(character(0), character(0)))
  setNames(as.character(tax_df[[col]]), trimws(rownames(tax_df)))
}

# ---------- Main ----------
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

biom_path <- file.path(STUDY_DIR, "tables", "study-seqs.biom")
meta_path <- file.path(STUDY_DIR, "metadata.tsv")

message("[factorial] Reading BIOM table: ", biom_path)
otu_all <- .read_biom_auto(biom_path)

meta_all <- read.delim(meta_path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
for (col in c("Treatment", "Radiation")) {
  if (!(col %in% colnames(meta_all)))
    stop("metadata.tsv has no '", col, "' column. It is derived from 'Group'; ",
         "see the Treatment/Radiation columns added for this analysis.")
}

common <- intersect(colnames(otu_all), rownames(meta_all))
if (length(common) < 2) stop("Too few overlapping samples between BIOM table and metadata.")
otu_all  <- otu_all[, common, drop = FALSE]
meta_all <- meta_all[common, , drop = FALSE]

meta_all$Treatment <- factor(meta_all$Treatment, levels = TREATMENT_LEVELS)
meta_all$Radiation <- factor(meta_all$Radiation, levels = RADIATION_LEVELS)
if (anyNA(meta_all$Treatment) || anyNA(meta_all$Radiation))
  stop("Treatment/Radiation contain values outside the declared levels: ",
       paste(c(setdiff(unique(meta_all$Treatment), TREATMENT_LEVELS),
               setdiff(unique(meta_all$Radiation), RADIATION_LEVELS)), collapse = ", "))

if (!is.null(timepoint_filter)) {
  if (!("Timepoint" %in% colnames(meta_all))) stop("--timepoint given but metadata has no 'Timepoint' column.")
  keep_s <- meta_all$Timepoint == timepoint_filter
  if (!any(keep_s)) stop("No samples with Timepoint == '", timepoint_filter, "'.")
  message("[factorial] Restricting to Timepoint == ", timepoint_filter,
          " (", sum(keep_s), " of ", nrow(meta_all), " samples)")
  meta_all <- meta_all[keep_s, , drop = FALSE]
  otu_all  <- otu_all[, rownames(meta_all), drop = FALSE]
  meta_all$Treatment <- droplevels(meta_all$Treatment)
  meta_all$Radiation <- droplevels(meta_all$Radiation)
}

if (!is.null(treatment_subset)) {
  keep_d <- meta_all$Treatment == treatment_subset
  if (!any(keep_d)) stop("No samples with Treatment == '", treatment_subset, "'.")
  message("[factorial] Restricting to Treatment == ", treatment_subset,
          " (", sum(keep_d), " of ", nrow(meta_all), " samples)")
  meta_all <- meta_all[keep_d, , drop = FALSE]
  otu_all  <- otu_all[, rownames(meta_all), drop = FALSE]
  meta_all$Treatment <- droplevels(meta_all$Treatment)
}

message("[factorial] Design (samples per cell):")
print(table(meta_all$Treatment, meta_all$Radiation))

# Presentation-only metadata columns. These encode plot colour and display
# order for the Group axis, so they are deterministic functions of Group and
# would always trip the balance check below without ever being real
# confounders. Excluded so the check reports only variables that could
# actually carry a biological effect.
IGNORE_COLUMNS <- c("PrimaryColor", "PrimaryOrder")

# A model term is only interpretable if it is not aliased with something absent
# from the formula. Any other categorical metadata column that is unevenly
# distributed across Treatment or Radiation is a confounder, and its effect will be
# attributed to whichever model term it tracks.
for (nuisance in setdiff(colnames(meta_all), c("Treatment", "Radiation", IGNORE_COLUMNS))) {
  vals <- meta_all[[nuisance]]
  if (is.numeric(vals) || length(unique(vals)) < 2 || length(unique(vals)) > 8) next
  for (fac in c("Treatment", "Radiation")) {
    tab <- table(meta_all[[fac]], vals)
    # Unbalanced == some level of `fac` never sees some level of `nuisance`.
    if (any(tab == 0)) {
      message("[factorial] WARNING: '", nuisance, "' is unbalanced across ", fac,
              " — its effect is aliased into the ", fac, " coefficients:")
      print(tab)
    }
  }
}

otu_raw <- otu_all[rowSums(otu_all, na.rm = TRUE) > 0, , drop = FALSE]

# With a single treatment arm there is no Treatment term to estimate and no interaction.
FIX_FORMULA <- if (nlevels(meta_all$Treatment) < 2) "Radiation" else "Treatment * Radiation"

# The fit is per-database now, not shared: each database collapses the same
# reads onto a different feature space, so they are genuinely different models
# rather than one model with two annotations.
for (db in dbs) {
  tax_path <- file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(db, "_taxonomy.tsv"))
  if (!file.exists(tax_path)) {
    message("[factorial] Skipping ", db, ": no taxonomy at ", tax_path)
    next
  }
  tax_map <- read_taxonomy_map(tax_path)
  level <- if (!is.null(collapse_override)) collapse_override else
           if (db %in% names(COLLAPSE_LEVEL)) COLLAPSE_LEVEL[[db]] else 0L

  counts_db <- collapse_to_level(otu_raw, tax_map, level)
  n_present <- rowSums(counts_db > 0)
  # Every cell of the design must see the taxon at least once. ANCOM-BC2
  # treats a zero as a missing log count, so a taxon absent from a whole cell
  # leaves that cell's coefficient unidentifiable, and the fit can return
  # plausible-looking numbers for it instead of an error. Not optional, so no
  # flag. With six cells this implies presence in >= 6 of 24 samples.
  cell <- if (nlevels(meta_all$Treatment) < 2) droplevels(meta_all$Radiation) else
          interaction(meta_all$Treatment, meta_all$Radiation, drop = TRUE)
  in_every_cell <- rowSums(sapply(levels(cell), function(k)
    rowSums(counts_db[, cell == k, drop = FALSE] > 0) == 0)) == 0
  keep_db <- n_present >= MIN_SAMPLE_PRESENCE & in_every_cell
  # Taxa common enough to pass the prevalence filter but absent from a whole
  # cell are recorded with their per-cell presence rather than dropped silently.
  excl <- n_present >= MIN_SAMPLE_PRESENCE & !in_every_cell
  per_cell <- sapply(levels(cell), function(k) rowSums(counts_db[excl, cell == k, drop = FALSE] > 0))
  per_cell <- matrix(per_cell, nrow = sum(excl), dimnames = list(NULL, levels(cell)))
  colnames(per_cell) <- sapply(levels(cell), function(k) {
    i <- which(cell == k)[1]; sprintf("%s (n=%d)", meta_all$Group[i], sum(cell == k)) })
  excl_path <- file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", db, out_suffix, "_untestable_absent_from_a_cell.tsv"))
  write.table(data.frame(taxon = rownames(counts_db)[excl], n_samples = n_present[excl], per_cell,
                         check.names = FALSE),
              excl_path, sep = "\t", quote = FALSE, row.names = FALSE)
  counts_db <- counts_db[keep_db, , drop = FALSE]
  message(sprintf("[factorial] %s: collapse level %d -> %d feature(s); %d in >= %d samples, %d of those in every cell are tested; %d untestable written",
                  db, level, length(keep_db), sum(n_present >= MIN_SAMPLE_PRESENCE), MIN_SAMPLE_PRESENCE,
                  nrow(counts_db), sum(excl)))
  if (nrow(counts_db) == 0) { message("[factorial] nothing left for ", db); next }

  # The empty-cell failure's only symptom is a vector-recycling warning; with
  # the per-cell rule above it must not appear, so it is an error here.
  fit <- function(treatment_levels) {
    m <- meta_all
    m$Treatment <- factor(m$Treatment, levels = treatment_levels)
    message("[factorial] Fitting: ", FIX_FORMULA, "  (group = NULL; reference ",
            treatment_levels[1], ", ", levels(m$Radiation)[1], ")")
    res <- withCallingHandlers(tryCatch(
      ANCOMBC::ancombc2(data = counts_db, taxa_are_rows = TRUE, meta_data = m,
               fix_formula = FIX_FORMULA, group = NULL,
               prv_cut = PRV_CUT, lib_cut = LIB_CUT, p_adj_method = P_ADJ_METHOD,
               alpha = ALPHA, pseudo_sens = PSEUDO_SENS, struc_zero = FALSE,
               neg_lb = FALSE, global = FALSE, pairwise = FALSE, dunnet = FALSE,
               trend = FALSE, verbose = FALSE),
      error = function(e) stop("ancombc2 failed for ", db, ": ", conditionMessage(e))),
      warning = function(w) {
        if (grepl("multiple of", conditionMessage(w)))
          stop("ancombc2 recycled a vector across columns (rank-deficient cell?): ", conditionMessage(w))
      })
    as.data.frame(res$res)
  }

  res_df <- fit(levels(meta_all$Treatment))
  q_cols <- setdiff(grep("^q_", colnames(res_df), value = TRUE), "q_(Intercept)")

  out <- res_df
  # Collapsed features ARE taxonomy paths, so the id is already the label.
  out$Taxon <- if (level > 0) as.character(out$taxon) else
               unname(tax_map[trimws(as.character(out$taxon))])
  out <- out[, c("taxon", "Taxon", setdiff(colnames(out), c("taxon", "Taxon")))]

  out_path <- file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", db, out_suffix, ".tsv"))
  write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[factorial] Wrote ", out_path)
  if (pairwise_mode) {
    pair_df <- factor_pairs(res_df, levels(meta_all$Treatment), "Treatment", fit, P_ADJ_METHOD, ALPHA)
    pair_path <- sub("\\.tsv$", "_pairs.tsv", out_path)
    write.table(pair_df, pair_path, sep = "\t", quote = FALSE, row.names = FALSE)
    message("[factorial] Wrote ", pair_path)
  }

  # ---- one .res per coefficient, for the shared effect-size/cladogram renderer ----
  # Each main-effect coefficient is a two-class contrast conditional on the
  # other factor's reference level, so the .res is written over exactly that
  # subpopulation: abundances shown then match the comparison being drawn
  # rather than pooling cells the coefficient never looked at. Interaction
  # terms get no .res -- a difference of differences has no "enriched class".
  res_taxa <- as.character(res_df$taxon)
  counts_res <- counts_db[res_taxa[res_taxa %in% rownames(counts_db)], , drop = FALSE]
  res_taxa <- rownames(counts_res)
  lineage_map <- lapply(res_taxa, function(t) {
    s <- if (level > 0) t else unname(tax_map[trimws(t)])
    if (length(s) == 0 || is.na(s)) s <- ""
    canonical_lineage(paste(trimws(strsplit(s, ";")[[1]]), collapse = "|"))
  })
  keys <- vapply(lineage_map, paste, character(1), collapse = ".")
  row_of <- match(res_taxa, as.character(res_df$taxon))

  cell_of <- function(d, r) {
    g <- unique(meta_all$Group[meta_all$Treatment == d & meta_all$Radiation == r])
    if (length(g) != 1) paste(d, r, sep = "/") else g
  }
  treatment_ref <- levels(meta_all$Treatment)[1]; rad_ref <- levels(meta_all$Radiation)[1]
  specs <- list()
  # In pairwise mode the treatment contrasts come from pair_df (joint BH) below,
  # not from the per-coefficient BH columns.
  if (!pairwise_mode)
  for (lvl in setdiff(levels(meta_all$Treatment), treatment_ref))
    specs[[paste0("Treatment", lvl)]] <- list(keep = meta_all$Radiation == rad_ref,
      pos = cell_of(lvl, rad_ref), neg = cell_of(treatment_ref, rad_ref))
  for (lvl in setdiff(levels(meta_all$Radiation), rad_ref))
    specs[[paste0("Radiation", lvl)]] <- list(keep = meta_all$Treatment == treatment_ref,
      pos = cell_of(treatment_ref, lvl), neg = cell_of(treatment_ref, rad_ref))

  for (coef in names(specs)) {
    sp <- specs[[coef]]
    qc <- paste0("q_", coef); lc <- paste0("lfc_", coef); sc <- paste0("passed_ss_", coef)
    if (!all(c(qc, lc) %in% colnames(res_df))) next
    q <- res_df[[qc]][row_of]; l <- res_df[[lc]][row_of]
    ss <- if (sc %in% colnames(res_df)) res_df[[sc]][row_of] else TRUE
    sig <- !is.na(q) & !is.na(l) & q < ALPHA & (is.na(ss) | ss %in% c(TRUE, "TRUE"))

    # ancombc2 fits lm() on log(count) with R's log(), so lfc is a NATURAL-log
    # fold change. Divide by log(2) to publish the log2 fold change the
    # effect-size axis names. Ranking and q-values are unaffected; only scale.
    markers <- if (any(sig)) data.frame(
      lineage = keys[sig], class = ifelse(l[sig] > 0, sp$pos, sp$neg),
      lda = abs(l[sig]) / log(2), pvalue = q[sig], stringsAsFactors = FALSE) else NULL

    samples <- rownames(meta_all)[sp$keep]
    samples <- samples[samples %in% colnames(counts_res)]
    cs <- counts_res[, samples, drop = FALSE]
    keep_depth <- colSums(cs, na.rm = TRUE) > 0
    cs <- cs[, keep_depth, drop = FALSE]
    rp <- file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", db, out_suffix, "_", coef, ".res"))
    write_res(cs, as.character(meta_all[colnames(cs), "Group"]), lineage_map, markers, rp)
    message(sprintf("[factorial] Wrote %s  (%d marker(s), %s vs %s, n=%d)",
                    basename(rp), if (is.null(markers)) 0L else nrow(markers),
                    sp$pos, sp$neg, ncol(cs)))
  }

  # ---- pairwise contrasts: significance is the jointly corrected q (diff_) ----
  # pair_df names a contrast "TreatmentB_TreatmentA" for B - A. The .res file
  # names it the way a reference-level refit would ("_refA_TreatmentB"), so the
  # plotting script titles it without knowing about pairwise mode.
  if (pairwise_mode) {
    prow <- match(res_taxa, as.character(pair_df$taxon))
    for (lc in grep("^lfc_", colnames(pair_df), value = TRUE)) {
      k <- sub("^lfc_", "", lc)
      lv <- sub("^Treatment", "", strsplit(k, "_")[[1]])
      b <- lv[1]; a <- if (length(lv) > 1) lv[2] else treatment_ref
      l <- pair_df[[lc]][prow]; d <- pair_df[[paste0("diff_", k)]][prow]
      ss <- pair_df[[paste0("passed_ss_", k)]][prow]
      sig <- !is.na(l) & d %in% c(TRUE, "TRUE") & (is.na(ss) | ss %in% c(TRUE, "TRUE"))
      pos <- cell_of(b, rad_ref); neg <- cell_of(a, rad_ref)
      markers <- if (any(sig)) data.frame(lineage = keys[sig], class = ifelse(l[sig] > 0, pos, neg),
        lda = abs(l[sig]) / log(2), pvalue = pair_df[[paste0("q_", k)]][prow][sig],
        stringsAsFactors = FALSE) else NULL
      samples <- rownames(meta_all)[meta_all$Radiation == rad_ref & meta_all$Treatment %in% c(a, b)]
      cs <- counts_res[, samples[samples %in% colnames(counts_res)], drop = FALSE]
      cs <- cs[, colSums(cs, na.rm = TRUE) > 0, drop = FALSE]
      tag <- if (a == treatment_ref) "" else paste0("_ref", a)
      rp <- file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", db, out_suffix, tag, "_Treatment", b, ".res"))
      write_res(cs, as.character(meta_all[colnames(cs), "Group"]), lineage_map, markers, rp)
      message(sprintf("[factorial] Wrote %s  (%d marker(s), joint BH, %s vs %s, n=%d)", basename(rp),
                      if (is.null(markers)) 0L else nrow(markers), pos, neg, ncol(cs)))
    }
  }

  # `passed_ss` is ancombc2's own pseudo-count sensitivity analysis: FALSE means
  # the call flips depending on the pseudo-count used, i.e. it is an artifact of
  # zero-handling rather than a real signal. Counting raw q < alpha without it
  # massively overstates the result, so the sensitivity-passing count is the
  # headline and the raw count is shown only for contrast.
  message("[factorial] Significant taxa per coefficient, ", db,
          " (q < ", ALPHA, " AND passed_ss | raw q < ", ALPHA, "):")
  for (qc in q_cols) {
    coef <- sub("^q_", "", qc)
    sig  <- !is.na(out[[qc]]) & out[[qc]] < ALPHA
    ss   <- paste0("passed_ss_", coef)
    n_ss <- if (ss %in% colnames(out)) sum(sig & out[[ss]] %in% c(TRUE, "TRUE")) else NA_integer_
    message(sprintf("    %-34s %4s | %4d", coef, n_ss, sum(sig)))
  }
}
