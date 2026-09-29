#!/usr/bin/env Rscript
# Readers, the per-database collapse policy, and the pairwise helper, shared by
# this study's project-local ANCOM-BC2 scripts.
#
# These were copied out of scripts/run_ancombc2.R once already, into
# run_ancombc2_factorial.R, because that file is a Snakemake `script:` target and
# expects `snakemake@` bindings at load time, so it cannot be sourced. Rather
# than copy them a third time for the Day17-vs-Day25 fit they live here.
# run_ancombc2_factorial.R still carries its own inline copy and should be
# migrated onto this file; doing so regenerates every figure it produced, so it
# wants a deliberate run rather than a drive-by edit.
#
# Sourcing this file defines names and reads nothing. `LOG_PREFIX` may be set by
# the caller before sourcing to tag the messages.

if (!exists("LOG_PREFIX")) LOG_PREFIX <- "[ancombc2]"
.say <- function(...) message(LOG_PREFIX, " ", ...)

# Taxonomic resolution is a property of the DATABASE, not of the method: the
# same policy run_lefse_huttenhower.py applies. Verified on this study's own
# taxonomy, Silva138 assigns no species at all here (0/2325 features, so levels 6
# and 7 are identical) while Greengenes2 resolves species for 1369/2325.
# Collapsing also makes every feature a taxonomy path, so the .res tree maps 1:1
# instead of hundreds of identically annotated ASVs piling onto one node.
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
  if (any(drop)) .say(sum(drop), " feature(s) had no taxonomy and were dropped")
  rowsum(counts[!drop, , drop = FALSE], group = paths[!drop], reorder = FALSE)
}

read_taxonomy_map <- function(path) {
  tax_df <- read.delim(path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
  col <- intersect(c("Taxon", "taxonomy"), colnames(tax_df))[1]
  if (is.na(col)) return(setNames(character(0), character(0)))
  setNames(as.character(tax_df[[col]]), trimws(rownames(tax_df)))
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

read_biom_auto <- function(path) {
  if (.is_hdf5_biom(path)) .read_hdf5_biom(path)
  else as.matrix(biomformat::biom_data(biomformat::read_biom(path)))
}

# ---------- ANCOMBC version ----------
# These scripts rely on ANCOMBC >= 2.14.0: interaction formulas written with `*`,
# a plain count matrix as input, and the refit-based pseudo-count sensitivity
# analysis. Older releases run but silently give different answers.
require_ancombc <- function(min_version) {
  v <- utils::packageVersion("ANCOMBC")
  if (v < min_version)
    stop("ANCOMBC ", v, " is installed; these scripts need >= ", min_version,
         ". Run them in containers/ancombc2-env.sif.")
  .say("ANCOMBC ", v)
}

# ---------- Every pair of a factor's levels, each as a direct coefficient ----------
# fit_ref(levels) refits the model with levels[1] as the reference and returns
# ancombc2's `res` table; `first` is that table for levels as given. The pair
# (a, b), a earlier in `levels`, is read off the fit whose reference is a, so
# every pair is one model coefficient rather than a difference of two: ANCOM-BC2
# bias-corrects each coefficient separately, and a difference of two corrected
# coefficients does not equal the directly fitted pair. One p.adjust pass
# covers every pair's p-values. Columns are named as ancombc2's own pairwise
# output names them: "<F>b" against the first level, "<F>b_<F>a" otherwise.
factor_pairs <- function(first, levels, factor_name, fit_ref, p_adj_method, alpha) {
  fits <- list(first)
  names(fits) <- levels[1]
  for (r in levels[-c(1, length(levels))]) fits[[r]] <- fit_ref(c(r, setdiff(levels, r)))
  pair_df <- data.frame(taxon = as.character(first$taxon))
  keys <- character(0)
  for (i in seq_len(length(levels) - 1)) for (b in levels[(i + 1):length(levels)]) {
    a <- levels[i]; d <- fits[[a]]; coef <- paste0(factor_name, b)
    k <- if (i == 1) coef else paste0(coef, "_", factor_name, a)
    r <- match(pair_df$taxon, as.character(d$taxon))
    for (s in c("lfc", "se", "W", "p", "passed_ss")) pair_df[[paste0(s, "_", k)]] <- d[[paste0(s, "_", coef)]][r]
    keys <- c(keys, k)
  }
  q <- matrix(p.adjust(unlist(pair_df[paste0("p_", keys)]), method = p_adj_method), ncol = length(keys))
  for (j in seq_along(keys)) {
    k <- keys[j]
    pair_df[[paste0("q_", k)]] <- q[, j]
    pair_df[[paste0("diff_", k)]] <- !is.na(q[, j]) & q[, j] <= alpha
    pair_df[[paste0("diff_robust_", k)]] <- pair_df[[paste0("diff_", k)]] & pair_df[[paste0("passed_ss_", k)]] %in% TRUE
  }
  pair_df[c("taxon", unlist(lapply(c("lfc", "se", "W", "p", "q", "diff", "passed_ss", "diff_robust"),
                                   function(s) paste0(s, "_", keys))))]
}
