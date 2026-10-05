#!/usr/bin/env Rscript
# Fit ANCOM-BC2 as a Diet x Condition factorial for tisha251224, with an
# interaction term. Greengenes2 only.
#
# Ported from ahmad260430/additional_comparison_scripts/run_ancombc2_factorial.R;
# that file's header explains why this lives outside the Snakemake pipeline.
# In short: the pipeline runs ancombc2 with `group = <column>`, and its
# multi-group tests cover main effects only when the formula has an interaction.
# Here every fit is `group = NULL` (per-coefficient inference). `--pairwise`
# covers every pair of diets, each as a direct coefficient of a fit whose
# reference is one of the pair.
#
# Needs ANCOMBC >= 2.14.0: run it in containers/ancombc2-env.sif.
#
# Design: 3 diets (Cellulose, Inulin, RS) x 2 conditions (Control = not
# irradiated, IR); 19 samples, 2-4 per cell. Timepoint is ignored by design.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_ancombc2_factorial.R [--min-presence N] [--diet-reference Inulin]
#       [--condition-reference IR] [--pairwise]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Condition-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
META_PATH <- file.path(OUT_DIR, "metadata_factorial.tsv")   # build_factorial_metadata.py
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PREFIX <- "[factorial]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
# canonical_lineage()/write_res(): the shared .res interface the pipeline's
# effect-size barplot and cladogram renderer consumes.
source(file.path(REPO_ROOT, "scripts", "write_res.R"))

suppressPackageStartupMessages({
  library(biomformat)
  library(ANCOMBC)
})
require_ancombc("2.14.0")

DB           <- "Greengenes2"
PRV_CUT      <- 0.05
LIB_CUT      <- 1000
P_ADJ_METHOD <- "BH"
ALPHA        <- 0.05
PSEUDO_SENS  <- TRUE

# Reference levels pinned explicitly: alphabetical order happens to agree here,
# but a relabel must never silently flip the sign of a coefficient.
DIET_LEVELS      <- c("Cellulose", "Inulin", "RS")
CONDITION_LEVELS <- c("Control", "IR")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
out_suffix <- ""

# Two prevalence rules, both applied:
#  * every Diet x Condition cell must see the taxon in at least one sample.
#    ANCOM-BC2 treats zeros as missing log values, so a taxon absent from
#    a whole cell leaves the interaction unidentifiable and the fit returns
#    plausible-looking nonsense instead of an error. Not optional, so no flag.
#    With six cells this alone implies presence in >= 6 of 19 samples.
#  * --min-presence N: the taxon must also be seen in >= N samples overall.
#    Runs at a non-default value get their own output file.
MIN_SAMPLE_PRESENCE <- as.integer(arg_val("--min-presence", "6"))
if (is.na(MIN_SAMPLE_PRESENCE)) stop("--min-presence must be an integer")
if (MIN_SAMPLE_PRESENCE != 6L) out_suffix <- paste0(out_suffix, "_minpres", MIN_SAMPLE_PRESENCE)

# Every coefficient is conditional on the other factor's reference level, so
# `ConditionIR` is the IR effect within the reference diet only. Refitting with
# another diet as the reference makes ancombc2 estimate IR within that diet
# directly, with its own standard error.
diet_ref <- arg_val("--diet-reference", DIET_LEVELS[1])
if (!(diet_ref %in% DIET_LEVELS)) stop("--diet-reference must be one of: ", paste(DIET_LEVELS, collapse = ", "))
if (!identical(diet_ref, DIET_LEVELS[1])) {
  DIET_LEVELS <- c(diet_ref, setdiff(DIET_LEVELS, diet_ref))
  out_suffix  <- paste0(out_suffix, "_ref", diet_ref)
}
cond_ref <- arg_val("--condition-reference", CONDITION_LEVELS[1])
if (!(cond_ref %in% CONDITION_LEVELS)) stop("--condition-reference must be one of: ", paste(CONDITION_LEVELS, collapse = ", "))
if (!identical(cond_ref, CONDITION_LEVELS[1])) {
  CONDITION_LEVELS <- c(cond_ref, setdiff(CONDITION_LEVELS, cond_ref))
  out_suffix <- paste0(out_suffix, "_condref", cond_ref)
}

# "--pairwise": every pair of diets at the condition reference level. A pair is
# the direct coefficient of a fit whose reference is one of its two diets
# (Inulin and RS against Cellulose from the Cellulose-reference fit, RS against
# Inulin from the Inulin-reference fit), and one Benjamini-Hochberg pass covers
# all pairs' p-values together. Writes <stem>_pairs.tsv.
pairwise_mode <- "--pairwise" %in% args
if (pairwise_mode && "--diet-reference" %in% args)
  stop("--pairwise fits each diet reference it needs; drop --diet-reference")
if (pairwise_mode) out_suffix <- paste0(out_suffix, "_pairwise")

# ---------- Main ----------
otu_all <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
meta <- read.delim(META_PATH, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
common <- intersect(colnames(otu_all), rownames(meta))
if (length(common) != nrow(meta))
  .say("WARNING: ", nrow(meta) - length(common), " metadata sample(s) missing from the BIOM table")
otu_all <- otu_all[, common, drop = FALSE]
meta    <- meta[common, , drop = FALSE]
meta$Diet      <- factor(meta$Diet, levels = DIET_LEVELS)
meta$Condition <- factor(meta$Condition, levels = CONDITION_LEVELS)
if (anyNA(meta$Diet) || anyNA(meta$Condition)) stop("Diet/Condition outside the declared levels")

.say("Design (samples per cell):")
print(table(meta$Diet, meta$Condition))

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(DB, "_taxonomy.tsv")))
counts <- collapse_to_level(otu_all[rowSums(otu_all) > 0, , drop = FALSE], tax_map, COLLAPSE_LEVEL)

cell <- interaction(meta$Diet, meta$Condition, drop = TRUE)
in_every_cell <- sapply(levels(cell), function(k) rowSums(counts[, cell == k, drop = FALSE] > 0) > 0)
n_present <- rowSums(counts > 0)
keep <- rowSums(!in_every_cell) == 0 & n_present >= MIN_SAMPLE_PRESENCE
.say(sprintf("%s: %d features after collapse; %d seen in every cell; %d of those in >= %d samples are tested",
             DB, nrow(counts), sum(rowSums(!in_every_cell) == 0), sum(keep), MIN_SAMPLE_PRESENCE))

# Taxa common enough to test but absent from at least one whole cell cannot
# enter the model. Where they are missing is a result in its own right, so
# record it rather than dropping them silently.
if (out_suffix == "") {
  excl <- n_present >= MIN_SAMPLE_PRESENCE & !keep
  per_cell <- sapply(levels(cell), function(k) rowSums(counts[excl, cell == k, drop = FALSE] > 0))
  colnames(per_cell) <- sapply(levels(cell), function(k) {
    i <- which(cell == k)[1]; sprintf("%s (n=%d)", meta$Group[i], sum(cell == k)) })
  excl_df <- data.frame(taxon = rownames(counts)[excl], n_samples = n_present[excl],
                        per_cell, check.names = FALSE)
  # Named for the presence analyses that read it, not for a DA method.
  excl_path <- file.path(TAB_DIR, paste0("presence_factorial_", DB, "_untestable_absent_from_a_cell.tsv"))
  write.table(excl_df, excl_path, sep = "\t", quote = FALSE, row.names = FALSE)
  .say("Wrote ", excl_path, " (", nrow(excl_df), " taxa)")
}
counts <- counts[keep, , drop = FALSE]
# The tested taxa after collapse and both filters, so LinDA fits exactly these
# features. "--export-only" stops here.
write.table(data.frame(taxon = rownames(counts), counts, check.names = FALSE),
            file.path(TAB_DIR, paste0("tested_counts_", DB, out_suffix, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
# The same taxa as a plain list: with the untestable table, the presence
# analyses' fixed taxon set (every taxon in >= --min-presence samples).
if (out_suffix == "")
  write.table(data.frame(taxon = rownames(counts)), file.path(TAB_DIR, paste0("presence_factorial_", DB, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
if ("--export-only" %in% args) quit(save = "no")

FIX_FORMULA <- "Diet * Condition"
# The empty-cell failure's only symptom is a recycling warning; with the
# per-cell rule above it must not appear, so turn it into an error.
fit <- function(diet_levels) {
  m <- meta
  m$Diet <- factor(m$Diet, levels = diet_levels)
  .say("Fitting: ", FIX_FORMULA, "  (group = NULL; reference ", diet_levels[1], ", ",
       levels(m$Condition)[1], ")")
  res <- withCallingHandlers(
    ANCOMBC::ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = m,
             fix_formula = FIX_FORMULA, group = NULL,
             prv_cut = PRV_CUT, lib_cut = LIB_CUT, p_adj_method = P_ADJ_METHOD,
             alpha = ALPHA, pseudo_sens = PSEUDO_SENS, struc_zero = FALSE,
             neg_lb = FALSE, global = FALSE, pairwise = FALSE, dunnet = FALSE,
             trend = FALSE, verbose = FALSE),
    warning = function(w) {
      if (grepl("multiple of", conditionMessage(w)))
        stop("ancombc2 recycled a vector across columns (rank-deficient cell?): ", conditionMessage(w))
    })
  as.data.frame(res$res)
}

res_df <- fit(levels(meta$Diet))
out <- res_df
out$Taxon <- as.character(out$taxon)
out <- out[, c("taxon", "Taxon", setdiff(colnames(out), c("taxon", "Taxon")))]
stem <- paste0("ANCOMBC2_factorial_", DB, out_suffix)
out_path <- file.path(TAB_DIR, paste0(stem, ".tsv"))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
.say("Wrote ", out_path)
if (pairwise_mode) {
  pair_df <- factor_pairs(res_df, levels(meta$Diet), "Diet", fit, P_ADJ_METHOD, ALPHA)
  pair_path <- file.path(TAB_DIR, paste0(stem, "_pairs.tsv"))
  write.table(pair_df, pair_path, sep = "\t", quote = FALSE, row.names = FALSE)
  .say("Wrote ", pair_path)
}

# ---- one .res per main-effect coefficient, for the shared renderer ----
# Each is a two-class contrast conditional on the other factor's reference
# level, so the .res carries exactly that subpopulation. Interaction terms get
# no .res: a difference of differences has no enriched class.
res_taxa <- intersect(as.character(res_df$taxon), rownames(counts))
counts_res <- counts[res_taxa, , drop = FALSE]
lineage_map <- lapply(res_taxa, function(t) canonical_lineage(paste(trimws(strsplit(t, ";")[[1]]), collapse = "|")))
keys <- vapply(lineage_map, paste, character(1), collapse = ".")
row_of <- match(res_taxa, as.character(res_df$taxon))

cell_of <- function(d, r) unique(meta$Group[meta$Diet == d & meta$Condition == r])
diet_ref <- levels(meta$Diet)[1]; cond_ref <- levels(meta$Condition)[1]

write_one <- function(rp, keep_samples, markers, pos, neg, note = "") {
  cs <- counts_res[, rownames(meta)[keep_samples], drop = FALSE]
  cs <- cs[, colSums(cs) > 0, drop = FALSE]
  write_res(cs, as.character(meta[colnames(cs), "Group"]), lineage_map, markers, rp)
  .say(sprintf("Wrote %s  (%d marker(s)%s, %s vs %s, n=%d)", basename(rp),
               if (is.null(markers)) 0L else nrow(markers), note, pos, neg, ncol(cs)))
}

# ancombc2 fits lm() on log(count) with R's log(), so lfc is a NATURAL-log fold
# change. Divide by log(2) for the log2 the effect-size axis names.
make_markers <- function(sig, l, q, pos, neg)
  if (any(sig)) data.frame(lineage = keys[sig], class = ifelse(l[sig] > 0, pos, neg),
                           lda = abs(l[sig]) / log(2), pvalue = q[sig], stringsAsFactors = FALSE) else NULL

specs <- list()
if (!pairwise_mode)
  for (lvl in setdiff(levels(meta$Diet), diet_ref))
    specs[[paste0("Diet", lvl)]] <- list(keep = meta$Condition == cond_ref & meta$Diet %in% c(lvl, diet_ref),
      pos = cell_of(lvl, cond_ref), neg = cell_of(diet_ref, cond_ref))
for (lvl in setdiff(levels(meta$Condition), cond_ref))
  specs[[paste0("Condition", lvl)]] <- list(keep = meta$Diet == diet_ref,
    pos = cell_of(diet_ref, lvl), neg = cell_of(diet_ref, cond_ref))

for (coef in names(specs)) {
  sp <- specs[[coef]]
  q <- res_df[[paste0("q_", coef)]][row_of]; l <- res_df[[paste0("lfc_", coef)]][row_of]
  ss <- res_df[[paste0("passed_ss_", coef)]][row_of]
  sig <- !is.na(q) & !is.na(l) & q < ALPHA & (is.na(ss) | ss %in% c(TRUE, "TRUE"))
  write_one(file.path(TAB_DIR, paste0(stem, "_", coef, ".res")), sp$keep,
            make_markers(sig, l, q, sp$pos, sp$neg), sp$pos, sp$neg)
}

# ---- pairwise contrasts: significance is the jointly corrected q (diff_) ----
# pair_df names a contrast "DietB_DietA" for B - A; the .res is named the way a
# reference-level refit would be ("_refA_DietB").
if (pairwise_mode) {
  prow <- match(res_taxa, as.character(pair_df$taxon))
  for (lc in grep("^lfc_", colnames(pair_df), value = TRUE)) {
    k <- sub("^lfc_", "", lc)
    lv <- sub("^Diet", "", strsplit(k, "_")[[1]])
    b <- lv[1]; a <- if (length(lv) > 1) lv[2] else diet_ref
    l <- pair_df[[lc]][prow]; d <- pair_df[[paste0("diff_", k)]][prow]
    ss <- pair_df[[paste0("passed_ss_", k)]][prow]
    sig <- !is.na(l) & d %in% c(TRUE, "TRUE") & (is.na(ss) | ss %in% c(TRUE, "TRUE"))
    pos <- cell_of(b, cond_ref); neg <- cell_of(a, cond_ref)
    tag <- if (a == diet_ref) "" else paste0("_ref", a)
    write_one(file.path(TAB_DIR, paste0(stem, tag, "_Diet", b, ".res")),
              meta$Condition == cond_ref & meta$Diet %in% c(a, b),
              make_markers(sig, l, pair_df[[paste0("q_", k)]][prow], pos, neg), pos, neg, ", joint BH")
  }
}

# Headline counts pass both q < alpha and ancombc2's pseudo-count sensitivity
# check; the raw q count is shown only for contrast.
.say("Significant taxa per coefficient (q < ", ALPHA, " AND passed_ss | raw q < ", ALPHA, "):")
for (qc in setdiff(grep("^q_", colnames(out), value = TRUE), "q_(Intercept)")) {
  coef <- sub("^q_", "", qc)
  sig  <- !is.na(out[[qc]]) & out[[qc]] < ALPHA
  n_ss <- sum(sig & out[[paste0("passed_ss_", coef)]] %in% c(TRUE, "TRUE"))
  message(sprintf("    %-28s %4d | %4d", coef, n_ss, sum(sig)))
}
