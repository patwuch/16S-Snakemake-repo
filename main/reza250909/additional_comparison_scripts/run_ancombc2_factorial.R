#!/usr/bin/env Rscript
# Fit ANCOM-BC2 as a Diet x Radiation factorial for reza250909, with an
# interaction term. Greengenes2 only.
#
# Ported from tisha251224/additional_comparison_scripts/run_ancombc2_factorial.R;
# ahmad260430's copy explains in its header why this lives outside the Snakemake
# pipeline. In short: the pipeline runs ancombc2 with `group = <column>`, and its
# multi-group tests cover main effects only when the formula has an interaction.
# This runs `group = NULL` (per-coefficient inference, BH across taxa).
#
# Needs ANCOMBC >= 2.14.0: run it in containers/ancombc2-env.sif.
#
# Design: 2 diets (Control, SCFAs) x 2 radiation states (None, IR), 8 samples
# per cell, 32 in all. With two levels per factor every question is one model
# coefficient, so the other scripts' --pairwise mode (every pair of three diets)
# has nothing to do here and is not carried over.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_ancombc2_factorial.R [--min-presence N] [--diet-reference SCFAs]
#       [--radiation-reference IR]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
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

# Reference levels pinned explicitly: a relabel must never silently flip the
# sign of a coefficient.
DIET_LEVELS      <- c("Control", "SCFAs")
RADIATION_LEVELS <- c("None", "IR")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
out_suffix <- ""

# Two prevalence rules, both applied:
#  * every Diet x Radiation cell must see the taxon in at least one sample.
#    ANCOM-BC2 treats zeros as missing log values, so a taxon absent from
#    a whole cell leaves the interaction unidentifiable and the fit returns
#    plausible-looking nonsense instead of an error. Not optional, so no flag.
#    With four cells this alone implies presence in >= 4 of 32 samples.
#  * --min-presence N: the taxon must also be seen in >= N samples overall.
#    Runs at a non-default value get their own output file.
DEFAULT_MIN_PRESENCE <- 4L
MIN_SAMPLE_PRESENCE <- as.integer(arg_val("--min-presence", as.character(DEFAULT_MIN_PRESENCE)))
if (is.na(MIN_SAMPLE_PRESENCE)) stop("--min-presence must be an integer")
if (MIN_SAMPLE_PRESENCE != DEFAULT_MIN_PRESENCE) out_suffix <- paste0(out_suffix, "_minpres", MIN_SAMPLE_PRESENCE)

# Every coefficient is conditional on the other factor's reference level, so
# `RadiationIR` is the IR effect within the reference diet only. Refitting with
# the other level as the reference makes ancombc2 estimate it directly, with
# its own standard error.
diet_ref <- arg_val("--diet-reference", DIET_LEVELS[1])
if (!(diet_ref %in% DIET_LEVELS)) stop("--diet-reference must be one of: ", paste(DIET_LEVELS, collapse = ", "))
if (!identical(diet_ref, DIET_LEVELS[1])) {
  DIET_LEVELS <- c(diet_ref, setdiff(DIET_LEVELS, diet_ref))
  out_suffix  <- paste0(out_suffix, "_ref", diet_ref)
}
rad_ref <- arg_val("--radiation-reference", RADIATION_LEVELS[1])
if (!(rad_ref %in% RADIATION_LEVELS)) stop("--radiation-reference must be one of: ", paste(RADIATION_LEVELS, collapse = ", "))
if (!identical(rad_ref, RADIATION_LEVELS[1])) {
  RADIATION_LEVELS <- c(rad_ref, setdiff(RADIATION_LEVELS, rad_ref))
  out_suffix <- paste0(out_suffix, "_radref", rad_ref)
}

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
meta$Radiation <- factor(meta$Radiation, levels = RADIATION_LEVELS)
if (anyNA(meta$Diet) || anyNA(meta$Radiation)) stop("Diet/Radiation outside the declared levels")

.say("Design (samples per cell):")
print(table(meta$Diet, meta$Radiation))

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(DB, "_taxonomy.tsv")))
counts <- collapse_to_level(otu_all[rowSums(otu_all) > 0, , drop = FALSE], tax_map, COLLAPSE_LEVEL)

cell <- interaction(meta$Diet, meta$Radiation, drop = TRUE)
in_every_cell <- sapply(levels(cell), function(k) rowSums(counts[, cell == k, drop = FALSE] > 0) > 0)
n_present <- rowSums(counts > 0)
keep <- rowSums(!in_every_cell) == 0 & n_present >= MIN_SAMPLE_PRESENCE
.say(sprintf("%s: %d features after collapse; %d seen in every cell; %d of those in >= %d samples are tested",
             DB, nrow(counts), sum(rowSums(!in_every_cell) == 0), sum(keep), MIN_SAMPLE_PRESENCE))
.say(sprintf("Tested taxa carry %.1f%% of reads", 100 * sum(counts[keep, ]) / sum(counts)))

# Taxa seen in at least four samples but absent from at least one whole cell
# cannot enter the model. Where they are missing is a result in its own right,
# so record it rather than dropping them silently.
if (out_suffix == "") {
  excl <- n_present >= DEFAULT_MIN_PRESENCE & !keep
  per_cell <- sapply(levels(cell), function(k) rowSums(counts[excl, cell == k, drop = FALSE] > 0))
  colnames(per_cell) <- sapply(levels(cell), function(k) {
    i <- which(cell == k)[1]; sprintf("%s (n=%d)", meta$Group[i], sum(cell == k)) })
  excl_df <- data.frame(taxon = rownames(counts)[excl], n_samples = n_present[excl],
                        per_cell, check.names = FALSE)
  excl_path <- file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", DB, "_untestable_absent_from_a_cell.tsv"))
  write.table(excl_df, excl_path, sep = "\t", quote = FALSE, row.names = FALSE)
  .say("Wrote ", excl_path, " (", nrow(excl_df), " taxa)")
}
counts <- counts[keep, , drop = FALSE]

FIX_FORMULA <- "Diet * Radiation"
.say("Fitting: ", FIX_FORMULA, "  (group = NULL)")
# The empty-cell failure's only symptom is this recycling warning; with the
# per-cell rule above it must not appear, so turn it into an error.
res <- withCallingHandlers(
  ANCOMBC::ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = meta,
           fix_formula = FIX_FORMULA, group = NULL,
           prv_cut = PRV_CUT, lib_cut = LIB_CUT, p_adj_method = P_ADJ_METHOD,
           alpha = ALPHA, pseudo_sens = PSEUDO_SENS, struc_zero = FALSE,
           neg_lb = FALSE, global = FALSE, pairwise = FALSE, dunnet = FALSE,
           trend = FALSE, verbose = FALSE),
  warning = function(w) {
    if (grepl("multiple of", conditionMessage(w)))
      stop("ancombc2 recycled a vector across columns (rank-deficient cell?): ", conditionMessage(w))
  })

res_df <- as.data.frame(res$res)
out <- res_df
out$Taxon <- as.character(out$taxon)
out <- out[, c("taxon", "Taxon", setdiff(colnames(out), c("taxon", "Taxon")))]
stem <- paste0("ANCOMBC2_factorial_", DB, out_suffix)
out_path <- file.path(TAB_DIR, paste0(stem, ".tsv"))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
.say("Wrote ", out_path)

# ---- one .res per main-effect coefficient, for the shared renderer ----
# Each is a two-class contrast conditional on the other factor's reference
# level, so the .res carries exactly that subpopulation. The interaction term
# gets no .res: a difference of differences has no enriched class.
res_taxa <- intersect(as.character(res_df$taxon), rownames(counts))
counts_res <- counts[res_taxa, , drop = FALSE]
lineage_map <- lapply(res_taxa, function(t) canonical_lineage(paste(trimws(strsplit(t, ";")[[1]]), collapse = "|")))
keys <- vapply(lineage_map, paste, character(1), collapse = ".")
row_of <- match(res_taxa, as.character(res_df$taxon))

cell_of <- function(d, r) unique(meta$Group[meta$Diet == d & meta$Radiation == r])
diet_ref <- levels(meta$Diet)[1]; rad_ref <- levels(meta$Radiation)[1]

# ancombc2 fits lm() on log(count) with R's log(), so lfc is a NATURAL-log fold
# change. Divide by log(2) for the log2 the effect-size axis names.
make_markers <- function(sig, l, q, pos, neg)
  if (any(sig)) data.frame(lineage = keys[sig], class = ifelse(l[sig] > 0, pos, neg),
                           lda = abs(l[sig]) / log(2), pvalue = q[sig], stringsAsFactors = FALSE) else NULL

diet_lvl <- levels(meta$Diet)[2]; rad_lvl <- levels(meta$Radiation)[2]
specs <- list()
specs[[paste0("Diet", diet_lvl)]] <- list(keep = meta$Radiation == rad_ref,
  pos = cell_of(diet_lvl, rad_ref), neg = cell_of(diet_ref, rad_ref))
specs[[paste0("Radiation", rad_lvl)]] <- list(keep = meta$Diet == diet_ref,
  pos = cell_of(diet_ref, rad_lvl), neg = cell_of(diet_ref, rad_ref))

for (coef in names(specs)) {
  sp <- specs[[coef]]
  q <- res_df[[paste0("q_", coef)]][row_of]; l <- res_df[[paste0("lfc_", coef)]][row_of]
  ss <- res_df[[paste0("passed_ss_", coef)]][row_of]
  sig <- !is.na(q) & !is.na(l) & q < ALPHA & (is.na(ss) | ss %in% c(TRUE, "TRUE"))
  markers <- make_markers(sig, l, q, sp$pos, sp$neg)
  cs <- counts_res[, rownames(meta)[sp$keep], drop = FALSE]
  cs <- cs[, colSums(cs) > 0, drop = FALSE]
  rp <- file.path(TAB_DIR, paste0(stem, "_", coef, ".res"))
  write_res(cs, as.character(meta[colnames(cs), "Group"]), lineage_map, markers, rp)
  .say(sprintf("Wrote %s  (%d marker(s), %s vs %s, n=%d)", basename(rp),
               if (is.null(markers)) 0L else nrow(markers), sp$pos, sp$neg, ncol(cs)))
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
