#!/usr/bin/env Rscript
# ANCOM-BC2 for Day17-vs-Day25 with both arms in one model:
#   fix_formula  = Timepoint * Arm (Timepoint + Arm + Timepoint:Arm)
#   rand_formula = (1|MouseID)
#
# Why a second model beside run_ancombc2_day17_day25.R. The per-arm fits answer
# "does each arm change" and stay the primary answer to that: each arm's shift is
# estimated from its own four mice, with its own bias correction and residual
# variance. They cannot answer "does the shift differ between arms", which is
# the Timepoint:Arm coefficient here. The joint fit's Timepoint coefficient also
# re-estimates each arm's shift with residual variance pooled across both arms,
# which is the cross-check the report reads against the per-arm counts.
#
# What the terms mean with repeated measures. Arm varies between mice and
# Timepoint within them. The random intercept removes each mouse's baseline from
# the within-mouse terms (Timepoint, Timepoint:Arm); the Arm main effect is a
# between-mouse Day17 comparison the pairing does not help, and the factorial
# analysis already tests it with more mice, so no .res is written for it.
#
# Two parameterisations, one per reference arm, so each arm's shift is a
# coefficient of its own. The interaction is the same contrast with its sign
# flipped between them, which is checked below rather than assumed.
#
# Singular fits. ANCOM-BC2 fits lmerTest::lmer per taxon and never says
# whether the mouse variance came out at zero, in which case that taxon's fit is
# the unpaired one. This refits the same formula on ANCOM-BC2's bias-corrected
# log table and records lme4::isSingular per taxon, so the report can say how
# often the pairing actually did anything.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_ancombc2_day17_day25_joint.R [--db Greengenes2|Silva138]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))

LOG_PREFIX <- "[day17_day25_joint]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
source(file.path(REPO_ROOT, "scripts", "write_res.R"))

suppressPackageStartupMessages({
  library(biomformat)
  library(ANCOMBC)
  library(lme4)
})
require_ancombc("2.14.0")

# Same settings as run_ancombc2_day17_day25.R, so the two models differ only in
# the formula and the samples each sees.
MIN_SAMPLE_PRESENCE <- 2
PRV_CUT             <- 0.05
LIB_CUT             <- 1000
P_ADJ_METHOD        <- "BH"
ALPHA               <- 0.05
PSEUDO_SENS         <- TRUE

ARMS       <- c(IFA = "IFA", IFA_IR = "IFA+IR")   # level -> Group value
DAY_LEVELS <- c("Day17", "Day25")
FIX_FORMULA  <- "Timepoint * Arm"
RAND_FORMULA <- "(1|MouseID)"

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
db <- arg_val("--db", "Greengenes2")

mouse_id <- function(sid) {
  m <- regmatches(sid, regexec("^E\\d+_(\\d+)\\.(\\d+)$", sid))[[1]]
  if (length(m) == 0) m <- regmatches(sid, regexec("^(\\d+)\\.(\\d+)colon$", sid))[[1]]
  if (length(m) == 0) stop("Cannot parse a mouse id out of sample '", sid, "'")
  paste(m[2], m[3], sep = "_")
}


otu_all  <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
meta_all <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), header = TRUE, row.names = 1,
                       sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
meta <- meta_all[meta_all$Group %in% ARMS & rownames(meta_all) %in% colnames(otu_all), , drop = FALSE]
meta$Arm       <- names(ARMS)[match(meta$Group, ARMS)]
meta$Timepoint <- factor(meta$Timepoint, levels = DAY_LEVELS)
meta$MouseID   <- factor(vapply(rownames(meta), mouse_id, character(1)))
if (anyNA(meta$Timepoint)) stop("Timepoint values outside ", paste(DAY_LEVELS, collapse = "/"))

pairing <- table(meta$MouseID, meta$Timepoint)
if (!all(pairing == 1))
  stop("Not cleanly paired; every mouse must appear exactly once per timepoint:\n",
       paste(capture.output(print(pairing)), collapse = "\n"))
if (any(rowSums(table(meta$MouseID, meta$Arm) > 0) != 1))
  stop("A mouse appears in more than one arm")
.say(nrow(meta), " samples, ", nlevels(meta$MouseID), " mice, arms: ",
     paste(names(table(meta$Arm)), table(meta$Arm), collapse = ", "))

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(db, "_taxonomy.tsv")))
level <- COLLAPSE_LEVEL[[db]]
counts <- otu_all[, rownames(meta), drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
counts <- collapse_to_level(counts, tax_map, level)

# The union of the two per-arm feature sets: a feature is tested if it passes
# the per-arm fit's presence rule in either arm. Filtering on all 16 samples
# instead would test a different set, and counts could then differ between the
# models for filter reasons alone.
in_arm <- sapply(names(ARMS), function(a)
  rowSums(counts[, meta$Arm == a, drop = FALSE] > 0) >= MIN_SAMPLE_PRESENCE)
# On top of that, every Arm x Timepoint cell must see the feature at least once.
# ANCOM-BC2 treats a zero as a missing log count, so a feature absent from a
# whole cell leaves the interaction unidentifiable. ANCOMBC then either
# returns no SE at all or, when lme4 drops the rank-deficient column, writes 0
# for the dropped coefficient and recycles per-taxon df across the wrong
# columns (the "number of columns of result is not a multiple" warning), which
# yields confident-looking numbers for a model that was never fitted.
cell <- interaction(meta$Arm, meta$Timepoint, drop = TRUE)
in_every_cell <- sapply(levels(cell), function(k) rowSums(counts[, cell == k, drop = FALSE] > 0) > 0)
keep <- rowSums(in_arm) > 0 & rowSums(!in_every_cell) == 0
.say(sprintf("%s: %d features pass the per-arm presence rule in either arm (%d in IFA, %d in IFA_IR); %d of those are seen in all four cells and are tested",
             db, sum(rowSums(in_arm) > 0), sum(in_arm[, "IFA"]), sum(in_arm[, "IFA_IR"]), sum(keep)))
counts <- counts[keep, , drop = FALSE]

tab_dir <- file.path(STUDY_DIR, "comparison-day17-day25", "tables")
artefact <- function(...) file.path(tab_dir, paste0(paste(c("ANCOMBC2_day17day25_joint", db, ...), collapse = "_")))

singular_flags <- function(bc, md) {
  f <- stats::as.formula(paste("y ~", FIX_FORMULA, "+", RAND_FORMULA))
  ctrl <- lme4::lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  one <- function(t) {
    d <- data.frame(y = as.numeric(bc[t, rownames(md)]), md)
    d <- d[is.finite(d$y), , drop = FALSE]
    fit <- tryCatch(suppressWarnings(suppressMessages(lme4::lmer(f, data = d, control = ctrl))),
                    error = function(e) NULL)
    if (is.null(fit)) return(data.frame(taxon = t, n_obs = nrow(d), singular = NA, mouse_var_share = NA))
    vc <- as.data.frame(lme4::VarCorr(fit))$vcov
    data.frame(taxon = t, n_obs = nrow(d), singular = lme4::isSingular(fit),
               mouse_var_share = vc[1] / sum(vc))
  }
  do.call(rbind, lapply(rownames(bc), one))
}

lineage_of <- function(t) canonical_lineage(paste(trimws(strsplit(t, ";")[[1]]), collapse = "|"))

fits <- list()
for (ref in names(ARMS)) {
  other <- setdiff(names(ARMS), ref)
  md <- meta
  md$Arm <- factor(md$Arm, levels = c(ref, other))
  .say("Fitting ", FIX_FORMULA, " + ", RAND_FORMULA, ", reference arm ", ref)
  res <- ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = md,
                  fix_formula = FIX_FORMULA, rand_formula = RAND_FORMULA,
                  group = NULL, prv_cut = PRV_CUT, lib_cut = LIB_CUT,
                  p_adj_method = P_ADJ_METHOD, alpha = ALPHA,
                  pseudo_sens = PSEUDO_SENS, struc_zero = FALSE, neg_lb = FALSE,
                  global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
                  verbose = FALSE)
  out <- as.data.frame(res$res)
  out$Taxon <- as.character(out$taxon)
  out <- out[, c("taxon", "Taxon", setdiff(colnames(out), c("taxon", "Taxon")))]
  tsv <- paste0(artefact(paste0("armref", ref)), ".tsv")
  write.table(out, tsv, sep = "\t", quote = FALSE, row.names = FALSE)
  .say("Wrote ", tsv)

  bc <- as.matrix(res$bias_correct_log_table)
  sing <- singular_flags(bc, md)
  sing_path <- paste0(artefact(paste0("armref", ref), "singular"), ".tsv")
  write.table(sing, sing_path, sep = "\t", quote = FALSE, row.names = FALSE)
  .say(sprintf("armref %s: %d of %d refits singular (mouse variance estimated at 0)",
               ref, sum(sing$singular, na.rm = TRUE), sum(!is.na(sing$singular))))

  # .res for this reference arm's Day17-vs-Day25 coefficient, drawn over that
  # arm's eight samples so the abundances shown are the comparison made.
  coef <- "TimepointDay25"
  q <- out[[paste0("q_", coef)]]; l <- out[[paste0("lfc_", coef)]]
  ss <- out[[paste0("passed_ss_", coef)]]
  sig <- !is.na(q) & !is.na(l) & q < ALPHA & (is.na(ss) | ss %in% c(TRUE, "TRUE"))
  taxa <- as.character(out$taxon)
  lineage_map <- lapply(taxa, lineage_of)
  keys <- vapply(lineage_map, paste, character(1), collapse = ".")
  markers <- if (any(sig)) data.frame(
    lineage = keys[sig], class = ifelse(l[sig] > 0, "Day25", "Day17"),
    lda = abs(l[sig]) / log(2), pvalue = q[sig], stringsAsFactors = FALSE) else NULL
  arm_samples <- rownames(md)[md$Arm == ref]
  cs <- counts[taxa, arm_samples, drop = FALSE]
  res_path <- paste0(artefact(paste0("armref", ref)), ".res")
  write_res(cs, as.character(md[arm_samples, "Timepoint"]), lineage_map, markers, res_path)

  icoef <- paste0("TimepointDay25:Arm", other)
  iq <- out[[paste0("q_", icoef)]]; iss <- out[[paste0("passed_ss_", icoef)]]
  .say(sprintf("armref %s: shift in %s %d taxa (q < %g and passed_ss); interaction %d taxa; of %d tested",
               ref, ref, sum(sig), ALPHA,
               sum(!is.na(iq) & iq < ALPHA & iss %in% c(TRUE, "TRUE")), nrow(out)))
  fits[[ref]] <- out
}

# The two parameterisations estimate one interaction with opposite signs. If they
# disagree beyond numerical noise, something other than the reference changed.
a <- fits[["IFA"]]; b <- fits[["IFA_IR"]][match(a$taxon, fits[["IFA_IR"]]$taxon), ]
d <- abs(a[["lfc_TimepointDay25:ArmIFA_IR"]] + b[["lfc_TimepointDay25:ArmIFA"]])
.say(sprintf("interaction sign-flip check: max |lfc_a + lfc_b| = %.2e over %d taxa",
             max(d, na.rm = TRUE), sum(!is.na(d))))
.say("Done.")
