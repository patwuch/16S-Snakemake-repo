#!/usr/bin/env Rscript
# MaAsLin3 (1.5.7) for Day17-vs-Day25 in both IFA arms. Settings were fixed in
# comparison-day17-day25/ANALYSIS_LOG.md ("MaAsLin3 as the primary method:
# settings fixed before fitting") before the first fit.
#
# Every fit includes read depth (`reads`), as MaAsLin3's tutorial advises.
# Rows and the fits that estimate each one directly:
#   time shift per arm, interaction : ~ Timepoint * Arm + strata(MouseID), one fit per reference arm
#   averaged main effects           : the same with both factors coded -1/2, +1/2
#   arm difference on one day       : ~ Arm on that day's 8 samples (one per mouse)
# strata(MouseID) is conditional logistic regression in the prevalence model and
# (1|MouseID) in the abundance model; Arm is constant within a mouse, so the
# prevalence model cannot estimate it and the averaged arm row is abundance only.
#
# Input: export_counts_day17_day25.R.
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/maaslin3-env.sif \
#     Rscript run_maaslin3_day17_day25.R [--db Greengenes2|Silva138]

suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
TAB_DIR  <- file.path(dirname(THIS_DIR), "comparison-day17-day25", "tables")

args <- commandArgs(trailingOnly = TRUE)
i <- match("--db", args)
DB <- if (is.na(i)) "Greengenes2" else args[i + 1L]
ALPHA <- 0.05
TOL   <- 0.5   # log2; median stability across taxon sets
# The median comparison simulates the median's null (rnorm, 10,000 draws), so
# p-values depend on the RNG state. Each fit is seeded; SEED_CHECK refits the
# primary set at other seeds.
SEED <- 1
SEED_CHECK <- 2:10

all_counts <- as.matrix(read.delim(file.path(TAB_DIR, sprintf("collapsed_counts_day17day25_%s.tsv", DB)),
                                   row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(TAB_DIR, "samples_day17day25.tsv"), row.names = 1,
                   stringsAsFactors = FALSE)[colnames(all_counts), ]
# Read depth as a covariate, as MaAsLin3's tutorial advises: the sample's total
# in the collapsed table, before any taxon filter.
meta$reads <- colSums(all_counts)

# U: present in >= 2 samples of either arm. S1-S3 as in run_linda_day17_day25.R.
pres <- all_counts > 0
in_arm <- sapply(c("IFA", "IFA_IR"), function(a) rowSums(pres[, meta$Arm == a]) >= 2)
cell <- paste(meta$Arm, meta$Timepoint)
every_cell <- sapply(unique(cell), function(k) rowSums(pres[, cell == k, drop = FALSE]) > 0)
U  <- rowSums(in_arm) > 0
S1 <- U & rowSums(!every_cell) == 0
SETS <- list(U = U, S1 = S1, S2 = S1 & rowSums(pres) >= 8, S3 = rowSums(pres) >= 4)
message(sprintf("[maaslin3] %s taxa: %s", DB, paste(names(SETS), sapply(SETS, sum), collapse = ", ")))

PAIRED <- "~ Timepoint * Arm + reads + strata(MouseID)"
FITS <- list(
  ref_IFA    = list(samples = TRUE, arm = c("IFA", "IFA_IR"), formula = PAIRED),
  ref_IFA_IR = list(samples = TRUE, arm = c("IFA_IR", "IFA"), formula = PAIRED),
  centred    = list(samples = TRUE, centred = TRUE, formula = "~ Time * Radiation + reads + strata(MouseID)"),
  Day17      = list(samples = meta$Timepoint == "Day17", arm = c("IFA", "IFA_IR"), formula = "~ Arm + reads"),
  Day25      = list(samples = meta$Timepoint == "Day25", arm = c("IFA", "IFA_IR"), formula = "~ Arm + reads"))
ROWS <- data.frame(
  row  = c("Day 17 to 25 in IFA", "Day 17 to 25 in IFA+IR", "IFA+IR vs IFA at Day 17",
           "IFA+IR vs IFA at Day 25", "Interaction (IFA+IR shift minus IFA shift)",
           "Average Day 17 to 25 shift", "Average IFA+IR vs IFA"),
  fit  = c("ref_IFA", "ref_IFA_IR", "Day17", "Day25", "ref_IFA", "centred", "centred"),
  term = c("TimepointDay25", "TimepointDay25", "ArmIFA_IR", "ArmIFA_IR",
           "TimepointDay25:ArmIFA_IR", "Time", "Radiation"),
  stringsAsFactors = FALSE)
TIME_ROWS <- ROWS$row[1:2]

run_fit <- function(counts, spec, formula = spec$formula, seed = SEED, ...) {
  set.seed(seed)
  m <- meta[spec$samples, , drop = FALSE]
  if (isTRUE(spec$centred)) {
    m$Time      <- ifelse(m$Timepoint == "Day25", 0.5, -0.5)
    m$Radiation <- ifelse(m$Arm == "IFA_IR", 0.5, -0.5)
  } else {
    m$Arm       <- factor(m$Arm, levels = spec$arm)
    m$Timepoint <- factor(m$Timepoint, levels = c("Day17", "Day25"))
  }
  y <- counts[, rownames(m), drop = FALSE]
  out <- tempfile("maaslin3_")
  o <- suppressWarnings(maaslin3(input_data = t(y), input_metadata = m, output = out, formula = formula,
                                 normalization = "TSS", transform = "LOG", min_prevalence = 0,
                                 plot_summary_plot = FALSE, plot_associations = FALSE,
                                 save_models = FALSE, verbosity = "ERROR", ...))
  unlink(out, recursive = TRUE)
  r <- rbind(o$fit_data_abundance$results, o$fit_data_prevalence$results)
  r$model <- unname(c(linear = "abundance", logistic = "prevalence")[r$model])
  stopifnot(!anyNA(r$model))
  r[, c("feature", "name", "model", "coef", "null_hypothesis", "stderr", "pval_individual", "error", "N_not_zero")]
}

# One row's table from a list of fits: both models, BH within the row over the
# p-values it has. MaAsLin3 skips the median comparison for any feature whose
# fit carries an error string, lme4's boundary-fit warning included, and keeps
# that feature's test against 0. Those p-values are flagged and left out of the
# BH family (decided after the first fit; see ANALYSIS_LOG).
row_table <- function(fits, rows = ROWS, median = TRUE) do.call(rbind, lapply(seq_len(nrow(rows)), function(k) {
  r <- fits[[rows$fit[k]]]
  r <- r[r$name == rows$term[k], ]
  if (nrow(r) == 0) stop("no term ", rows$term[k], " in fit ", rows$fit[k])
  r$skipped_median <- median & r$model == "abundance" & is.na(r$null_hypothesis) & !is.na(r$pval_individual)
  r$q <- NA_real_
  ok <- !is.na(r$pval_individual) & !r$skipped_median
  r$q[ok] <- p.adjust(r$pval_individual[ok], "BH")
  data.frame(row = rows$row[k], r, row.names = NULL, check.names = FALSE)
}))

primary <- list(); med <- list()
for (s in names(SETS)) {
  counts <- all_counts[SETS[[s]], , drop = FALSE]
  fits <- lapply(FITS, run_fit, counts = counts)
  res <- row_table(fits); res$set <- s
  primary[[s]] <- res
  med[[s]] <- do.call(rbind, lapply(ROWS$row, function(r) {
    a <- res[res$row == r & res$model == "abundance", ]
    data.frame(set = s, n_taxa = nrow(counts), row = r,
               median_subtracted = unique(na.omit(a$null_hypothesis))[1],
               abundance_tested = sum(!is.na(a$pval_individual)),
               prevalence_tested = sum(!is.na(res$pval_individual[res$row == r & res$model == "prevalence"])),
               calls_abundance = sum(a$q < ALPHA, na.rm = TRUE),
               calls_prevalence = sum(res$q[res$row == r & res$model == "prevalence"] < ALPHA, na.rm = TRUE))
  }))
}
res <- do.call(rbind, primary)
med <- do.call(rbind, med)

# Sensitivity on U: no median comparison; for the time rows, a random intercept
# in place of strata, and no mouse term at all.
counts <- all_counts[SETS$U, , drop = FALSE]
key <- function(d) paste(d$row, d$feature, d$model)
u <- res[res$set == "U", ]
sens <- list(
  no_median = row_table(lapply(FITS, run_fit, counts = counts, median_comparison_abundance = FALSE), median = FALSE),
  random_intercept = row_table(list(ref_IFA = run_fit(counts, FITS$ref_IFA, "~ Timepoint * Arm + reads + (1|MouseID)", small_random_effects = TRUE),
                                    ref_IFA_IR = run_fit(counts, FITS$ref_IFA_IR, "~ Timepoint * Arm + reads + (1|MouseID)", small_random_effects = TRUE)),
                               ROWS[ROWS$row %in% TIME_ROWS, ]),
  unpaired = row_table(list(ref_IFA = run_fit(counts, FITS$ref_IFA, "~ Timepoint * Arm + reads"),
                            ref_IFA_IR = run_fit(counts, FITS$ref_IFA_IR, "~ Timepoint * Arm + reads")),
                       ROWS[ROWS$row %in% TIME_ROWS, ]))
for (s in names(sens)) {
  d <- sens[[s]]
  u[[paste0("coef_", s)]] <- d$coef[match(key(u), key(d))]
  u[[paste0("q_", s)]]    <- d$q[match(key(u), key(d))]
}

# Seed check on U: calls per row at each other seed.
seed_calls <- do.call(rbind, lapply(SEED_CHECK, function(sd) {
  d <- row_table(lapply(FITS, run_fit, counts = counts, seed = sd))
  data.frame(seed = sd, row = ROWS$row,
             calls_abundance = sapply(ROWS$row, function(r) sum(d$q[d$row == r & d$model == "abundance"] < ALPHA, na.rm = TRUE)),
             calls_prevalence = sapply(ROWS$row, function(r) sum(d$q[d$row == r & d$model == "prevalence"] < ALPHA, na.rm = TRUE)),
             taxa = sapply(ROWS$row, function(r) paste(sort(d$feature[d$row == r & !is.na(d$q) & d$q < ALPHA]), collapse = "|")),
             row.names = NULL)
}))

TAG <- paste0("MaAsLin3_day17day25_", DB)
write.table(seed_calls, file.path(TAB_DIR, paste0(TAG, "_seed_check.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(u, file.path(TAB_DIR, paste0(TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res, file.path(TAB_DIR, paste0(TAG, "_all_sets.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

verdict <- do.call(rbind, lapply(ROWS$row, function(r) {
  d <- med[med$row == r, ]
  x <- setNames(d$median_subtracted, d$set)
  data.frame(row = r, t(x), median_range = diff(range(x)), pass_stable = diff(range(x)) <= TOL, check.names = FALSE)
}))
write.table(med, file.path(TAB_DIR, paste0(TAG, "_median_check.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(verdict, file.path(TAB_DIR, paste0(TAG, "_verdict.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

calls <- function(q) if (all(is.na(q))) NA_integer_ else sum(q < ALPHA, na.rm = TRUE)
summ <- do.call(rbind, lapply(ROWS$row, function(r) {
  d <- u[u$row == r, ]; a <- d$model == "abundance"; p <- !a
  data.frame(row = r, abund = calls(d$q[a]), prev = calls(d$q[p]),
             abund_no_median = calls(d$q_no_median[a]),
             abund_random = calls(d$q_random_intercept[a]), prev_random = calls(d$q_random_intercept[p]),
             abund_unpaired = calls(d$q_unpaired[a]), prev_unpaired = calls(d$q_unpaired[p]))
}))
options(width = 200)
print(med, digits = 3, right = FALSE)
print(verdict, digits = 3, right = FALSE)
print(summ, right = FALSE)
print(aggregate(cbind(calls_abundance, calls_prevalence) ~ row, seed_calls, function(x) paste(range(x), collapse = "-")), right = FALSE)
