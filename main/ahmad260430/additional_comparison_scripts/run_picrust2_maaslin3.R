#!/usr/bin/env Rscript
# Per-feature tests on PICRUSt2's predictions (picrust2/) with MaAsLin3 1.5.7, for both
# ahmad260430 reports, replacing LinDA, whose mode assumption failed
# check_picrust2_linda_mode.R (factorial 24 of 27 contrasts, day17-day25 9 of 9; 2026-10-02).
# Ported from tisha251224/additional_comparison_scripts/run_picrust2_maaslin3.R.
#
# Both models are tested against zero (median comparison off): abundance = log share of the
# sample's predicted total, over samples where the feature is predicted; prevalence = predicted
# or not. Features predicted in >= 2 of the report's samples, filtered once. An abundance
# coefficient counts only if the feature is predicted in every cell its row compares.
# robust = q < 0.05 with the same sign in every leave-one-mouse-out refit (a mouse's samples
# leave together).
#
# factorial   (24 Day 17 mice, 4 per group): ~ Treatment * Radiation + log10_reads, one fit per
#             pair of reference levels; every row a direct coefficient: the three treatment pairs
#             within each radiation state, radiation within each treatment, the three
#             interaction contrasts. BH within each family, as tisha251224: the treatment pairs of
#             one radiation state together, each radiation row alone, the interaction rows together.
#             Depth: the mouse's full 16S library (tables/sample_reads_day17.tsv, from the biom).
# day17-day25 (IFA and IFA+IR arms, 4 mice each, both days): the rows and fits of
#             run_maaslin3_day17_day25.R: ~ Timepoint * Arm + reads + strata(MouseID) per reference
#             arm, centred +/-1/2 coding for the averaged effects, ~ Arm + reads on each day's 8
#             samples. BH within each row, as that script does. Under strata the prevalence model
#             cannot estimate Arm (constant within a mouse). Depth: as that script, the sample's
#             total in collapsed_counts_day17day25_Greengenes2.tsv.
#
# Writes <report>/tables/PICRUSt2_<level>_maaslin3_rows[_noloo].tsv and _rows_summary[_noloo].tsv.
#
# Usage, in containers/maaslin3-env.sif (one level at a time keeps memory down):
#   Rscript run_picrust2_maaslin3.R factorial|day17-day25 MetaCyc|KO|EC [--no-loo] [--cores N]

suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
PIC       <- file.path(STUDY_DIR, "picrust2")
ALPHA <- 0.05
MIN_SAMPLES <- 2L
LEVELS <- list(
  MetaCyc = c("pipeline/pathways_out/path_abun_unstrat.tsv.gz", "described/pathway_abun_unstrat_described.tsv.gz"),
  KO      = c("pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/KO_metagenome_unstrat_described.tsv.gz"),
  EC      = c("pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/EC_metagenome_unstrat_described.tsv.gz"))

args <- commandArgs(trailingOnly = TRUE)
REPORT <- args[1]; LV <- args[2]
if (is.na(REPORT) || is.na(LV) || !REPORT %in% c("factorial", "day17-day25") || !LV %in% names(LEVELS))
  stop("usage: run_picrust2_maaslin3.R factorial|day17-day25 MetaCyc|KO|EC [--no-loo] [--cores N]")
LOO <- !"--no-loo" %in% args
i <- match("--cores", args)
CORES <- if (is.na(i)) 1L else as.integer(args[i + 1L])

if (REPORT == "factorial") {
  OUT <- "comparison-Treatment-x-Radiation-factorial"
  TAB_DIR <- file.path(STUDY_DIR, OUT, "tables")
  TREATS <- c("Control", "IFA", "Inosine"); RADS <- c("None", "IR")
  m <- read.delim(file.path(TAB_DIR, "metadata_day17.tsv"), row.names = 1, stringsAsFactors = FALSE)
  rd <- read.delim(file.path(TAB_DIR, "sample_reads_day17.tsv"), row.names = 1)
  meta <- data.frame(Treatment = m$Treatment, Radiation = m$Radiation, MouseID = rownames(m),
                     log10_reads = log10(rd[rownames(m), "reads"]), row.names = rownames(m), stringsAsFactors = FALSE)
  stopifnot(!anyNA(meta$log10_reads),
            identical(as.vector(table(factor(meta$Treatment, TREATS), factor(meta$Radiation, RADS))), rep(4L, 6)))
  rad_name <- c(None = "unirradiated", IR = "irradiated")
  pair <- function(a, b, rad) list(row = paste(a, "vs", b), stratum = rad, family = paste0("treatment_", rad),
                                   fit = paste(b, rad), term = paste0("Treatment", a),
                                   cells = data.frame(x = c(a, b), y = rad, stringsAsFactors = FALSE))
  within <- function(t) list(row = paste("IR within", t), stratum = "", family = paste0("IR_", t), fit = paste(t, "None"),
                             term = "RadiationIR", cells = data.frame(x = t, y = RADS, stringsAsFactors = FALSE))
  inter <- function(a, b) list(row = paste("Interaction:", a, "vs", b), stratum = "", family = "interaction",
                               fit = paste(b, "None"), term = paste0("Treatment", a, ":RadiationIR"),
                               cells = expand.grid(x = c(a, b), y = RADS, stringsAsFactors = FALSE))
  ROWS <- list(pair("IFA", "Control", "None"), pair("Inosine", "Control", "None"), pair("IFA", "Inosine", "None"),
               pair("IFA", "Control", "IR"), pair("Inosine", "Control", "IR"), pair("IFA", "Inosine", "IR"),
               within("Control"), within("IFA"), within("Inosine"),
               inter("IFA", "Control"), inter("Inosine", "Control"), inter("IFA", "Inosine"))
  cell_of <- function(rows) paste(meta$Treatment[rows], meta$Radiation[rows])
  prep <- function(fit_name, rows) {
    ref <- strsplit(fit_name, " ")[[1]]
    m <- meta[rows, , drop = FALSE]
    m$Treatment <- relevel(factor(m$Treatment, levels = TREATS), ref[1])
    m$Radiation <- relevel(factor(m$Radiation, levels = RADS), ref[2])
    list(m = m, formula = "~ Treatment * Radiation + log10_reads")
  }
} else {
  OUT <- "comparison-day17-day25"
  TAB_DIR <- file.path(STUDY_DIR, OUT, "tables")
  s <- read.delim(file.path(TAB_DIR, "samples_day17day25.tsv"), row.names = 1, stringsAsFactors = FALSE)
  tot <- read.delim(file.path(TAB_DIR, "collapsed_counts_day17day25_Greengenes2.tsv"), row.names = 1, check.names = FALSE)
  meta <- data.frame(Arm = s$Arm, Timepoint = s$Timepoint, MouseID = s$MouseID,
                     reads = colSums(tot)[rownames(s)], row.names = rownames(s), stringsAsFactors = FALSE)
  stopifnot(!anyNA(meta$reads))
  ARMS <- c("IFA", "IFA_IR"); DAYS <- c("Day17", "Day25")
  cells2 <- function(a, d) expand.grid(x = a, y = d, stringsAsFactors = FALSE)
  row_ <- function(row, fit, term, cells) list(row = row, stratum = "", family = row, fit = fit, term = term, cells = cells)
  ROWS <- list(row_("Day 17 to 25 in IFA", "ref_IFA", "TimepointDay25", cells2("IFA", DAYS)),
               row_("Day 17 to 25 in IFA+IR", "ref_IFA_IR", "TimepointDay25", cells2("IFA_IR", DAYS)),
               row_("IFA+IR vs IFA at Day 17", "Day17", "ArmIFA_IR", cells2(ARMS, "Day17")),
               row_("IFA+IR vs IFA at Day 25", "Day25", "ArmIFA_IR", cells2(ARMS, "Day25")),
               row_("Interaction (IFA+IR shift minus IFA shift)", "ref_IFA", "TimepointDay25:ArmIFA_IR", cells2(ARMS, DAYS)),
               row_("Average Day 17 to 25 shift", "centred", "Time", cells2(ARMS, DAYS)),
               row_("Average IFA+IR vs IFA", "centred", "Radiation", cells2(ARMS, DAYS)))
  cell_of <- function(rows) paste(meta$Arm[rows], meta$Timepoint[rows])
  prep <- function(fit_name, rows) {
    m <- meta[rows, , drop = FALSE]
    if (fit_name == "centred") {
      m$Time <- ifelse(m$Timepoint == "Day25", 0.5, -0.5); m$Radiation <- ifelse(m$Arm == "IFA_IR", 0.5, -0.5)
      return(list(m = m, formula = "~ Time * Radiation + reads + strata(MouseID)"))
    }
    if (fit_name %in% DAYS) {
      m <- m[m$Timepoint == fit_name, , drop = FALSE]
      m$Arm <- factor(m$Arm, levels = ARMS)
      return(list(m = m, formula = "~ Arm + reads"))
    }
    m$Arm <- factor(m$Arm, levels = if (fit_name == "ref_IFA") ARMS else rev(ARMS))
    m$Timepoint <- factor(m$Timepoint, levels = DAYS)
    list(m = m, formula = "~ Timepoint * Arm + reads + strata(MouseID)")
  }
}
FITS <- unique(vapply(ROWS, `[[`, "", "fit"))
ROWTAB <- do.call(rbind, lapply(ROWS, function(r) data.frame(row = r$row, stratum = r$stratum, family = r$family,
                                                             fit = r$fit, term = r$term, stringsAsFactors = FALSE)))

ab <- read.delim(gzfile(file.path(PIC, LEVELS[[LV]][1])), row.names = 1, check.names = FALSE)[, rownames(meta)]
ab <- ab[rowSums(ab > 0) >= MIN_SAMPLES, , drop = FALSE]
counts <- as.data.frame(t(ab))
ids <- data.frame(feature = make.names(colnames(counts), unique = TRUE), id = colnames(counts), stringsAsFactors = FALSE)
colnames(counts) <- ids$feature
n <- nrow(counts)
message(sprintf("[maaslin3] %s %s: %d features x %d samples, loo %s, cores %d", REPORT, LV, ncol(counts), n, LOO, CORES))

fit <- function(fit_name, rows) {
  set.seed(1)
  p <- prep(fit_name, rows)
  r <- maaslin3(input_data = counts[rownames(p$m), , drop = FALSE], input_metadata = p$m, output = tempfile(),
                formula = p$formula, min_prevalence = 0, normalization = "TSS", transform = "LOG",
                median_comparison_abundance = FALSE, median_comparison_prevalence = FALSE,
                augment = TRUE, warn_prevalence = TRUE, cores = CORES,
                plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE, verbosity = "ERROR")
  list(abundance = r$fit_data_abundance$results, prevalence = r$fit_data_prevalence$results)
}
fit_all <- function(rows = seq_len(n)) setNames(lapply(FITS, fit, rows = rows), FITS)
identifiable <- function(feature, k, rows) {
  cl <- ROWS[[k]]$cells
  ok <- rep(TRUE, length(feature))
  for (j in seq_len(nrow(cl))) {
    s <- rows[cell_of(rows) == paste(cl$x[j], cl$y[j])]
    ok <- ok & colSums(counts[s, feature, drop = FALSE] > 0) > 0
  }
  ok
}
rows_table <- function(f, rows = seq_len(n)) {
  out <- do.call(rbind, lapply(seq_along(ROWS), function(k) do.call(rbind, lapply(c("abundance", "prevalence"), function(model) {
    r <- f[[ROWS[[k]]$fit]][[model]]
    x <- r[r$name == ROWS[[k]]$term, ]
    if (!nrow(x)) return(NULL)   # e.g. Arm in the prevalence model under strata
    data.frame(ROWTAB[k, c("row", "stratum", "family")], model = model, feature = x$feature,
               coef = x$coef, stderr = x$stderr, p = x$pval_individual, error = x$error,
               N_not_zero = x$N_not_zero, row.names = NULL, stringsAsFactors = FALSE)
  }))))
  bad <- rep(FALSE, nrow(out))
  for (k in seq_along(ROWS)) {
    i <- which(out$row == ROWS[[k]]$row & out$stratum == ROWS[[k]]$stratum & out$model == "abundance" & !is.na(out$p))
    bad[i] <- !identifiable(out$feature[i], k, rows)
  }
  out$error[bad] <- "absent from a compared cell: coefficient not identifiable"
  out$p[bad] <- NA
  fitted <- !is.na(out$p) & (is.na(out$error) | grepl("^Prevalence association possibly induced", out$error))
  out$p[!fitted] <- NA
  out$q <- ave(out$p, out$model, out$family, FUN = function(p) p.adjust(p, "BH"))
  out
}

t0 <- Sys.time()
res <- rows_table(fit_all())
message(sprintf("[maaslin3] %s %s main fit: %.1f min", REPORT, LV, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res$call <- !is.na(res$q) & res$q < ALPHA
mice <- unique(meta$MouseID)
if (LOO) {
  id <- paste(res$row, res$stratum, res$model, res$feature)
  hold <- integer(nrow(res))
  for (k in seq_along(mice)) {
    keep <- which(meta$MouseID != mice[k])
    r <- rows_table(fit_all(keep), keep)
    r <- r[match(id, paste(r$row, r$stratum, r$model, r$feature)), ]
    hold <- hold + ifelse(!is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(res$coef), 1L, 0L)
    message(sprintf("[maaslin3] %s %s refit %d/%d done (%.1f min)", REPORT, LV, k, length(mice),
                    as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res$loo_holds_in <- hold
  res$robust <- res$call & hold == length(mice)
} else {
  res$loo_holds_in <- NA_integer_
  res$robust <- NA
}
res$induced_flag <- grepl("^Prevalence association possibly induced", res$error)

desc <- read.delim(gzfile(file.path(PIC, LEVELS[[LV]][2])), check.names = FALSE, quote = "", stringsAsFactors = FALSE)
res$id <- ids$id[match(res$feature, ids$feature)]
stopifnot(!anyNA(res$id))
res$description <- desc$description[match(res$id, desc[[1]])]
res$n_predicted <- unname(colSums(counts > 0)[res$feature])
res$level <- LV
res <- res[order(match(paste(res$row, res$stratum), paste(ROWTAB$row, ROWTAB$stratum)), res$model, res$p), ]
res <- res[, c("level", "row", "stratum", "family", "model", "id", "description", "n_predicted", "coef", "stderr",
               "p", "q", "call", "loo_holds_in", "robust", "induced_flag", "error", "N_not_zero")]
suffix <- if (LOO) "" else "_noloo"
write.table(res, file.path(TAB_DIR, paste0("PICRUSt2_", LV, "_maaslin3_rows", suffix, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
summ <- do.call(rbind, lapply(seq_len(nrow(ROWTAB)), function(k) {
  d <- res[res$row == ROWTAB$row[k] & res$stratum == ROWTAB$stratum[k], ]
  a <- d[d$model == "abundance", ]; p <- d[d$model == "prevalence", ]
  data.frame(level = LV, row = ROWTAB$row[k], stratum = ROWTAB$stratum[k],
             abund_fitted = sum(!is.na(a$p)), abund_nominal = sum(a$p < ALPHA, na.rm = TRUE),
             abund_calls = sum(a$call), abund_robust = sum(a$robust, na.rm = TRUE),
             abund_up = sum(a$call & a$coef > 0), abund_down = sum(a$call & a$coef < 0),
             median_share_shift = if (any(!is.na(a$p))) median(a$coef[!is.na(a$p)]) else NA,
             prev_fitted = sum(!is.na(p$p)), prev_calls = sum(p$call), prev_robust = sum(p$robust, na.rm = TRUE))
}))
write.table(summ, file.path(TAB_DIR, paste0("PICRUSt2_", LV, "_maaslin3_rows_summary", suffix, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(summ, digits = 2, right = FALSE)
message(sprintf("[maaslin3] %s %s done: %.1f min", REPORT, LV, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
