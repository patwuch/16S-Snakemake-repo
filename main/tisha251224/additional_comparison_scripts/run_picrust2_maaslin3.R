#!/usr/bin/env Rscript
# Per-feature tests on PICRUSt2's predictions with MaAsLin3 1.5.7, replacing ggpicrust2's
# LinDA (whose mode assumption failed check_picrust2_linda_mode.R in 23 of 27 contrasts,
# 2026-10-02). Same model, rows and rules as the taxa (run_maaslin3_rows.R):
#   ~ Diet * Condition + log10_reads, TSS + log, one fit per pair of reference levels, each
#   report row read as a direct coefficient (factorial_rows.R), so the nine pairwise rows and
#   the three interaction rows are all tested per feature;
#   abundance  - log share of the sample's predicted total, over samples where the feature is
#                predicted, tested against zero (median comparison off: no assumption that the
#                typical feature is unchanged, so an effect is a change in share);
#   prevalence - predicted or not, over all samples;
#   features predicted in >= 2 of 19 samples (MaAsLin3's min_prevalence 0.10), filtered once;
#   an abundance coefficient counts only if the feature is predicted in every cell its row
#   compares (otherwise R drops an aliased column and the coefficient changes meaning);
#   BH within each (model, family) of the report, not MaAsLin3's pooled q;
#   robust = q < 0.05 with the same sign in all 19 leave-one-sample-out refits.
# Input: PICRUSt2's unstratified predictions (picrust2/pipeline), unrounded.
# Descriptions from PICRUSt2's add_descriptions.py output (picrust2/described).
#
# Writes tables/PICRUSt2_<level>_maaslin3_rows.tsv and _rows_summary.tsv.
#
# Usage, in containers/maaslin3-env.sif (one level at a time keeps memory down):
#   Rscript run_picrust2_maaslin3.R MetaCyc|KO|EC [--no-loo] [--cores N]

suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
source(file.path(THIS_DIR, "factorial_rows.R"))
OUT_DIR <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
TAB_DIR <- file.path(OUT_DIR, "tables")
PIC     <- file.path(OUT_DIR, "picrust2")
ALPHA <- 0.05
MIN_SAMPLES <- 2L
FORMULA <- "~ Diet * Condition + log10_reads"
LEVELS <- list(
  MetaCyc = c("pipeline/pathways_out/path_abun_unstrat.tsv.gz", "described/pathway_abun_unstrat_described.tsv.gz"),
  KO      = c("pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/KO_metagenome_unstrat_described.tsv.gz"),
  EC      = c("pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/EC_metagenome_unstrat_described.tsv.gz"))

args <- commandArgs(trailingOnly = TRUE)
LV <- args[1]
if (is.na(LV) || !LV %in% names(LEVELS)) stop("first argument: MetaCyc, KO or EC")
LOO <- !"--no-loo" %in% args
i <- match("--cores", args)
CORES <- if (is.na(i)) 1L else as.integer(args[i + 1L])

meta <- read.delim(file.path(TAB_DIR, "maaslin3", "metadata.tsv"), row.names = 1, stringsAsFactors = FALSE)
ab <- read.delim(gzfile(file.path(PIC, LEVELS[[LV]][1])), row.names = 1, check.names = FALSE)[, rownames(meta)]
ab <- ab[rowSums(ab > 0) >= MIN_SAMPLES, , drop = FALSE]
# MaAsLin3 wants samples in rows; feature IDs made syntactic and mapped back after.
counts <- as.data.frame(t(ab))
ids <- data.frame(feature = make.names(colnames(counts), unique = TRUE), id = colnames(counts), stringsAsFactors = FALSE)
colnames(counts) <- ids$feature
n <- nrow(counts)
stopifnot(identical(as.vector(table(factor(meta$Diet, DIETS), factor(meta$Condition, CONDS))), c(3L, 3L, 2L, 4L, 4L, 3L)))
message(sprintf("[maaslin3] %s: %d features x %d samples, loo %s, cores %d", LV, ncol(counts), n, LOO, CORES))

row_fits <- unique(ROWS[, c("d_ref", "c_ref")])
key <- function(d, c) paste(d, c)
fit <- function(d_ref, c_ref, rows) {
  m <- meta[rows, , drop = FALSE]
  m$Diet      <- relevel(factor(m$Diet, levels = DIETS), d_ref)
  m$Condition <- relevel(factor(m$Condition, levels = CONDS), c_ref)
  r <- maaslin3(input_data = counts[rows, , drop = FALSE], input_metadata = m, output = tempfile(),
                formula = FORMULA, min_prevalence = 0, normalization = "TSS", transform = "LOG",
                median_comparison_abundance = FALSE, median_comparison_prevalence = FALSE,
                augment = TRUE, warn_prevalence = TRUE, cores = CORES,
                plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE, verbosity = "ERROR")
  list(abundance = r$fit_data_abundance$results, prevalence = r$fit_data_prevalence$results)
}
fit_all <- function(rows = seq_len(n)) {
  f <- lapply(seq_len(nrow(row_fits)), function(k) fit(row_fits$d_ref[k], row_fits$c_ref[k], rows))
  names(f) <- key(row_fits$d_ref, row_fits$c_ref); f
}
identifiable <- function(feature, k, rows) {
  cells <- row_cells(k)
  ok <- rep(TRUE, length(feature))
  for (j in seq_len(nrow(cells))) {
    s <- rows[meta$Diet[rows] == cells$diet[j] & meta$Condition[rows] == cells$cond[j]]
    ok <- ok & colSums(counts[s, feature, drop = FALSE] > 0) > 0
  }
  ok
}
rows_table <- function(f, rows = seq_len(n)) {
  out <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) do.call(rbind, lapply(c("abundance", "prevalence"), function(model) {
    r <- f[[key(ROWS$d_ref[k], ROWS$c_ref[k])]][[model]]
    x <- r[r$name == ROWS$coef[k], ]
    data.frame(ROWS[k, c("row", "stratum", "family")], model = model, feature = x$feature,
               coef = x$coef, stderr = x$stderr, p = x$pval_individual, error = x$error,
               N_not_zero = x$N_not_zero, row.names = NULL, stringsAsFactors = FALSE)
  }))))
  bad <- rep(FALSE, nrow(out))
  for (k in seq_len(nrow(ROWS))) {
    i <- which(out$row == ROWS$row[k] & out$stratum == ROWS$stratum[k] & out$model == "abundance" & !is.na(out$p))
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
message(sprintf("[maaslin3] %s main fit: %.1f min", LV, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
res$call <- !is.na(res$q) & res$q < ALPHA
if (LOO) {
  id <- paste(res$row, res$stratum, res$model, res$feature)
  hold <- integer(nrow(res))
  for (s in seq_len(n)) {
    keep <- setdiff(seq_len(n), s)
    r <- rows_table(fit_all(keep), keep)
    r <- r[match(id, paste(r$row, r$stratum, r$model, r$feature)), ]
    hold <- hold + ifelse(!is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(res$coef), 1L, 0L)
    message(sprintf("[maaslin3] %s refit %d/%d done (%.1f min)", LV, s, n, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  }
  res$loo_holds_in <- hold
  res$robust <- res$call & hold == n
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
res <- res[order(match(paste(res$row, res$stratum), paste(ROWS$row, ROWS$stratum)), res$model, res$p), ]
res <- res[, c("level", "row", "stratum", "family", "model", "id", "description", "n_predicted", "coef", "stderr",
               "p", "q", "call", "loo_holds_in", "robust", "induced_flag", "error", "N_not_zero")]
suffix <- if (LOO) "" else "_noloo"
write.table(res, file.path(TAB_DIR, paste0("PICRUSt2_", LV, "_maaslin3_rows", suffix, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  d <- res[res$row == ROWS$row[k] & res$stratum == ROWS$stratum[k], ]
  a <- d[d$model == "abundance", ]; p <- d[d$model == "prevalence", ]
  data.frame(level = LV, row = ROWS$row[k], stratum = ROWS$stratum[k],
             abund_fitted = sum(!is.na(a$p)), abund_nominal = sum(a$p < ALPHA, na.rm = TRUE),
             abund_calls = sum(a$call), abund_robust = sum(a$robust, na.rm = TRUE),
             abund_up = sum(a$call & a$coef > 0), abund_down = sum(a$call & a$coef < 0),
             median_share_shift = median(a$coef[!is.na(a$p)]),
             prev_fitted = sum(!is.na(p$p)), prev_calls = sum(p$call), prev_robust = sum(p$robust, na.rm = TRUE))
}))
write.table(summ, file.path(TAB_DIR, paste0("PICRUSt2_", LV, "_maaslin3_rows_summary", suffix, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, digits = 2, right = FALSE)
message(sprintf("[maaslin3] %s done: %.1f min", LV, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
