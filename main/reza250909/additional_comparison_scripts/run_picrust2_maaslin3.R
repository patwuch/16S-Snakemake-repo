#!/usr/bin/env Rscript
# Per-feature tests on reza250909's PICRUSt2 predictions with MaAsLin3 1.5.7, replacing
# ggpicrust2's LinDA, whose mode assumption failed check_picrust2_linda_mode.R in 8 of 12 contrasts
# (2026-10-03), the KO and EC SCFAs+IR contrasts among them; the user chose MaAsLin3 for all four
# reports (2026-10-04). Ported from ahmad260430/additional_comparison_scripts/run_picrust2_maaslin3.R.
#
# Both models are tested against zero (median comparison off): abundance = log share of the
# mouse's predicted total, over mice where the feature is predicted; prevalence = predicted or not.
# Features predicted in >= 2 of 32 mice, filtered once. An abundance coefficient counts only if the
# feature is predicted in every cell its row compares. robust = q < 0.05 with the same sign in every
# leave-one-mouse-out refit.
#
# ~ Diet * Radiation + log10_reads, one fit per pair of reference levels, each of the report's five
# rows a direct coefficient: SCFAs vs Control (unirradiated), SCFAs+IR vs IR, IR vs Control (control
# diet), SCFAs+IR vs SCFAs, and the interaction. BH within each row, as the taxon table does.
# Depth: the mouse's full 16S library (tables/sample_reads.tsv, from the biom).
#
# Writes tables/PICRUSt2_<level>_maaslin3_rows[_noloo].tsv and _rows_summary[_noloo].tsv.
#
# Usage, in containers/maaslin3-env.sif, one core (cores > 1 multiplies memory):
#   Rscript run_picrust2_maaslin3.R MetaCyc|KO|EC [--no-loo]

suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
ALPHA <- 0.05
MIN_SAMPLES <- 2L
LEVELS <- list(
  MetaCyc = c("pipeline/pathways_out/path_abun_unstrat.tsv.gz", "described/pathway_abun_unstrat_described.tsv.gz"),
  KO      = c("pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/KO_metagenome_unstrat_described.tsv.gz"),
  EC      = c("pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz", "described/EC_metagenome_unstrat_described.tsv.gz"))

args <- commandArgs(trailingOnly = TRUE)
REPORT <- "reza250909"; LV <- args[1]
if (is.na(LV) || !LV %in% names(LEVELS)) stop("usage: run_picrust2_maaslin3.R MetaCyc|KO|EC [--no-loo]")
LOO <- !"--no-loo" %in% args
i <- match("--cores", args)
CORES <- if (is.na(i)) 1L else as.integer(args[i + 1L])

OUT <- "comparison-Diet-x-Radiation-factorial"
TAB_DIR <- file.path(STUDY_DIR, OUT, "tables")
PIC <- file.path(STUDY_DIR, OUT, "picrust2")
DIETS <- c("Control", "SCFAs"); RADS <- c("None", "IR")
m <- read.delim(file.path(STUDY_DIR, OUT, "metadata_factorial.tsv"), row.names = 1, stringsAsFactors = FALSE)
rd <- read.delim(file.path(TAB_DIR, "sample_reads.tsv"), row.names = 1)
meta <- data.frame(Diet = m$Diet, Radiation = m$Radiation, MouseID = rownames(m),
                   log10_reads = log10(rd[rownames(m), "reads"]), row.names = rownames(m), stringsAsFactors = FALSE)
stopifnot(!anyNA(meta$log10_reads),
          identical(as.vector(table(factor(meta$Diet, DIETS), factor(meta$Radiation, RADS))), rep(8L, 4)))
row_ <- function(row, fit, term, cells) list(row = row, stratum = "", family = row, fit = fit, term = term,
                                             cells = cells)
cl <- function(d, r) expand.grid(x = d, y = r, stringsAsFactors = FALSE)
ROWS <- list(row_("SCFAs vs Control", "Control None", "DietSCFAs", cl(DIETS, "None")),
             row_("SCFAs+IR vs IR", "Control IR", "DietSCFAs", cl(DIETS, "IR")),
             row_("IR vs Control", "Control None", "RadiationIR", cl("Control", RADS)),
             row_("SCFAs+IR vs SCFAs", "SCFAs None", "RadiationIR", cl("SCFAs", RADS)),
             row_("Interaction: SCFAs vs Control", "Control None", "DietSCFAs:RadiationIR", cl(DIETS, RADS)))
cell_of <- function(rows) paste(meta$Diet[rows], meta$Radiation[rows])
prep <- function(fit_name, rows) {
  ref <- strsplit(fit_name, " ")[[1]]
  m <- meta[rows, , drop = FALSE]
  m$Diet <- relevel(factor(m$Diet, levels = DIETS), ref[1])
  m$Radiation <- relevel(factor(m$Radiation, levels = RADS), ref[2])
  list(m = m, formula = "~ Diet * Radiation + log10_reads")
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
