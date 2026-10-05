#!/usr/bin/env Rscript
# MaAsLin3 1.5.7 on the Diet x Condition factorial, every report row as a
# direct coefficient of a fit whose reference levels make it one (the diets
# are peers). Model ~ Diet * Condition + log10_reads, TSS + log2.
#
# Two models per feature, both tested against zero:
#   abundance  - log2 share of reads, over the samples where the feature is
#                present. Median comparison is off: check_maaslin3_median.R
#                found the median too weakly pinned to serve as "no change",
#                so an effect here is a change in share of reads, not in
#                absolute amount. The median itself is carried as context.
#   prevalence - presence/absence over all samples (logistic, augmented).
#                MaAsLin3's warn_prevalence flag is kept as a column: a
#                presence effect possibly induced by a stronger abundance one.
#
# An abundance coefficient counts only if the feature is present in at least
# one sample of every cell the row compares. MaAsLin3 fits abundance on present
# samples only; a feature absent from a whole cell makes the model rank
# deficient, R drops one aliased column, and the coefficients left silently
# change meaning (an irradiated diet contrast becomes the unirradiated one)
# with no error. Such rows are set to NA before correction; that taxon's
# absence is a structural zero, reported by the presence tests.
#
# BH within the report's families (per model), not MaAsLin3's pooled q.
# Robustness: every call is refitted without each sample in turn; it holds
# only if q < ALPHA with the same sign in all 19 refits.
#
# Inputs from prep_maaslin3_input.R. Needs containers/maaslin3-env.sif:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/maaslin3-env.sif \
#     Rscript run_maaslin3_rows.R [--min-samples 2]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "maaslin3_common.R"))
ALPHA <- 0.05

fit_all <- function(rows = seq_len(n)) {
  f <- lapply(seq_len(nrow(row_fits)), function(k) fit(row_fits$d_ref[k], row_fits$c_ref[k], rows))
  names(f) <- key(row_fits$d_ref, row_fits$c_ref); f
}

# One long table: every row x model x feature, with q from BH within
# (model, family) over features that fit.
# TRUE where `feature` is present in every cell row k compares, among `rows`.
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
  out[order(match(paste(out$row, out$stratum), paste(ROWS$row, ROWS$stratum)), out$model, out$feature), ]
}

res <- rows_table(fit_all())
message(sprintf("[maaslin3] abundance coefficients set to NA as not identifiable: %d",
                sum(res$error %in% "absent from a compared cell: coefficient not identifiable")))
id <- paste(res$row, res$stratum, res$model, res$feature)

loo_q <- sapply(seq_len(n), function(s) {
  keep <- setdiff(seq_len(n), s)
  r <- rows_table(fit_all(keep), keep)
  r <- r[match(id, paste(r$row, r$stratum, r$model, r$feature)), ]
  ifelse(!is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(res$coef), 1L, 0L)
})
res$call <- !is.na(res$q) & res$q < ALPHA
res$loo_holds_in <- rowSums(loo_q)
res$robust <- res$call & res$loo_holds_in == n
res$induced_flag <- grepl("^Prevalence association possibly induced", res$error)

key_tab <- read.delim(file.path(MA_DIR, "feature_key.tsv"), stringsAsFactors = FALSE)
res$taxon <- key_tab$taxon[match(res$feature, key_tab$feature)]
stopifnot(!anyNA(res$taxon))
write.table(res, file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_rows.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Per row: features fitted, calls, calls holding without any one mouse, and
# the median abundance coefficient over the row's valid features: the typical
# taxon's change in share, context for reading any single taxon's.
summ <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  d <- res[res$row == ROWS$row[k] & res$stratum == ROWS$stratum[k], ]
  a <- d[d$model == "abundance", ]; p <- d[d$model == "prevalence", ]
  data.frame(row = ROWS$row[k], stratum = ROWS$stratum[k],
             abund_fitted = sum(!is.na(a$p)), abund_calls = sum(a$call), abund_robust = sum(a$robust),
             abund_up = sum(a$call & a$coef > 0), abund_down = sum(a$call & a$coef < 0),
             median_share_shift = median(a$coef[!is.na(a$p)]),
             prev_fitted = sum(!is.na(p$p)), prev_calls = sum(p$call), prev_robust = sum(p$robust),
             prev_induced_flag = sum(p$call & p$induced_flag))
}))
write.table(summ, file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_rows_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, digits = 2, right = FALSE)
