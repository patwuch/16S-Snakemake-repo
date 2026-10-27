#!/usr/bin/env Rscript
# MaAsLin3 1.5.7 on the Diet x Radiation factorial, Greengenes2 species. Every
# report row is a direct coefficient of a fit whose reference levels make it
# one (maaslin3_common.R). Replaces ANCOM-BC2 and ALDEx3 for the taxa
# (2026-10-22). Model ~ Diet * Radiation + log10_reads, TSS + log2.
#
# Two models per feature:
#   abundance  - log2 share of reads, over the mice where the feature is
#                present. Against zero, or against the row's median feature
#                with --median-comparison (only if check_maaslin3_median.R
#                passes the row).
#   prevalence - present or not, over all mice (logistic, augmented), against
#                zero. A feature missing from a whole group is fitted here.
#                warn_prevalence's flag is kept as a column: a presence effect
#                possibly induced by a stronger abundance one.
#
# An abundance coefficient counts only if the feature is present in at least
# one mouse of every cell the row compares (see identifiable()); otherwise it
# is NA before correction.
#
# BH within each row and model. A call (q < 0.05) is robust when it stays at
# q < 0.05 with the same sign in every one of 32 leave-one-mouse-out refits.
#
# Inputs from prep_maaslin3_input.R. Needs containers/maaslin3-env.sif:
#   Rscript run_maaslin3_rows.R [--min-samples 4] [--median-comparison]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "maaslin3_common.R"))
MEDIAN <- "--median-comparison" %in% args
if (MEDIAN) TAG <- paste0(TAG, "_median")

fit_all <- function(rows = seq_len(n)) setNames(lapply(FITS, fit, rows = rows, median_comparison = MEDIAN), FITS)

rows_table <- function(f, rows = seq_len(n)) {
  out <- do.call(rbind, lapply(seq_along(ROWS), function(k) do.call(rbind, lapply(c("abundance", "prevalence"), function(model) {
    r <- f[[ROWS[[k]]$fit]][[model]]
    x <- r[r$name == ROWS[[k]]$term, ]
    data.frame(ROWTAB[k, c("row", "stratum", "family")], model = model, feature = x$feature,
               coef = x$coef, stderr = x$stderr, p = x$pval_individual, error = x$error,
               N_not_zero = x$N_not_zero, row.names = NULL, stringsAsFactors = FALSE)
  }))))
  bad <- rep(FALSE, nrow(out))
  for (k in seq_along(ROWS)) {
    i <- which(out$row == ROWS[[k]]$row & out$model == "abundance" & !is.na(out$p))
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
message(sprintf("[maaslin3] %s: abundance coefficients set to NA as not identifiable: %d", TAG,
                sum(res$error %in% "absent from a compared cell: coefficient not identifiable")))
res$call <- !is.na(res$q) & res$q < ALPHA
id <- paste(res$row, res$model, res$feature)
hold <- integer(nrow(res))
for (s in seq_len(n)) {
  keep <- setdiff(seq_len(n), s)
  r <- rows_table(fit_all(keep), keep)
  r <- r[match(id, paste(r$row, r$model, r$feature)), ]
  hold <- hold + ifelse(!is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(res$coef), 1L, 0L)
  message(sprintf("[maaslin3] %s refit %d/%d (%.1f min)", TAG, s, n, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
}
res$loo_holds_in <- hold
res$robust <- res$call & hold == n
res$induced_flag <- grepl("^Prevalence association possibly induced", res$error)

key_tab <- read.delim(file.path(MA_DIR, "feature_key.tsv"), stringsAsFactors = FALSE)
res$taxon <- key_tab$taxon[match(res$feature, key_tab$feature)]
stopifnot(!anyNA(res$taxon))
# Mice carrying the feature in each group, for reading presence calls.
for (cell in unique(cell_of(seq_len(n)))) {
  res[[paste0("n_", sub(" ", "_", cell))]] <- unname(colSums(counts[cell_of(seq_len(n)) == cell, , drop = FALSE] > 0)[res$feature])
}
res <- res[order(match(res$row, ROWTAB$row), res$model, res$p), ]
write.table(res, file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_rows.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(ROWTAB$row, function(rw) {
  d <- res[res$row == rw, ]
  a <- d[d$model == "abundance", ]; p <- d[d$model == "prevalence", ]
  data.frame(row = rw,
             abund_fitted = sum(!is.na(a$p)), abund_nominal = sum(a$p < ALPHA, na.rm = TRUE),
             abund_calls = sum(a$call), abund_robust = sum(a$robust),
             abund_up = sum(a$call & a$coef > 0), abund_down = sum(a$call & a$coef < 0),
             median_share_shift = median(a$coef[!is.na(a$p)]),
             prev_fitted = sum(!is.na(p$p)), prev_nominal = sum(p$p < ALPHA, na.rm = TRUE),
             prev_calls = sum(p$call), prev_robust = sum(p$robust),
             prev_induced_flag = sum(p$call & p$induced_flag))
}))
write.table(summ, file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_rows_summary.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(summ, digits = 2, right = FALSE)
message(sprintf("[maaslin3] %s done: %.1f min", TAG, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
