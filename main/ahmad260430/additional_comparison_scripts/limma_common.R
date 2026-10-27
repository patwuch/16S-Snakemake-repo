# Shared by run_limma_rows.R and calibrate_shuffle.R: the limma version of
# run_maaslin3_rows.R's abundance model. The caller sources maaslin3_common.R
# first; the fits read its globals (counts, meta, ROWS, FITS), so a caller may
# change meta's labels between calls.
suppressPackageStartupMessages(library(limma))

N_TERMS <- 7L  # intercept, 2 treatment, radiation, 2 interaction, log10_reads
logshare <- t(log2(as.matrix(counts) / rowSums(counts)))  # features x samples
logshare[!is.finite(logshare)] <- NA

# Taxa with an abundance estimate: in every group, with a residual df left.
fittable <- function(rows) {
  ok <- identifiable(colnames(counts), NA, rows)
  ok & colSums(counts[rows, , drop = FALSE] > 0) > N_TERMS
}

limma_rows <- function(rows) {
  keep <- colnames(counts)[fittable(rows)]
  y <- logshare[keep, rows, drop = FALSE]
  fits <- lapply(setNames(FITS, FITS), function(fn) {
    ref <- strsplit(fn, " ")[[1]]
    m <- meta[rows, , drop = FALSE]
    m$Treatment <- relevel(factor(m$Treatment, levels = TREATS), ref[1])
    m$Radiation <- relevel(factor(m$Radiation, levels = RADS), ref[2])
    eBayes(lmFit(y, model.matrix(~ Treatment * Radiation + log10_reads, m)), trend = TRUE)
  })
  s2 <- sapply(fits, function(f) f$sigma)
  stopifnot(max(abs(s2 - s2[, 1])) < 1e-8)
  out <- do.call(rbind, lapply(seq_along(ROWS), function(k) {
    f <- fits[[ROWS[[k]]$fit]]; term <- ROWS[[k]]$term
    se_ord <- f$stdev.unscaled[, term] * f$sigma
    data.frame(ROWTAB[k, c("row", "stratum", "family")], feature = keep,
               coef = unname(f$coefficients[, term]),
               p_ordinary = unname(2 * pt(-abs(f$coefficients[, term] / se_ord), f$df.residual)),
               se = unname(f$stdev.unscaled[, term] * sqrt(f$s2.post)),
               df_total = unname(f$df.total), p = unname(f$p.value[, term]),
               s2_own = unname(f$sigma^2), s2_prior = unname(f$s2.prior), df_prior = f$df.prior[1],
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  out$q <- ave(out$p, out$family, FUN = function(p) p.adjust(p, "BH"))
  out
}

