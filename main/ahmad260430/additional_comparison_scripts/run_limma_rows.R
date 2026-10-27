#!/usr/bin/env Rscript
# limma check of the MaAsLin3 abundance calls (rule fixed in ANALYSIS_LOG.md,
# 2026-10-23, before this ran). Same taxa, same response and same model as
# run_maaslin3_rows.R's abundance model: log2 share of reads (TSS over the
# features maaslin3_common.R keeps), zeros as NA, ~ Treatment * Radiation +
# log10_reads, each report row the direct coefficient of a fit whose reference
# levels make it one. Only the variance differs: eBayes(trend = TRUE) shrinks
# each taxon's residual variance toward an abundance-dependent prior learned
# from all taxa.
#
# One fit per reference pair rather than contrasts.fit: with NAs every taxon
# has its own design rows, and contrasts.fit's standard errors are then only
# approximate. The refits are reparameterisations of one model, so residual
# variances, and the prior, are identical across them (asserted).
#
# Asserts that limma's estimates equal MaAsLin3's coefficients and that its
# unmoderated p-values equal MaAsLin3's, so the two differ in variance alone.
# BH within the MaAsLin3 families; robust = q < 0.05 with the same sign in all
# 24 leave-one-mouse-out refits.
#
# Needs containers/maaslin3-env.sif with limma (maaslin3-env.limma-layer.def):
#   Rscript run_limma_rows.R

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "maaslin3_common.R"))
source(file.path(THIS_DIR, "limma_common.R"))

res <- limma_rows(seq_len(n))
key <- function(d) paste(d$row, d$stratum, d$feature)

# Same estimate and same ordinary test as MaAsLin3's abundance model.
m3 <- read.delim(file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_rows.tsv")), stringsAsFactors = FALSE)
m3 <- m3[m3$model == "abundance" & !is.na(m3$p), ]
m3$stratum[is.na(m3$stratum)] <- ""
stopifnot(setequal(key(res), key(m3)))
i <- match(key(res), key(m3))
d_coef <- max(abs(res$coef - m3$coef[i])); d_p <- max(abs(res$p_ordinary - m3$p[i]))
message(sprintf("[limma] vs MaAsLin3: max |coef diff| %.2g, max |ordinary p diff| %.2g over %d estimates", d_coef, d_p, nrow(res)))
stopifnot(d_coef < 1e-8, d_p < 1e-8)
res$maaslin3_q <- m3$q[i]; res$maaslin3_call <- m3$call[i]; res$maaslin3_robust <- m3$robust[i]
message(sprintf("[limma] prior: df %.1f; prior variance %.2f to %.2f (trend); own variance %.3g to %.3g",
                res$df_prior[1], min(res$s2_prior), max(res$s2_prior), min(res$s2_own), max(res$s2_own)))

res$call <- res$q < ALPHA
hold <- integer(nrow(res))
for (s in seq_len(n)) {
  r <- limma_rows(setdiff(seq_len(n), s))
  r <- r[match(key(res), key(r)), ]
  hold <- hold + ifelse(!is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(res$coef), 1L, 0L)
}
res$loo_holds_in <- hold
res$robust <- res$call & hold == n

key_tab <- read.delim(file.path(MA_DIR, "feature_key.tsv"), stringsAsFactors = FALSE)
res$taxon <- key_tab$taxon[match(res$feature, key_tab$feature)]
res <- res[order(match(paste(res$row, res$stratum), paste(ROWTAB$row, ROWTAB$stratum)), res$p), ]
write.table(res, file.path(MA_DIR, paste0("limma_factorial_", TAG, "_rows.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(seq_len(nrow(ROWTAB)), function(k) {
  d <- res[res$row == ROWTAB$row[k] & res$stratum == ROWTAB$stratum[k], ]
  data.frame(row = ROWTAB$row[k], stratum = ROWTAB$stratum[k], tested = nrow(d),
             limma_calls = sum(d$call), limma_robust = sum(d$robust),
             maaslin3_calls = sum(d$maaslin3_call), maaslin3_robust = sum(d$maaslin3_robust),
             both_robust = sum(d$robust & d$maaslin3_robust))
}))
write.table(summ, file.path(MA_DIR, paste0("limma_factorial_", TAG, "_rows_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
options(width = 200)
print(summ, right = FALSE)
