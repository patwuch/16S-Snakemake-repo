#!/usr/bin/env Rscript
# Checks MaAsLin3's median assumption for every report row before any call is
# read. Ported from reza250909 (itself from tisha251224). MaAsLin3's abundance model (log2 TSS, non-zero
# samples only) can test each feature's coefficient against the median
# coefficient of that term, taken over features that fit without error and
# have p < 0.95; this assumes the median feature is unchanged. This script
# looks only at the coefficients and that median, never at p-values against
# it or at calls.
#
# Same four criteria and thresholds as tisha251224, fixed before the first
# run. A row passes only if all four hold at both filters:
#   (a) unimodal: Hartigan's dip test p >= 0.05 on the coefficients the median
#       is taken over;
#   (b) precise: the 95% feature-bootstrap interval of the median is narrower
#       than 1 log2;
#   (c) no single mouse: dropping any one sample and refitting moves the
#       median by less than 0.5 log2;
#   (d) filter-stable: the median moves by less than 0.5 log2 between the two
#       filters (checked when both tables exist).
# Filters: a feature non-zero in >= 3 of 24 samples (MaAsLin3's min_prevalence
# 0.10) and >= 6 (the dot map's and ALDEx3's filter).
#
# Inputs from prep_maaslin3_input.R. Needs containers/maaslin3-env.sif (1.5.7):
#   Rscript check_maaslin3_median.R [--min-samples 3]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "maaslin3_common.R"))
suppressPackageStartupMessages({ library(diptest); library(ggplot2) })

DIP_P      <- 0.05
WIDTH_MAX  <- 1
LOO_MAX    <- 0.5
FILTER_MAX <- 0.5
FILTERS    <- c(3L, 6L)
B <- 2000

# The coefficients MaAsLin3 takes the median over: fitted without error, p < 0.95.
eligible <- function(res, term) {
  x <- res[res$name == term, ]
  x$coef[is.na(x$error) & !is.na(x$pval_individual) & x$pval_individual < 0.95]
}
fit_all <- function(rows = seq_len(n)) setNames(lapply(FITS, function(f) fit(f, rows, abundance_only = TRUE)$abundance), FITS)
coefs_of <- function(f) lapply(ROWS, function(r) eligible(f[[r$fit]], r$term))

full <- fit_all()
co <- coefs_of(full)
med0 <- sapply(co, median)

# Guard: the median recomputed here is the one MaAsLin3 tests against.
for (fn in FITS) {
  r <- fit(fn, median_comparison = TRUE, abundance_only = TRUE)$abundance
  for (j in which(ROWTAB$fit == fn)) {
    nh <- unique(r$null_hypothesis[r$name == ROWTAB$term[j] & !is.na(r$null_hypothesis)])
    stopifnot(length(nh) == 1, isTRUE(all.equal(nh, med0[j])))
  }
}

set.seed(42)
dip_p <- sapply(co, function(x) dip.test(x)$p.value)
boot <- sapply(co, function(x) replicate(B, median(sample(x, replace = TRUE))))
ci <- apply(boot, 2, quantile, c(0.025, 0.975))

loo <- sapply(seq_len(n), function(s) sapply(coefs_of(fit_all(setdiff(seq_len(n), s))), median))
colnames(loo) <- rownames(counts)
loo_shift <- loo - med0
worst <- apply(abs(loo_shift), 1, which.max)

LABEL <- trimws(paste(ROWTAB$row, ifelse(ROWTAB$stratum == "", "", paste0("(", ROWTAB$stratum, ")"))))
out <- data.frame(
  row = ROWTAB$row, stratum = ROWTAB$stratum, features = ncol(counts), median_over = lengths(co),
  median = med0, dip_p = dip_p,
  boot_lo = ci[1, ], boot_hi = ci[2, ], boot_width = ci[2, ] - ci[1, ],
  loo_max_shift = abs(loo_shift)[cbind(seq_along(ROWS), worst)],
  loo_worst_sample = colnames(loo)[worst], stringsAsFactors = FALSE)
out$pass_a_unimodal <- out$dip_p >= DIP_P
out$pass_b_precise  <- out$boot_width < WIDTH_MAX
out$pass_c_loo      <- out$loo_max_shift < LOO_MAX
write.table(out, file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_median_check.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(row = rep(LABEL, each = n), sample = colnames(loo),
                       median_without = as.vector(t(loo)), shift = as.vector(t(loo_shift))),
            file.path(MA_DIR, paste0("MaAsLin3_factorial_", TAG, "_median_check_loo.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

pd <- data.frame(panel = factor(LABEL, levels = rev(LABEL)), median = med0, lo = ci[1, ], hi = ci[2, ],
                 loo_lo = med0 + apply(loo_shift, 1, min), loo_hi = med0 + apply(loo_shift, 1, max))
p <- ggplot(pd, aes(y = panel)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 3, colour = "grey80") +
  geom_errorbar(aes(xmin = loo_lo, xmax = loo_hi), width = 0.35, orientation = "y", colour = "#2471A3") +
  geom_point(aes(x = median), colour = "#C0392B", size = 2.2) +
  labs(x = "log2", y = NULL,
       title = sprintf("How well MaAsLin3's median is pinned, features in >= %d samples", MIN_SAMPLES),
       subtitle = "red: median tested against; grey bar: 95% feature bootstrap\nblue: range over leave-one-sample-out refits") +
  theme_bw(base_size = 9)
ggsave(file.path(FIG_DIR, paste0("MaAsLin3_factorial_", TAG, "_median_check.png")), p, width = 8, height = 5, dpi = 150)

options(width = 200)
print(out[, c("row", "stratum", "median_over", "median", "dip_p", "boot_width", "loo_max_shift",
              "loo_worst_sample", "pass_a_unimodal", "pass_b_precise", "pass_c_loo")],
      digits = 2, right = FALSE)

paths <- file.path(MA_DIR, paste0("MaAsLin3_factorial_Greengenes2_Day17_min", FILTERS, "_median_check.tsv"))
if (all(file.exists(paths))) {
  a <- read.delim(paths[1], stringsAsFactors = FALSE)
  b <- read.delim(paths[2], stringsAsFactors = FALSE)
  stopifnot(identical(a$row, b$row), identical(a$stratum, b$stratum))
  v <- data.frame(row = a$row, stratum = a$stratum, filter_shift = abs(b$median - a$median),
                  a_unimodal = a$pass_a_unimodal & b$pass_a_unimodal,
                  b_precise = a$pass_b_precise & b$pass_b_precise,
                  c_loo = a$pass_c_loo & b$pass_c_loo)
  v$d_filter <- v$filter_shift < FILTER_MAX
  v$pass <- v$a_unimodal & v$b_precise & v$c_loo & v$d_filter
  write.table(v, file.path(MA_DIR, "MaAsLin3_factorial_Greengenes2_Day17_median_check_verdict.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nVerdict (both filters):\n")
  print(v, digits = 2, right = FALSE)
}
