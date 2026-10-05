#!/usr/bin/env Rscript
# Checks LinDA's mode assumption for every report row before any call is read.
# LinDA assumes most taxa do not change, so the uncorrected coefficients of a
# row pile up at one peak whose location is the compositional bias; it
# subtracts that mode from every taxon in the row. This script looks only at
# the uncorrected coefficients and the mode, never at p-values or calls.
#
# Criteria, fixed before the first run (2026-09-30). A row passes only if all
# four hold at both prevalence filters:
#   (a) unimodal: Hartigan's dip test p >= 0.05 on the uncorrected
#       coefficients, under the primary zero handling and every pseudo-count;
#   (b) precise: the 95% taxa-bootstrap interval of the mode is narrower than
#       1 log2 (the mode shifts every fold change in the row by the same
#       amount, so +/- 0.5 log2 is a 1.4-fold error on all of them);
#   (c) no single mouse: dropping any one sample and refitting moves the
#       mode by less than 0.5 log2;
#   (d) filter-stable: the mode moves by less than 0.5 log2 between the >= 6
#       and >= 12 sample filters (checked when both tables exist).
# With 31-35 taxa the dip test has little power, so (a) passing is weak
# evidence; (b)-(d) carry more weight.
#
# Needs containers/linda-env.sif with diptest (envs/linda-env.yaml):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/linda-env.sif \
#     Rscript check_linda_mode.R [--min-presence 6]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "linda_common.R"))
suppressPackageStartupMessages(library(diptest))

DIP_P     <- 0.05
WIDTH_MAX <- 1
LOO_MAX   <- 0.5
FILTER_MAX <- 0.5
B <- 2000
set.seed(42)

n <- ncol(counts)
lindamode <- function(x) suppressMessages(mlv(sqrt(n) * x, method = "meanshift", kernel = "gaussian")) / sqrt(n)
row_fits <- unique(ROWS[, c("d_ref", "c_ref")])
key <- function(d, c) paste(d, c)
raw_of <- function(f, k) f$wald$beta[ROWS$coef[k], ]
mode_of <- function(f, k) f$wald$bias[match(ROWS$coef[k], rownames(f$wald$beta))]
label <- trimws(paste(ROWS$row, ifelse(ROWS$stratum == "", "", paste0("(", ROWS$stratum, ")"))))

fit_all <- function(zero, cols = seq_len(n)) {
  f <- lapply(seq_len(nrow(row_fits)), function(i) fit(row_fits$d_ref[i], row_fits$c_ref[i], zero, cols))
  names(f) <- key(row_fits$d_ref, row_fits$c_ref); f
}
fits <- lapply(settings, fit_all)
names(fits) <- c("primary", paste0("pc", PSEUDO))

# (a) dip test under every zero handling.
dip <- sapply(names(fits), function(s) sapply(seq_len(nrow(ROWS)), function(k)
  dip.test(raw_of(fits[[s]][[key(ROWS$d_ref[k], ROWS$c_ref[k])]], k))$p.value))

# (b) taxa bootstrap of the mode, one resample of taxa shared by all rows.
prim <- fits$primary
raw <- sapply(seq_len(nrow(ROWS)), function(k) raw_of(prim[[key(ROWS$d_ref[k], ROWS$c_ref[k])]], k))
mode0 <- sapply(seq_len(nrow(ROWS)), function(k) mode_of(prim[[key(ROWS$d_ref[k], ROWS$c_ref[k])]], k))
stopifnot(isTRUE(all.equal(mode0, apply(raw, 2, lindamode))))
# A resample repeats taxa, and modeest warns ("encountered a tie ... could be
# multimodal") when the repeats trip its tie rule. Those draws are counted
# and left out of the interval rather than trusted.
boot_mode <- function(x) tryCatch(lindamode(x), warning = function(w) NA_real_)
boot <- t(replicate(B, { i <- sample.int(nrow(raw), replace = TRUE); apply(raw[i, , drop = FALSE], 2, boot_mode) }))
ci <- apply(boot, 2, quantile, c(0.025, 0.975), na.rm = TRUE)

# (c) leave one sample out, refit, recompute each row's mode.
loo <- sapply(seq_len(n), function(s) {
  f <- fit_all(PRIMARY, setdiff(seq_len(n), s))
  sapply(seq_len(nrow(ROWS)), function(k) mode_of(f[[key(ROWS$d_ref[k], ROWS$c_ref[k])]], k))
})
colnames(loo) <- colnames(counts)
loo_shift <- loo - mode0
worst <- apply(abs(loo_shift), 1, which.max)

out <- data.frame(
  row = ROWS$row, stratum = ROWS$stratum, taxa = nrow(counts), zero_handling = PRIMARY,
  mode = mode0, median = apply(raw, 2, median),
  dip_p_min = apply(dip, 1, min), dip_p_primary = dip[, "primary"],
  boot_lo = ci[1, ], boot_hi = ci[2, ], boot_width = ci[2, ] - ci[1, ],
  boot_tie_draws = colSums(is.na(boot)),
  loo_max_shift = abs(loo_shift)[cbind(seq_len(nrow(ROWS)), worst)],
  loo_worst_sample = colnames(counts)[worst],
  stringsAsFactors = FALSE)
out$pass_a_unimodal <- out$dip_p_min >= DIP_P
out$pass_b_precise  <- out$boot_width < WIDTH_MAX
out$pass_c_loo      <- out$loo_max_shift < LOO_MAX
write.table(out, file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, "_mode_check.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(row = rep(label, each = n), sample = colnames(counts),
                       mode_without = as.vector(t(loo)), shift = as.vector(t(loo_shift))),
            file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, "_mode_check_loo.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

pd <- data.frame(panel = factor(label, levels = rev(label)), mode = mode0, median = out$median,
                 lo = ci[1, ], hi = ci[2, ], loo_lo = mode0 + apply(loo_shift, 1, min),
                 loo_hi = mode0 + apply(loo_shift, 1, max))
p <- ggplot(pd, aes(y = panel)) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, orientation = "y", linewidth = 3, colour = "grey80") +
  geom_errorbar(aes(xmin = loo_lo, xmax = loo_hi), width = 0.35, orientation = "y", colour = "#2471A3") +
  geom_point(aes(x = median), shape = 4, size = 2) +
  geom_point(aes(x = mode), colour = "#C0392B", size = 2.2) +
  labs(x = "log2", y = NULL,
       title = sprintf("How well LinDA's mode is pinned, %d taxa", nrow(counts)),
       subtitle = "red: mode; grey bar: 95% taxa bootstrap; blue: range over leave-one-sample-out refits; x: median") +
  theme_bw(base_size = 9)
ggsave(file.path(FIG_DIR, paste0("LinDA_factorial_", TAG, "_mode_check.png")), p, width = 8, height = 5, dpi = 150)

print(out[, c("row", "stratum", "mode", "median", "dip_p_min", "boot_width", "boot_tie_draws", "loo_max_shift",
              "loo_worst_sample", "pass_a_unimodal", "pass_b_precise", "pass_c_loo")],
      digits = 2, right = FALSE)

# (d) and the verdict, once both filters have been checked.
paths <- file.path(TAB_DIR, paste0("LinDA_factorial_", DB, c("", "_minpres12"), "_mode_check.tsv"))
if (all(file.exists(paths))) {
  a <- read.delim(paths[1], stringsAsFactors = FALSE)
  b <- read.delim(paths[2], stringsAsFactors = FALSE)
  stopifnot(identical(a$row, b$row), identical(a$stratum, b$stratum))
  v <- data.frame(row = a$row, stratum = a$stratum,
                  filter_shift = abs(b$mode - a$mode),
                  a_unimodal = a$pass_a_unimodal & b$pass_a_unimodal,
                  b_precise = a$pass_b_precise & b$pass_b_precise,
                  c_loo = a$pass_c_loo & b$pass_c_loo)
  v$d_filter <- v$filter_shift < FILTER_MAX
  v$pass <- v$a_unimodal & v$b_precise & v$c_loo & v$d_filter
  write.table(v, file.path(TAB_DIR, paste0("LinDA_factorial_", DB, "_mode_check_verdict.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nVerdict (both filters):\n")
  print(v, digits = 2, right = FALSE)
}
