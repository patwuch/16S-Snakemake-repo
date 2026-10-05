#!/usr/bin/env Rscript
# Second half of check_picrust2_linda_mode.R: criterion (a), Hartigan's dip test on each
# contrast's uncorrected LinDA coefficients, and the verdict over (a)-(d). Runs in
# containers/maaslin3-env.sif, the image that carries diptest (linda-env.sif and
# ggpicrust2-env.sif do not).
#
# Reads and rewrites tables/PICRUSt2_linda_mode_check.tsv; writes
# tables/PICRUSt2_linda_mode_check_verdict.tsv (one row per level x contrast).
#
# Usage:  Rscript picrust2_mode_dip.R

suppressPackageStartupMessages(library(diptest))
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
TAB_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial", "tables")
DIP_P <- 0.05
FILTER_MAX <- 0.5

chk <- read.delim(file.path(TAB_DIR, "PICRUSt2_linda_mode_check.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
cf  <- read.delim(file.path(TAB_DIR, "PICRUSt2_linda_mode_coefs.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
key <- function(d) paste(d$level, d$min_presence, d$contrast)
dip <- vapply(split(cf$coef, key(cf)), function(x) dip.test(x)$p.value, 0)
chk$dip_p <- unname(dip[key(chk)])
chk$pass_a_unimodal <- chk$dip_p >= DIP_P
write.table(chk, file.path(TAB_DIR, "PICRUSt2_linda_mode_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

a <- chk[chk$min_presence == 0, ]
b <- chk[chk$min_presence != 0, ]
b <- b[match(paste(a$level, a$contrast), paste(b$level, b$contrast)), ]
v <- data.frame(level = a$level, contrast = a$contrast, mode_all = a$mode, mode_filtered = b$mode,
                filter_shift = abs(b$mode - a$mode),
                a_unimodal = a$pass_a_unimodal & b$pass_a_unimodal,
                b_precise = a$pass_b_precise & b$pass_b_precise,
                c_loo = a$pass_c_loo & b$pass_c_loo)
v$d_filter <- v$filter_shift < FILTER_MAX
v$pass <- v$a_unimodal & v$b_precise & v$c_loo & v$d_filter
write.table(v, file.path(TAB_DIR, "PICRUSt2_linda_mode_check_verdict.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(chk[, c("level", "min_presence", "contrast", "features", "mode", "median", "dip_p", "boot_width", "loo_max_shift", "loo_worst_sample")],
      digits = 2, right = FALSE)
cat("\nVerdict:\n"); print(v, digits = 2, right = FALSE)
cat(sprintf("\n%d of %d contrasts pass all four criteria\n", sum(v$pass), nrow(v)))
