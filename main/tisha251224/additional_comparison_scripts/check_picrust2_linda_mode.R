#!/usr/bin/env Rscript
# Checks LinDA's mode assumption for the PICRUSt2 per-feature tests (run_ggpicrust2_linda.R)
# before any call is read. LinDA assumes most features do not change, so a contrast's
# uncorrected coefficients pile up at one peak whose location is the compositional bias,
# which it subtracts from every feature. This looks only at the uncorrected coefficients
# and that mode, never at p-values or calls.
#
# The fits are ggpicrust2's own (pathway_daa -> MicrobiomeStat::linda, formula ~ Group,
# count input, adaptive, pseudo-count 0.5, no filter beyond dropping all-zero features),
# refitted here because ggpicrust2 does not return the bias. Each contrast is the
# coefficient of its second-named group's reference fit, as in run_ggpicrust2_linda.R.
#
# Criteria: the four of check_linda_mode.R (taxa, 2026-09-30), fixed before this first
# run (2026-10-01). A contrast passes only if all hold:
#   (a) unimodal: Hartigan's dip test p >= 0.05 on the uncorrected coefficients, at both
#       filters (picrust2_mode_dip.R, maaslin3-env, which has diptest);
#   (b) precise: the 95% feature-bootstrap interval of the mode is narrower than 1 log2,
#       at both filters;
#   (c) no single sample: refitting without any one sample moves the mode by less than
#       0.5 log2, at both filters;
#   (d) filter-stable: the mode moves by less than 0.5 log2 between all features and the
#       features predicted in >= 10 of 19 samples.
# One difference from the taxa check: ggpicrust2 fixes the zero handling (pseudo-count
# 0.5), so (a) is not repeated across pseudo-counts.
#
# Writes tables/PICRUSt2_linda_mode_check.tsv (one row per level x filter x contrast),
# tables/PICRUSt2_linda_mode_coefs.tsv (uncorrected coefficients, for the dip test) and
# tables/PICRUSt2_linda_mode_loo.tsv.
#
# Usage, in containers/ggpicrust2-env.sif, then picrust2_mode_dip.R in maaslin3-env.sif:
#   Rscript check_picrust2_linda_mode.R

suppressPackageStartupMessages({ library(MicrobiomeStat); library(modeest); library(readr) })
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
TAB_DIR  <- file.path(OUT_DIR, "tables")
PIC      <- file.path(OUT_DIR, "picrust2", "pipeline")

WIDTH_MAX <- 1
LOO_MAX   <- 0.5
FILTER_MAX <- 0.5
FILTER    <- 10L
B         <- 2000L
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
want <- list(Cellulose = c("Inulin", "RS", "Cellulose+IR"), Inulin = c("RS", "Inulin+IR"), RS = "RS+IR",
             `Cellulose+IR` = c("Inulin+IR", "RS+IR"), `Inulin+IR` = "RS+IR")

meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE, na.strings = character(0))
LEVELS <- unique(meta$Group[order(meta$PrimaryOrder)])
md <- data.frame(Group = meta$Group, row.names = meta$sampleid)

fit <- function(ab, ref, cols) {
  x <- ab[, cols, drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]            # ggpicrust2's validate_daa_input(filter_zero = TRUE)
  m <- md[cols, , drop = FALSE]
  m$Group <- relevel(factor(m$Group, levels = LEVELS), ref = ref)
  o <- suppressMessages(suppressWarnings(linda(x, m, "~ Group", feature.dat.type = "count", prev.filter = 0,
         mean.abund.filter = 0, adaptive = TRUE, zero.handling = "pseudo-count", pseudo.cnt = 0.5,
         p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE)))
  o
}
coef_of <- function(o, grp) paste0("Group", grp)
raw_of <- function(o, grp) { v <- coef_of(o, grp); o$output[[v]]$log2FoldChange + o$bias[match(v, o$variables)] }
mode_of <- function(o, grp) o$bias[match(coef_of(o, grp), o$variables)]
lindamode <- function(x, n) as.numeric(suppressMessages(mlv(sqrt(n) * x, method = "meanshift", kernel = "gaussian"))) / sqrt(n)
boot_mode <- function(x, n) tryCatch(lindamode(x, n), warning = function(w) NA_real_)

set.seed(42)
rows <- list(); coefs <- list(); loos <- list()
for (lv in names(FILES)) {
  ab0 <- as.data.frame(read_delim(file.path(PIC, FILES[[lv]]), delim = "\t", show_col_types = FALSE))
  rownames(ab0) <- ab0[[1]]
  ab0 <- as.matrix(ab0[, rownames(md)])
  for (mp in c(0L, FILTER)) {
    ab <- if (mp) ab0[rowSums(ab0 > 0) >= mp, , drop = FALSE] else ab0
    n <- ncol(ab)
    for (ref in names(want)) {
      o <- fit(ab, ref, rownames(md))
      for (grp in want[[ref]]) {
        k <- paste(grp, "vs", ref)
        raw <- raw_of(o, grp); m0 <- mode_of(o, grp)
        stopifnot(isTRUE(all.equal(m0, lindamode(raw, n))))
        bt <- replicate(B, boot_mode(sample(raw, replace = TRUE), n))
        ci <- quantile(bt, c(0.025, 0.975), na.rm = TRUE)
        loo <- vapply(rownames(md), function(s) mode_of(fit(ab, ref, setdiff(rownames(md), s)), grp), 0)
        sh <- loo - m0
        w <- which.max(abs(sh))
        rows[[length(rows) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, features = length(raw),
          mode = m0, median = median(raw), boot_lo = ci[[1]], boot_hi = ci[[2]], boot_width = ci[[2]] - ci[[1]],
          boot_tie_draws = sum(is.na(bt)), loo_max_shift = abs(sh[[w]]), loo_worst_sample = names(sh)[w])
        coefs[[length(coefs) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, coef = raw)
        loos[[length(loos) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, sample = names(loo),
                                               mode_without = unname(loo), shift = unname(sh))
        message(sprintf("[mode] %-7s >=%2d %-26s mode %+.2f median %+.2f boot width %.2f loo %.2f (%s)",
                        lv, mp, k, m0, median(raw), ci[[2]] - ci[[1]], abs(sh[[w]]), names(sh)[w]))
      }
    }
  }
}
out <- do.call(rbind, rows)
out$pass_b_precise <- out$boot_width < WIDTH_MAX
out$pass_c_loo <- out$loo_max_shift < LOO_MAX
write.table(out, file.path(TAB_DIR, "PICRUSt2_linda_mode_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, coefs), file.path(TAB_DIR, "PICRUSt2_linda_mode_coefs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, loos), file.path(TAB_DIR, "PICRUSt2_linda_mode_loo.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
message("wrote tables/PICRUSt2_linda_mode_check.tsv; next: picrust2_mode_dip.R in maaslin3-env")
