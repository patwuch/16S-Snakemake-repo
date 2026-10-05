#!/usr/bin/env Rscript
# Checks LinDA's mode assumption for the PICRUSt2 per-feature tests of one ahmad260430
# report (run_picrust2_linda.R) before any call is read. LinDA assumes most features do
# not change, so a contrast's uncorrected coefficients pile up at one peak whose location
# is the compositional bias, which it subtracts from every feature. This looks only at
# the uncorrected coefficients and that mode, never at p-values or calls.
#
# The fits are run_picrust2_linda.R's (MicrobiomeStat::linda with ggpicrust2's settings:
# count input, adaptive, pseudo-count 0.5, all-zero features dropped), refitted because
# neither ggpicrust2 nor that script keeps the bias:
#   factorial    ~ Group, one fit per reference group, as ggpicrust2's pathway_daa;
#   day17-day25  ~ Timepoint + (1 | MouseID) within each arm, and
#                ~ Arm * Timepoint + (1 | MouseID) for the shift's difference.
#
# Criteria: those of tisha251224's check_picrust2_linda_mode.R (and of its taxa check,
# 2026-09-30), fixed before this first run (2026-10-02). A contrast passes only if all hold:
#   (a) unimodal: Hartigan's dip test p >= 0.05 on the uncorrected coefficients, at both
#       filters (picrust2_mode_dip.R, maaslin3-env, which has diptest);
#   (b) precise: the 95% feature-bootstrap interval of the mode is narrower than 1 log2;
#   (c) no single sample: refitting without any one sample moves the mode by less than
#       0.5 log2;
#   (d) filter-stable: the mode moves by less than 0.5 log2 between all features and the
#       report's main refit (>= 12 of 24 samples, factorial; >= 8 of 16, day17-day25).
#
# Writes <report>/tables/PICRUSt2_linda_mode_check.tsv, _mode_coefs.tsv and _mode_loo.tsv.
#
# Usage, in containers/ggpicrust2-env.sif, then picrust2_mode_dip.R <report> in maaslin3-env.sif:
#   Rscript check_picrust2_linda_mode.R factorial|day17-day25

suppressPackageStartupMessages({ library(MicrobiomeStat); library(modeest); library(readr) })
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
PIC       <- file.path(STUDY_DIR, "picrust2", "pipeline")
WIDTH_MAX <- 1
LOO_MAX   <- 0.5
B         <- 2000L
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
study_meta <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = character(0))

# A fit: its samples, its formula, how to set the meta table up, and contrast -> coefficient.
group_fit <- function(meta, ref, groups) list(
  samples = rownames(meta), formula = "~ Group",
  prep = function(m) { m$Group <- relevel(factor(m$Group), ref = ref); m },
  terms = setNames(as.list(paste0("Group", groups)), paste(groups, "vs", ref)))
REPORTS <- list(
  factorial = list(dir = "comparison-Treatment-x-Radiation-factorial", filter = 12L,
    fits = function() {
      m <- study_meta[study_meta$Timepoint == "Day17", ]
      meta <- data.frame(Group = m$Group, row.names = m$sampleid, stringsAsFactors = FALSE)
      want <- list(Control = c("IFA", "Inosine", "IR"), Inosine = c("IFA", "Inosine+IR"), IR = c("IFA+IR", "Inosine+IR"),
                   "Inosine+IR" = "IFA+IR", IFA = "IFA+IR")
      list(meta = meta, fits = lapply(names(want), function(r) group_fit(meta, r, want[[r]])))
    }),
  "day17-day25" = list(dir = "comparison-day17-day25", filter = 8L,
    fits = function() {
      m <- read.delim(file.path(STUDY_DIR, "comparison-day17-day25", "metadata_joint_cells.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE, na.strings = character(0))
      meta <- data.frame(MouseID = m$MouseID, Arm = factor(sub(" Day.*", "", m$Group), c("IFA", "IFA+IR")),
                         Timepoint = factor(sub(".* ", "", m$Group), c("Day17", "Day25")),
                         row.names = m$sampleid, stringsAsFactors = FALSE)
      arm <- function(a) list(samples = rownames(meta)[meta$Arm == a], formula = "~ Timepoint + (1 | MouseID)",
                              prep = identity, terms = setNames(list("TimepointDay25"), paste(a, "Day25 vs", a, "Day17")))
      list(meta = meta, fits = list(arm("IFA"), arm("IFA+IR"),
        list(samples = rownames(meta), formula = "~ Arm * Timepoint + (1 | MouseID)", prep = identity,
             terms = list("Shift, IFA+IR vs IFA" = "ArmIFA+IR:TimepointDay25"))))
    }))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% names(REPORTS)) stop("first argument: ", paste(names(REPORTS), collapse = " or "))
report <- REPORTS[[args[1]]]
TAB_DIR <- file.path(STUDY_DIR, report$dir, "tables")
spec <- report$fits()
md <- spec$meta

linda_fit <- function(ab, f, cols) {
  x <- ab[, cols, drop = FALSE]
  x <- x[rowSums(x) > 0, , drop = FALSE]
  m <- f$prep(md[cols, , drop = FALSE])
  suppressMessages(suppressWarnings(linda(x, m, f$formula, feature.dat.type = "count", prev.filter = 0,
    mean.abund.filter = 0, adaptive = TRUE, zero.handling = "pseudo-count", pseudo.cnt = 0.5,
    p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE)))
}
raw_of <- function(o, v) o$output[[v]]$log2FoldChange + o$bias[match(v, o$variables)]
mode_of <- function(o, v) o$bias[match(v, o$variables)]
lindamode <- function(x, n) as.numeric(suppressMessages(mlv(sqrt(n) * x, method = "meanshift", kernel = "gaussian"))) / sqrt(n)
boot_mode <- function(x, n) tryCatch(lindamode(x, n), warning = function(w) NA_real_)

set.seed(42)
rows <- list(); coefs <- list(); loos <- list()
for (lv in names(FILES)) {
  ab0 <- as.data.frame(read_delim(file.path(PIC, FILES[[lv]]), delim = "\t", show_col_types = FALSE))
  rownames(ab0) <- ab0[[1]]
  ab0 <- as.matrix(ab0[, rownames(md)])
  for (mp in c(0L, report$filter)) {
    ab <- if (mp) ab0[rowSums(ab0 > 0) >= mp, , drop = FALSE] else ab0
    for (f in spec$fits) {
      o <- linda_fit(ab, f, f$samples)
      n <- length(f$samples)
      for (k in names(f$terms)) {
        v <- f$terms[[k]]
        raw <- raw_of(o, v); m0 <- mode_of(o, v)
        stopifnot(isTRUE(all.equal(m0, lindamode(raw, n))))
        bt <- replicate(B, boot_mode(sample(raw, replace = TRUE), n))
        ci <- quantile(bt, c(0.025, 0.975), na.rm = TRUE)
        loo <- vapply(f$samples, function(s) mode_of(linda_fit(ab, f, setdiff(f$samples, s)), v), 0)
        sh <- loo - m0
        w <- which.max(abs(sh))
        rows[[length(rows) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, features = length(raw),
          mode = m0, median = median(raw), boot_lo = ci[[1]], boot_hi = ci[[2]], boot_width = ci[[2]] - ci[[1]],
          boot_tie_draws = sum(is.na(bt)), loo_max_shift = abs(sh[[w]]), loo_worst_sample = names(sh)[w])
        coefs[[length(coefs) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, coef = raw)
        loos[[length(loos) + 1]] <- data.frame(level = lv, min_presence = mp, contrast = k, sample = names(loo),
                                               mode_without = unname(loo), shift = unname(sh))
        message(sprintf("[mode] %-7s >=%2d %-30s mode %+.2f median %+.2f boot width %.2f loo %.2f (%s)",
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
message("wrote ", file.path(report$dir, "tables", "PICRUSt2_linda_mode_check.tsv"), "; next: picrust2_mode_dip.R ", args[1])
