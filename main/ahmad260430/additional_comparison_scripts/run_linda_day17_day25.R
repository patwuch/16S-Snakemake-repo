#!/usr/bin/env Rscript
# LinDA pilot (zhouhj1994/LinDA 0.2.0) for Day17-vs-Day25 in both IFA arms:
#   ~ Timepoint * Arm + (1|MouseID)
# Settings and pass criteria were fixed in comparison-day17-day25/ANALYSIS_LOG.md
# ("LinDA pilot: settings and pass criteria fixed before fitting") before the
# first fit.
#
# Every report row is a direct coefficient of a fit whose reference levels make
# it one (LinDA subtracts a separate mode from each coefficient). The averaged
# main effects come from a fit with both factors coded -1/2, +1/2.
#
# Input: export_counts_day17_day25.R (collapsed table, 16 samples).
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/linda-env.sif \
#     Rscript run_linda_day17_day25.R [--db Greengenes2|Silva138]

suppressPackageStartupMessages({ library(LinDA); library(modeest); library(ggplot2) })
stopifnot(packageVersion("LinDA") == "0.2.0")

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR   <- file.path(dirname(THIS_DIR), "comparison-day17-day25")
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")

args <- commandArgs(trailingOnly = TRUE)
i <- match("--db", args)
DB <- if (is.na(i)) "Greengenes2" else args[i + 1L]
ALPHA  <- 0.05
PSEUDO <- c(0.1, 0.5, 1)
TOL    <- 0.5   # log2; both pass criteria
FORMULA <- "~ Timepoint * Arm + (1|MouseID)"

all_counts <- as.matrix(read.delim(file.path(TAB_DIR, sprintf("collapsed_counts_day17day25_%s.tsv", DB)),
                                   row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(TAB_DIR, "samples_day17day25.tsv"), row.names = 1,
                   stringsAsFactors = FALSE)[colnames(all_counts), ]

# Taxon sets. S1 is the old joint fit's rule: present in >= 2 samples of either
# arm, and seen in every Arm x Timepoint cell.
pres <- all_counts > 0
in_arm <- sapply(c("IFA", "IFA_IR"), function(a) rowSums(pres[, meta$Arm == a]) >= 2)
cell <- paste(meta$Arm, meta$Timepoint)
every_cell <- sapply(unique(cell), function(k) rowSums(pres[, cell == k, drop = FALSE]) > 0)
S1 <- rowSums(in_arm) > 0 & rowSums(!every_cell) == 0
SETS <- list(S1 = S1, S2 = S1 & rowSums(pres) >= 8, S3 = rowSums(pres) >= 4)
message(sprintf("[linda] %s taxa: %s", DB, paste(names(SETS), sapply(SETS, sum), collapse = ", ")))

# Fits: name -> reference arm, reference day, or "centred".
FITS <- list(IFA_Day17 = c("IFA", "Day17"), IFA_IR_Day17 = c("IFA_IR", "Day17"),
             IFA_Day25 = c("IFA", "Day25"), centred = NULL)
ROWS <- data.frame(
  row  = c("Day 17 to 25 in IFA", "Day 17 to 25 in IFA+IR", "IFA+IR vs IFA at Day 17",
           "IFA+IR vs IFA at Day 25", "Interaction (IFA+IR shift minus IFA shift)",
           "Average Day 17 to 25 shift", "Average IFA+IR vs IFA"),
  fit  = c("IFA_Day17", "IFA_IR_Day17", "IFA_Day17", "IFA_Day25", "IFA_Day17", "centred", "centred"),
  coef = c("TimepointDay25", "TimepointDay25", "ArmIFA_IR", "ArmIFA_IR",
           "TimepointDay25:ArmIFA_IR", "Time", "Radiation"),
  stringsAsFactors = FALSE)

fit_one <- function(counts, name, zero) {
  m <- meta
  if (name == "centred") {
    m$Time      <- ifelse(m$Timepoint == "Day25", 0.5, -0.5)
    m$Radiation <- ifelse(m$Arm == "IFA_IR", 0.5, -0.5)
    f <- "~ Time * Radiation + (1|MouseID)"
  } else {
    m$Arm       <- relevel(factor(m$Arm), FITS[[name]][1])
    m$Timepoint <- relevel(factor(m$Timepoint), FITS[[name]][2])
    f <- FORMULA
  }
  msg <- capture.output(o <- if (identical(zero, "adaptive"))
    linda(counts, m, f, type = "count", adaptive = TRUE, prev.cut = 0, lib.cut = 1, winsor.quan = NULL)
  else
    linda(counts, m, f, type = "count", adaptive = FALSE, imputation = FALSE, pseudo.cnt = zero,
          prev.cut = 0, lib.cut = 1, winsor.quan = NULL))
  o$zero_handling <- if (!identical(zero, "adaptive")) paste("pseudo-count", zero) else
    if (any(grepl("Imputation", msg))) "imputation" else "pseudo-count 0.5"
  o
}

row_table <- function(fits) do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  f <- fits[[ROWS$fit[k]]]
  o <- f$output[[ROWS$coef[k]]]
  # wald$beta has unnamed columns in output order; LinDA's estimate is the raw
  # coefficient minus the row's mode, which pins the alignment.
  raw <- f$wald$beta[ROWS$coef[k], ]
  stopifnot(length(raw) == nrow(o),
            isTRUE(all.equal(unname(raw) - f$wald$bias[match(ROWS$coef[k], rownames(f$wald$beta))],
                             o$log2FoldChange)))
  data.frame(row = ROWS$row[k], taxon = rownames(o), log2FC = o$log2FoldChange, lfcSE = o$lfcSE,
             df = o$df, p = o$pvalue, q = p.adjust(o$pvalue, "BH"), raw = unname(raw),
             row.names = NULL, stringsAsFactors = FALSE)
}))

results <- list(); diag <- list()
for (s in names(SETS)) {
  counts <- all_counts[SETS[[s]], , drop = FALSE]
  n <- ncol(counts)
  runs <- c(list(primary = "adaptive"), setNames(as.list(PSEUDO), paste0("pc", PSEUDO)))
  fits <- lapply(runs, function(z) setNames(lapply(names(FITS), fit_one, counts = counts, zero = z), names(FITS)))
  res <- row_table(fits$primary)
  for (r in names(runs)[-1]) {
    x <- row_table(fits[[r]]); stopifnot(identical(x$taxon, res$taxon))
    res[[paste0("q_", r)]] <- x$q
  }
  res$passed_pseudo <- Reduce(`&`, lapply(c("q", paste0("q_pc", PSEUDO)), function(c) res[[c]] < ALPHA))
  res$set <- s
  # The adaptive switch decides per fit, so zero handling can differ between rows.
  res$zero_handling <- sapply(fits$primary, `[[`, "zero_handling")[ROWS$fit[match(res$row, ROWS$row)]]
  results[[s]] <- res
  # One diagnostic row per run: the primary (adaptive) run decides the verdict;
  # the fixed pseudo-count runs show whether the adaptive switch moved the mode.
  diag[[s]] <- do.call(rbind, lapply(names(runs), function(run) do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
    f <- fits[[run]][[ROWS$fit[k]]]
    raw <- f$wald$beta[ROWS$coef[k], ]
    data.frame(set = s, run = run, n_taxa = nrow(counts), row = ROWS$row[k],
               zero_handling = f$zero_handling,
               linda_mode = f$wald$bias[match(ROWS$coef[k], rownames(f$wald$beta))],
               median = median(raw), hsm = mlv(sqrt(n) * raw, method = "hsm") / sqrt(n),
               sd_raw = sd(raw), calls = sum(res$row == ROWS$row[k] & res$q < ALPHA),
               calls_pass_pseudo = sum(res$row == ROWS$row[k] & res$passed_pseudo))
  }))))
}
res  <- do.call(rbind, results)
diag <- do.call(rbind, diag)
diag_all <- diag
diag <- diag[diag$run == "primary", ]
TAG <- paste0("LinDA_day17day25_", DB)
write.table(res[, setdiff(names(res), "raw")], file.path(TAB_DIR, paste0(TAG, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(res[, c("set", "row", "taxon", "raw")], file.path(TAB_DIR, paste0(TAG, "_raw_coefficients.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Pass criteria (fixed before fitting).
verdict <- do.call(rbind, lapply(ROWS$row, function(r) {
  d <- diag[diag$row == r, ]; s1 <- d[d$set == "S1", ]
  data.frame(row = r, mode_S1 = s1$linda_mode, mode_S2 = d$linda_mode[d$set == "S2"],
             mode_S3 = d$linda_mode[d$set == "S3"], mode_range = diff(range(d$linda_mode)),
             S1_mode_minus_median = s1$linda_mode - s1$median,
             S1_mode_minus_hsm = s1$linda_mode - s1$hsm,
             pass_stable = diff(range(d$linda_mode)) <= TOL,
             pass_peak = abs(s1$linda_mode - s1$median) <= TOL & abs(s1$linda_mode - s1$hsm) <= TOL)
}))
verdict$pass <- verdict$pass_stable & verdict$pass_peak
write.table(diag_all, file.path(TAB_DIR, paste0(TAG, "_mode_diagnostic.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(verdict, file.path(TAB_DIR, paste0(TAG, "_verdict.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

dens <- res[res$set == "S1", c("row", "raw")]
dens$row <- factor(dens$row, levels = ROWS$row)
d1 <- diag[diag$set == "S1", ]; d1$row <- factor(d1$row, levels = ROWS$row)
p <- ggplot(dens, aes(raw)) +
  geom_density(fill = "grey85", colour = "grey30") + geom_rug(alpha = 0.6) +
  geom_vline(data = diag, aes(xintercept = linda_mode, linetype = set), colour = "#C0392B") +
  geom_vline(data = d1, aes(xintercept = median), colour = "#2471A3", linetype = "dashed") +
  facet_wrap(~ factor(row, levels = ROWS$row), scales = "free", ncol = 2) +
  labs(x = "uncorrected coefficient (log2)", y = "density",
       title = sprintf("LinDA coefficients before bias correction, %s, S1 = %d taxa", DB, sum(SETS$S1)),
       subtitle = "red: mode LinDA subtracts (S1 solid, S2 dashed, S3 dotted); blue dashed: S1 median") +
  scale_linetype_manual(values = c(S1 = "solid", S2 = "dashed", S3 = "dotted")) +
  theme_bw(base_size = 9)
ggsave(file.path(FIG_DIR, paste0(TAG, "_mode_diagnostic.png")), p, width = 10, height = 11, dpi = 150)

options(width = 200)
print(diag[, c("set", "n_taxa", "row", "zero_handling", "linda_mode", "median", "hsm", "sd_raw", "calls", "calls_pass_pseudo")],
      digits = 3, right = FALSE)
print(verdict, digits = 3, right = FALSE)
# The same two criteria on each fixed pseudo-count run (supplementary; not the verdict).
for (run in paste0("pc", PSEUDO)) {
  d <- diag_all[diag_all$run == run, ]
  v <- sapply(ROWS$row, function(r) { x <- d[d$row == r, ]; s1 <- x[x$set == "S1", ]
    diff(range(x$linda_mode)) <= TOL & abs(s1$linda_mode - s1$median) <= TOL & abs(s1$linda_mode - s1$hsm) <= TOL })
  message(sprintf("[linda] %s: %d of %d rows pass (%s)", run, sum(v), length(v), paste(which(v), collapse = ",")))
  print(reshape(d[, c("set", "row", "linda_mode")], idvar = "row", timevar = "set", direction = "wide"), digits = 3, right = FALSE)
}
