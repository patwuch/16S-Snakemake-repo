# Shared by run_linda_factorial.R and check_linda_mode.R: inputs, design
# checks, the zero-handling decision, one LinDA fit, and the report rows.
# The caller sets THIS_DIR and `args` before sourcing.

suppressPackageStartupMessages({ library(LinDA); library(modeest); library(ggplot2) })
stopifnot(packageVersion("LinDA") == "0.2.0")

STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Condition-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
DB      <- arg_val("--db", "Greengenes2")
MINPRES <- as.integer(arg_val("--min-presence", "6"))
ALPHA   <- 0.05
PSEUDO  <- c(0.1, 0.5, 1)
FORMULA <- "~ Diet * Condition"
source(file.path(THIS_DIR, "factorial_rows.R"))
TAG     <- paste0(DB, if (MINPRES != 6L) paste0("_minpres", MINPRES) else "")

counts <- as.matrix(read.delim(file.path(TAB_DIR, paste0("tested_counts_", TAG, ".tsv")),
                               row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)[colnames(counts), c("Group", "Diet", "Condition")]
design <- table(factor(meta$Diet, DIETS), factor(meta$Condition, CONDS))
stopifnot(!anyNA(meta$Diet), !anyNA(meta$Condition),
          identical(as.vector(design), c(3L, 3L, 2L, 4L, 4L, 3L)))
stopifnot(all(sapply(split(seq_len(ncol(counts)), paste(meta$Diet, meta$Condition)),
                     function(i) all(rowSums(counts[, i, drop = FALSE] > 0) > 0))))
message(sprintf("[linda] %s: %d taxa x %d samples", TAG, nrow(counts), ncol(counts)))

# LinDA's adaptive switch imputes when any coefficient of log depth on the
# design has p <= corr.cut. Which coefficients exist depends on the reference
# coding, so the switch picks differently across the six fits. Decide once
# over every pair of cells instead, and fix that choice for all fits.
CORR_CUT <- 0.1
depth_p <- pairwise.t.test(log(colSums(counts)), interaction(meta$Diet, meta$Condition),
                           p.adjust.method = "none", pool.sd = TRUE)$p.value
PRIMARY <- if (min(depth_p, na.rm = TRUE) <= CORR_CUT) "imputation" else 0.5
message(sprintf("[linda] depth vs cell: min pairwise p = %.3g -> primary zero handling: %s",
                min(depth_p, na.rm = TRUE), PRIMARY))

# `zero` is "imputation" or a fixed pseudo-count; the adaptive switch is off.
# `cols` drops samples for the leave-one-out check; the default fits all.
fit <- function(d_ref, c_ref, zero, cols = seq_len(ncol(counts))) {
  m <- meta[cols, , drop = FALSE]
  m$Diet      <- relevel(factor(m$Diet, levels = DIETS), d_ref)
  m$Condition <- relevel(factor(m$Condition, levels = CONDS), c_ref)
  imp <- identical(zero, "imputation")
  capture.output(o <- linda(counts[, cols, drop = FALSE], m, FORMULA, type = "count", adaptive = FALSE,
                            imputation = imp, pseudo.cnt = if (imp) 0.5 else zero,
                            prev.cut = 0, lib.cut = 1, winsor.quan = NULL))
  o$zero_handling <- if (imp) "imputation" else paste("pseudo-count", zero)
  o
}

# The same rows as contrasts of the Cellulose/Control fit's coefficients
# (Intercept, Inulin, RS, IR, Inulin:IR, RS:IR), for the Wald comparison.
L <- rbind(c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0), c(0, -1, 1, 0, 0, 0),
           c(0, 1, 0, 0, 1, 0), c(0, 0, 1, 0, 0, 1), c(0, -1, 1, 0, -1, 1),
           c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 1, 1, 0), c(0, 0, 0, 1, 0, 1),
           c(0, 0, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 1), c(0, 0, 0, 0, -1, 1))

# Every diet as reference under both conditions: the rows need five of these,
# the omnibus symmetry check needs all six.
refs <- expand.grid(d_ref = DIETS, c_ref = CONDS, stringsAsFactors = FALSE)
settings <- c(list(PRIMARY), as.list(PSEUDO))
