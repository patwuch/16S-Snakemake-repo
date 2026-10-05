# Shared by run_maaslin3_rows.R and check_maaslin3_median.R: inputs, the
# feature filter, and one MaAsLin3 1.5.7 fit per reference pair.
# The caller sets THIS_DIR and `args` before sourcing.
suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")
source(file.path(THIS_DIR, "factorial_rows.R"))

OUT_DIR <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
MA_DIR  <- file.path(OUT_DIR, "tables", "maaslin3")
FIG_DIR <- file.path(OUT_DIR, "figures")

i <- match("--min-samples", args)
MIN_SAMPLES <- if (is.na(i)) 2L else as.integer(args[i + 1L])
TAG <- sprintf("Greengenes2_min%d", MIN_SAMPLES)
FORMULA <- "~ Diet * Condition + log10_reads"

counts <- read.delim(file.path(MA_DIR, "counts_Greengenes2.tsv"), row.names = 1, check.names = FALSE)
meta   <- read.delim(file.path(MA_DIR, "metadata.tsv"), row.names = 1, stringsAsFactors = FALSE)
stopifnot(identical(rownames(counts), rownames(meta)),
          identical(as.vector(table(factor(meta$Diet, DIETS), factor(meta$Condition, CONDS))),
                    c(3L, 3L, 2L, 4L, 4L, 3L)))
# Applied once to the full data (>= 2 of 19 is MaAsLin3's min_prevalence 0.10),
# so leave-one-out refits keep the same features.
counts <- counts[, colSums(counts > 0) >= MIN_SAMPLES, drop = FALSE]
n <- nrow(counts)
message(sprintf("[maaslin3] %s: %d features x %d samples", TAG, ncol(counts), n))

row_fits <- unique(ROWS[, c("d_ref", "c_ref")])
key <- function(d, c) paste(d, c)

# One fit with the given reference levels. Returns MaAsLin3's results tables,
# `abundance` and (unless abundance_only) `prevalence`.
fit <- function(d_ref, c_ref, rows = seq_len(n), median_comparison = FALSE, abundance_only = FALSE) {
  m <- meta[rows, , drop = FALSE]
  m$Diet      <- relevel(factor(m$Diet, levels = DIETS), d_ref)
  m$Condition <- relevel(factor(m$Condition, levels = CONDS), c_ref)
  r <- maaslin3(input_data = counts[rows, , drop = FALSE], input_metadata = m, output = tempfile(),
                formula = FORMULA, min_prevalence = 0, normalization = "TSS", transform = "LOG",
                evaluate_only = if (abundance_only) "abundance" else NULL,
                median_comparison_abundance = median_comparison, median_comparison_prevalence = FALSE,
                augment = TRUE, warn_prevalence = !abundance_only,
                plot_summary_plot = FALSE, plot_associations = FALSE,
                save_models = FALSE, verbosity = "ERROR")
  list(abundance = r$fit_data_abundance$results, prevalence = r$fit_data_prevalence$results)
}
