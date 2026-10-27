# Shared by run_maaslin3_rows.R and check_maaslin3_median.R: inputs, the
# feature filter, the report's twelve rows and one MaAsLin3 1.5.7 fit per pair
# of reference levels. The rows and families are run_aldex3_factorial.R's.
# The caller sets THIS_DIR and `args` before sourcing.
suppressPackageStartupMessages(library(maaslin3))
stopifnot(packageVersion("maaslin3") == "1.5.7")

OUT_DIR <- file.path(dirname(THIS_DIR), "comparison-Treatment-x-Radiation-factorial")
MA_DIR  <- file.path(OUT_DIR, "tables", "maaslin3")
FIG_DIR <- file.path(OUT_DIR, "figures")
ALPHA <- 0.05
TREATS <- c("Control", "Inosine", "IFA"); RADS <- c("None", "IR")
FORMULA <- "~ Treatment * Radiation + log10_reads"

# >= 3 of 24 is MaAsLin3's min_prevalence 0.10; it keeps at least 98.8% of every
# mouse's reads, so TSS is close to the full library. Applied once to the full
# data, so leave-one-out refits keep the same features.
i <- match("--min-samples", args)
MIN_SAMPLES <- if (is.na(i)) 3L else as.integer(args[i + 1L])
TAG <- sprintf("Greengenes2_Day17_min%d", MIN_SAMPLES)

counts <- read.delim(file.path(MA_DIR, "counts_Greengenes2.tsv"), row.names = 1, check.names = FALSE)
meta   <- read.delim(file.path(MA_DIR, "metadata.tsv"), row.names = 1, stringsAsFactors = FALSE)
stopifnot(identical(rownames(counts), rownames(meta)),
          identical(as.vector(table(factor(meta$Treatment, TREATS), factor(meta$Radiation, RADS))), rep(4L, 6)))
counts <- counts[, colSums(counts > 0) >= MIN_SAMPLES, drop = FALSE]
n <- nrow(counts)
message(sprintf("[maaslin3] %s: %d features x %d samples", TAG, ncol(counts), n))

# Each report row is the direct coefficient `term` of the fit whose reference
# levels are "<treatment> <radiation>"; `cells` are the groups it compares.
# BH runs within `family`, over its rows together.
row_ <- function(row, stratum, family, fit, term, cells) list(row = row, stratum = stratum, family = family,
                                                              fit = fit, term = term, cells = cells)
cl <- function(t, r) expand.grid(x = t, y = r, stringsAsFactors = FALSE)
ROWS <- list(
  row_("Inosine vs Control", "unirradiated", "treatment_unirradiated", "Control None", "TreatmentInosine", cl(c("Inosine", "Control"), "None")),
  row_("IFA vs Control", "unirradiated", "treatment_unirradiated", "Control None", "TreatmentIFA", cl(c("IFA", "Control"), "None")),
  row_("IFA vs Inosine", "unirradiated", "treatment_unirradiated", "Inosine None", "TreatmentIFA", cl(c("IFA", "Inosine"), "None")),
  row_("Inosine vs Control", "irradiated", "treatment_irradiated", "Control IR", "TreatmentInosine", cl(c("Inosine", "Control"), "IR")),
  row_("IFA vs Control", "irradiated", "treatment_irradiated", "Control IR", "TreatmentIFA", cl(c("IFA", "Control"), "IR")),
  row_("IFA vs Inosine", "irradiated", "treatment_irradiated", "Inosine IR", "TreatmentIFA", cl(c("IFA", "Inosine"), "IR")),
  row_("Radiation within Control", "", "radiation_Control", "Control None", "RadiationIR", cl("Control", RADS)),
  row_("Radiation within Inosine", "", "radiation_Inosine", "Inosine None", "RadiationIR", cl("Inosine", RADS)),
  row_("Radiation within IFA", "", "radiation_IFA", "IFA None", "RadiationIR", cl("IFA", RADS)),
  row_("Interaction: Inosine vs Control", "", "interaction", "Control None", "TreatmentInosine:RadiationIR", cl(c("Inosine", "Control"), RADS)),
  row_("Interaction: IFA vs Control", "", "interaction", "Control None", "TreatmentIFA:RadiationIR", cl(c("IFA", "Control"), RADS)),
  row_("Interaction: IFA vs Inosine", "", "interaction", "Inosine None", "TreatmentIFA:RadiationIR", cl(c("IFA", "Inosine"), RADS)))
FITS <- unique(vapply(ROWS, `[[`, "", "fit"))
ROWTAB <- do.call(rbind, lapply(ROWS, function(r) data.frame(row = r$row, stratum = r$stratum, family = r$family,
                                                             fit = r$fit, term = r$term, stringsAsFactors = FALSE)))
cell_of <- function(rows) paste(meta$Treatment[rows], meta$Radiation[rows])

fit <- function(fit_name, rows = seq_len(n), median_comparison = FALSE, abundance_only = FALSE) {
  set.seed(1)
  ref <- strsplit(fit_name, " ")[[1]]
  m <- meta[rows, , drop = FALSE]
  m$Treatment <- relevel(factor(m$Treatment, levels = TREATS), ref[1])
  m$Radiation <- relevel(factor(m$Radiation, levels = RADS), ref[2])
  r <- maaslin3(input_data = counts[rows, , drop = FALSE], input_metadata = m, output = tempfile(),
                formula = FORMULA, min_prevalence = 0, normalization = "TSS", transform = "LOG",
                evaluate_only = if (abundance_only) "abundance" else NULL,
                median_comparison_abundance = median_comparison, median_comparison_prevalence = FALSE,
                augment = TRUE, warn_prevalence = !abundance_only, cores = 1,
                plot_summary_plot = FALSE, plot_associations = FALSE, save_models = FALSE, verbosity = "ERROR")
  list(abundance = r$fit_data_abundance$results, prevalence = r$fit_data_prevalence$results)
}

# TRUE where `feature` is present in every one of the six cells, among `rows`.
# The abundance model is fitted on present samples only, so a feature absent
# from a whole cell makes it rank deficient: R drops aliased columns, and with
# four mice a group what is left can be one mouse against one on a single
# residual degree of freedom. Requiring every cell keeps each fit full rank;
# `k` is unused, kept for the shared rows_table() call.
identifiable <- function(feature, k, rows) {
  ok <- rep(TRUE, length(feature))
  for (cell in unique(cell_of(rows))) {
    ok <- ok & colSums(counts[rows[cell_of(rows) == cell], feature, drop = FALSE] > 0) > 0
  }
  ok
}
