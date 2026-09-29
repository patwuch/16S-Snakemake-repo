#!/usr/bin/env Rscript
# MaAsLin 3 on the tisha251224 Diet x Condition design, Greengenes2 species.
# Exploratory: run alongside the ANCOM-BC2 factorial, not yet part of report.html.
#
# MaAsLin 3 fits two models per feature: a logistic model of presence/absence
# (prevalence) over all samples, and a linear model of log2 relative abundance
# over the samples where the feature is present (abundance). The prevalence
# model is what ANCOM-BC2 cannot do here: test taxa that are missing from whole
# groups.
#
# Two fits:
#   interaction  ~ Diet * Condition + log10_reads   (mirrors the ANCOM-BC2 model)
#   additive     ~ Diet + Condition + log10_reads   (irradiation averaged over
#                diets, 9 vs 10 samples: the best-powered prevalence question)
# Read depth is a covariate in both, as the authors require for prevalence
# models. Settings otherwise follow MaAsLin 3's defaults: TSS + log2, the
# abundance model tested against the median coefficient (compositionality),
# the prevalence model against zero, and augmentation against separation.
#
# Inputs from prep_maaslin3_input.R (run that first, in the core container).
#
# Usage:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/maaslin3-env.sif \
#     Rscript run_maaslin3_factorial.R

suppressPackageStartupMessages(library(maaslin3))

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
MA_DIR   <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial", "tables", "maaslin3")

# A taxon must be seen in at least 3 of 19 samples: "more than 15%" in
# MaAsLin 3's terms. One or two detections carry no prevalence contrast.
MIN_PREVALENCE <- 0.15
ALPHA <- 0.05

counts <- read.delim(file.path(MA_DIR, "counts_Greengenes2.tsv"), row.names = 1, check.names = FALSE)
meta   <- read.delim(file.path(MA_DIR, "metadata.tsv"), row.names = 1, stringsAsFactors = FALSE)
meta$Diet      <- factor(meta$Diet, levels = c("Cellulose", "Inulin", "RS"))
meta$Condition <- factor(meta$Condition, levels = c("Control", "IR"))
stopifnot(identical(rownames(counts), rownames(meta)))

FITS <- c(interaction = "~ Diet * Condition + log10_reads",
          additive    = "~ Diet + Condition + log10_reads")

for (fit in names(FITS)) {
  out <- file.path(MA_DIR, paste0("fit_", fit))
  message("[maaslin3] ", fit, ": ", FITS[[fit]])
  maaslin3(input_data = counts, input_metadata = meta, output = out,
           formula = FITS[[fit]], min_prevalence = MIN_PREVALENCE,
           normalization = "TSS", transform = "LOG",
           median_comparison_abundance = TRUE, median_comparison_prevalence = FALSE,
           augment = TRUE, warn_prevalence = TRUE, max_significance = ALPHA,
           plot_summary_plot = TRUE, plot_associations = FALSE,
           verbosity = "WARN")
}
