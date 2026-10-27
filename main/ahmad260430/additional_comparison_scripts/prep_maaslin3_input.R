#!/usr/bin/env Rscript
# Write the inputs maaslin3_common.R reads, for the Day 17 Treatment x Radiation
# factorial: the study BIOM table collapsed to Greengenes2 species, as plain
# TSVs (maaslin3-env has no BIOM reader, so this runs in the qiime2 core
# container). Ported from reza250909.
#   tables/maaslin3/counts_Greengenes2.tsv    samples x features, raw counts
#   tables/maaslin3/metadata.tsv              Treatment, Radiation, log10_reads
#   tables/maaslin3/feature_key.tsv           feature id -> Greengenes2 lineage
#
# No prevalence filter here: maaslin3_common.R sets it, so the filter is
# visible where the model is.
#
# Usage:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/qiime2-2025.10-amplicon-core.sif \
#     Rscript prep_maaslin3_input.R

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
MA_DIR    <- file.path(OUT_DIR, "tables", "maaslin3")
dir.create(MA_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_PREFIX <- "[maaslin3-prep]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages(library(biomformat))

otu  <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
meta <- read.delim(file.path(OUT_DIR, "tables", "metadata_day17.tsv"), row.names = 1,
                   check.names = FALSE, stringsAsFactors = FALSE)
otu  <- otu[, rownames(meta), drop = FALSE]

# MaAsLin 3's authors ask for read depth as a covariate: deeper samples detect
# more taxa, which the prevalence model would otherwise read as an effect.
# Depth is the full ASV library, before collapsing drops unannotated reads.
meta$log10_reads <- log10(colSums(otu))

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"))
counts  <- collapse_to_level(otu[rowSums(otu) > 0, , drop = FALSE], tax_map, COLLAPSE_LEVEL[["Greengenes2"]])

# Taxonomy strings carry ";" and spaces that R rewrites in column names, so
# features travel as F001... ids with a key back to the lineage.
ids <- sprintf("F%03d", seq_len(nrow(counts)))
write.table(data.frame(feature = ids, taxon = rownames(counts)),
            file.path(MA_DIR, "feature_key.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tab <- t(counts); colnames(tab) <- ids
write.table(data.frame(sample = rownames(tab), tab, check.names = FALSE),
            file.path(MA_DIR, "counts_Greengenes2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(sample = rownames(meta), meta[, c("Group", "Treatment", "Radiation", "log10_reads")]),
            file.path(MA_DIR, "metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
.say(sprintf("wrote %d samples x %d features to %s", ncol(counts), nrow(counts), MA_DIR))
