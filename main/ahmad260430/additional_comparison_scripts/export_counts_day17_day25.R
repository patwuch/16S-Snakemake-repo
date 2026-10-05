#!/usr/bin/env Rscript
# Writes the collapsed Day17/Day25 count table (both IFA arms, 16 samples) and
# its sample table, for methods whose containers cannot read BIOM (LinDA,
# MaAsLin3). Collapse level and mouse pairing as in the old ANCOM-BC2 scripts.
#
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript export_counts_day17_day25.R [--db Greengenes2|Silva138]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
LOG_PREFIX <- "[export_day17_day25]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages(library(biomformat))

args <- commandArgs(trailingOnly = TRUE)
i <- match("--db", args)
db <- if (is.na(i)) "Greengenes2" else args[i + 1L]

ARMS <- c(IFA = "IFA", IFA_IR = "IFA+IR")
mouse_id <- function(sid) {
  m <- regmatches(sid, regexec("^E\\d+_(\\d+)\\.(\\d+)$", sid))[[1]]
  if (length(m) == 0) m <- regmatches(sid, regexec("^(\\d+)\\.(\\d+)colon$", sid))[[1]]
  if (length(m) == 0) stop("Cannot parse a mouse id out of sample '", sid, "'")
  paste(m[2], m[3], sep = "_")
}

otu_all  <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
meta_all <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE, quote = "", comment.char = "")
meta <- meta_all[meta_all$Group %in% ARMS & rownames(meta_all) %in% colnames(otu_all), , drop = FALSE]
meta <- data.frame(sampleid = rownames(meta), Arm = names(ARMS)[match(meta$Group, ARMS)],
                   Timepoint = meta$Timepoint,
                   MouseID = vapply(rownames(meta), mouse_id, character(1)))
stopifnot(all(table(meta$MouseID, meta$Timepoint) == 1), nrow(meta) == 16)

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(db, "_taxonomy.tsv")))
counts <- otu_all[, meta$sampleid, drop = FALSE]
counts <- counts[rowSums(counts) > 0, , drop = FALSE]
counts <- collapse_to_level(counts, tax_map, COLLAPSE_LEVEL[[db]])

tab_dir <- file.path(STUDY_DIR, "comparison-day17-day25", "tables")
write.table(data.frame(taxon = rownames(counts), counts, check.names = FALSE),
            file.path(tab_dir, sprintf("collapsed_counts_day17day25_%s.tsv", db)),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(meta, file.path(tab_dir, "samples_day17day25.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
.say(sprintf("%s: %d taxa x %d samples", db, nrow(counts), ncol(counts)))
