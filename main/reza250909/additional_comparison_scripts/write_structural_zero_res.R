#!/usr/bin/env Rscript
# Write a figure-only copy of a DA panel's .res with its structural zeros added.
#
# A panel's .res holds the taxa its fit tested. Structural zeros (taxa absent
# from every sample of one group) were never tested, so they have no node in
# that tree. For each panel in the plan this rebuilds the tree over the tested
# taxa plus the panel's structural zeros, with the same write_res() and
# canonical_lineage() the fits use, so node keys match by construction; carries
# the original markers over unchanged; and writes a sidecar naming the
# structural-zero nodes for the renderer's struc_zeros option.
#
# The original .res is not modified: tables and report checks read it, and the
# augmented copy exists only to draw figures.
#
# Plan: TSV, one row per (panel, structural-zero taxon), columns
#   res_in res_out sidecar_out fit_tsv meta filter select_col select_values
#   class_col first second taxon present_in annotation
# `filter` is "" or Column=Value; `select_values` is comma-separated. A panel
# with no structural zeros has no rows and is not written.
#
# Usage, inside the qiime2 container:
#   Rscript write_structural_zero_res.R <plan.tsv>

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))
LOG_PREFIX <- "[sz-res]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
source(file.path(REPO_ROOT, "scripts", "write_res.R"))

plan <- read.delim(commandArgs(trailingOnly = TRUE)[1], check.names = FALSE, quote = "", colClasses = "character",
                   stringsAsFactors = FALSE, na.strings = character(0))
counts_all <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
level <- if (length(COLLAPSE_LEVEL) > 1) COLLAPSE_LEVEL[["Greengenes2"]] else COLLAPSE_LEVEL
tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"))
collapsed <- collapse_to_level(counts_all[rowSums(counts_all) > 0, , drop = FALSE], tax_map, level)
key_of <- function(t) canonical_lineage(paste(trimws(strsplit(t, ";")[[1]]), collapse = "|"))

for (out in unique(plan$res_out)) {
  p <- plan[plan$res_out == out, , drop = FALSE]
  s <- p[1, ]
  meta <- read.delim(file.path(STUDY_DIR, s$meta), row.names = 1, check.names = FALSE,
                     stringsAsFactors = FALSE, quote = "", na.strings = character(0))
  if (nzchar(s$filter)) {
    kv <- strsplit(s$filter, "=", fixed = TRUE)[[1]]
    meta <- meta[meta[[kv[1]]] == kv[2], , drop = FALSE]
  }
  meta <- meta[meta[[s$select_col]] %in% strsplit(s$select_values, ",", fixed = TRUE)[[1]], , drop = FALSE]
  tested <- read.delim(file.path(STUDY_DIR, s$fit_tsv), check.names = FALSE, quote = "")$taxon
  # A fit that tests a taxon missing from a whole group (a single-factor per-arm
  # fit does, and returns no SE for it) already has its node; add only the rest.
  taxa <- unique(c(tested, p$taxon))
  cs <- collapsed[taxa, rownames(meta), drop = FALSE]
  cs <- cs[, colSums(cs) > 0, drop = FALSE]
  lineage_map <- lapply(taxa, key_of)

  # Markers exactly as the fit wrote them.
  orig <- read.delim(file.path(STUDY_DIR, s$res_in), header = FALSE, quote = "", colClasses = "character",
                     col.names = c("lineage", "value", "class", "lda", "pvalue"), na.strings = character(0))
  mk <- orig[nzchar(orig$class), c("lineage", "class", "lda", "pvalue"), drop = FALSE]
  mk$pvalue[mk$pvalue == "-"] <- NA
  res_out <- file.path(STUDY_DIR, out)
  write_res(cs, as.character(meta[colnames(cs), s$class_col]), lineage_map,
            if (nrow(mk)) mk else NULL, res_out)

  # The tested part of the tree must be the fit's own: same nodes, same markers.
  new <- read.delim(res_out, header = FALSE, quote = "", colClasses = "character",
                    col.names = c("lineage", "value", "class", "lda", "pvalue"), na.strings = character(0))
  if (!all(orig$lineage %in% new$lineage))
    stop(out, ": augmented tree lost nodes of ", s$res_in)
  m_new <- new[nzchar(new$class), c("lineage", "class", "lda", "pvalue")]
  m_old <- orig[nzchar(orig$class), c("lineage", "class", "lda", "pvalue")]
  if (!isTRUE(all.equal(m_new[order(m_new$lineage), ], m_old[order(m_old$lineage), ], check.attributes = FALSE)))
    stop(out, ": markers changed in the augmented copy")

  keys <- vapply(p$taxon, function(t) paste(key_of(t), collapse = "."), character(1))
  if (any(keys %in% m_new$lineage))
    stop(out, ": a structural zero is also a marked hit: ", keys[keys %in% m_new$lineage][1])
  writeLines(c(paste("#classes", s$first, s$second, sep = "\t"),
               paste(keys, p$present_in, p$annotation, sep = "\t")),
             file.path(STUDY_DIR, s$sidecar_out))
  .say(sprintf("%s: %d tested + %d structural-zero taxa, %d markers kept", basename(out),
               length(tested), nrow(p), nrow(m_new)))
}
