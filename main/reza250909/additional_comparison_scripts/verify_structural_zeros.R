#!/usr/bin/env Rscript
# Check that the taxa a factorial report lists as missing from a whole group are
# exactly ANCOM-BC2's own structural zeros.
#
# The factorial fits run with group = NULL, which switches off ancombc2's
# structural-zero detection. This refits with the flattened Group column as the
# only term, struc_zero = TRUE and neg_lb = FALSE (the rule the package
# recommends below ~30 samples per group: zero counts in every sample of the
# group), on exactly the taxa the report's prevalence rule keeps. It then
# compares ancombc2's zero_ind flags, group by group, with the
# <stem>_untestable_absent_from_a_cell.tsv table. Only the flags are used; the
# refit's estimates are discarded.
#
# Usage, inside containers/ancombc2-env.sif (ANCOMBC >= 2.14.0):
#   Rscript verify_structural_zeros.R --meta <tsv> --tested <tsv> --untestable <tsv>
#     [--filter Column=Value]
# Paths are relative to the study directory.

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
LOG_PREFIX <- "[struc_zero]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages(library(ANCOMBC))
require_ancombc("2.14.0")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
path_of <- function(p) file.path(STUDY_DIR, p)

meta <- read.delim(path_of(arg_val("--meta")), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = "")
filt <- arg_val("--filter")
if (!is.null(filt)) {
  kv <- strsplit(filt, "=", fixed = TRUE)[[1]]
  meta <- meta[meta[[kv[1]]] == kv[2], , drop = FALSE]
}
tested <- read.delim(path_of(arg_val("--tested")), check.names = FALSE, quote = "")$taxon
untest <- read.delim(path_of(arg_val("--untestable")), check.names = FALSE, quote = "")
taxa <- c(tested, untest$taxon)

counts <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))[, rownames(meta), drop = FALSE]
level <- if (length(COLLAPSE_LEVEL) > 1) COLLAPSE_LEVEL[["Greengenes2"]] else COLLAPSE_LEVEL
tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"))
counts <- collapse_to_level(counts[rowSums(counts) > 0, , drop = FALSE], tax_map, level)
missing <- setdiff(taxa, rownames(counts))
if (length(missing)) stop(length(missing), " report taxa not found after collapsing, e.g. ", missing[1])
counts <- counts[taxa, , drop = FALSE]

meta$Group <- factor(meta$Group)
res <- ANCOMBC::ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = meta,
                         fix_formula = "Group", group = "Group",
                         struc_zero = TRUE, neg_lb = FALSE, prv_cut = 0, lib_cut = 0,
                         pseudo_sens = FALSE, global = FALSE, pairwise = FALSE, dunnet = FALSE,
                         trend = FALSE, verbose = FALSE)
zi <- as.data.frame(res$zero_ind)
rownames(zi) <- zi$taxon
zi$taxon <- NULL
colnames(zi) <- sub("^structural_zero \\(Group = (.*)\\)$", "\\1", colnames(zi))

# The report's own flags: zero in every sample of the group.
cell_cols <- grep(" \\(n=\\d+\\)$", colnames(untest), value = TRUE)
ours <- matrix(FALSE, nrow = length(taxa), ncol = length(cell_cols),
               dimnames = list(taxa, sub(" \\(n=\\d+\\)$", "", cell_cols)))
ours[untest$taxon, ] <- as.matrix(untest[, cell_cols]) == 0
theirs <- as.matrix(zi[taxa, colnames(ours)]) == TRUE

agree <- identical(unname(ours), unname(theirs))
.say(sprintf("%d taxa, %d groups; ANCOM-BC2 flags %d taxon-group structural zeros in %d taxa; the report's table has %d in %d taxa; identical: %s",
             length(taxa), ncol(ours), sum(theirs), sum(rowSums(theirs) > 0), sum(ours), sum(rowSums(ours) > 0), agree))
if (!agree) {
  d <- which(ours != theirs, arr.ind = TRUE)
  print(head(data.frame(taxon = rownames(ours)[d[, 1]], group = colnames(ours)[d[, 2]],
                        ours = ours[d], ancombc2 = theirs[d]), 20))
  quit(status = 1)
}
