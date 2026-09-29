#!/usr/bin/env Rscript
# Does the interaction show up in which taxa are present, rather than in how
# abundant they are?
#
# The ANCOM-BC2 interaction (the forest plot) can only use taxa seen in every
# cell of the design: a taxon absent from a whole cell has no log abundance
# there. This asks the same interaction question of presence/absence, over the
# DA analysis's full prevalence-filtered taxon set: the taxa ANCOM-BC2 tests
# plus the ones it had to leave out. The set is fixed before looking at where
# anything is missing, so the test is not conditioned on its own outcome.
#
#  * PERMANOVA on binary Jaccard distances, A * B, Type II (each term entering
#    last, as run_diversity_factorial.R does it), 9,999 permutations, seed 42.
#  * PERMDISP across the cells of A x B, since PERMANOVA cannot separate a shift
#    in which taxa are present from a difference in how variable that is.
#  * --firth: per taxon, Firth-penalised logistic regression of presence on
#    A * B (logistf), likelihood-ratio p for the interaction, BH across taxa.
#    Firth, because a cell with the taxon in no sample separates ordinary
#    logistic regression. Only worth it with enough samples per cell.
#
# Identical in every study; the study-specific parts are the arguments.
#
# Usage, inside containers/ancombc2-env.sif (from the study's scripts folder):
#   Rscript presence_interaction.R --comparison <folder> --meta <tsv> [--filter Column=Value]
#     --factors A,B --levels-a L1,L2[,L3] --levels-b L1,L2 --stem <table stem> [--firth]
# Paths are relative to the study directory; --stem names the tested-taxa table
# (<stem>.tsv) and the untestable table (<stem>_untestable_absent_from_a_cell.tsv).

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
LOG_PREFIX <- "[presence]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages({ library(vegan); library(biomformat) })

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
NPERM <- 9999L
comp  <- arg_val("--comparison")
TAB   <- file.path(STUDY_DIR, comp, "tables")
stem  <- arg_val("--stem")
fac   <- strsplit(arg_val("--factors"), ",")[[1]]
lv_a  <- strsplit(arg_val("--levels-a"), ",")[[1]]
lv_b  <- strsplit(arg_val("--levels-b"), ",")[[1]]
firth <- "--firth" %in% args

meta <- read.delim(file.path(STUDY_DIR, arg_val("--meta")), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = "")
filt <- arg_val("--filter")
if (!is.null(filt)) {
  kv <- strsplit(filt, "=", fixed = TRUE)[[1]]
  meta <- meta[meta[[kv[1]]] == kv[2], , drop = FALSE]
}
meta[[fac[1]]] <- factor(meta[[fac[1]]], levels = lv_a)
meta[[fac[2]]] <- factor(meta[[fac[2]]], levels = lv_b)
if (anyNA(meta[fac])) stop("factor values outside the declared levels")

tested <- read.delim(file.path(TAB, paste0(stem, ".tsv")), check.names = FALSE, quote = "")$taxon
untest <- read.delim(file.path(TAB, paste0(stem, "_untestable_absent_from_a_cell.tsv")), check.names = FALSE, quote = "")$taxon
taxa <- c(tested, untest)

counts <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
counts <- counts[, rownames(meta), drop = FALSE]
level <- if (length(COLLAPSE_LEVEL) > 1) COLLAPSE_LEVEL[["Greengenes2"]] else COLLAPSE_LEVEL
tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"))
counts <- collapse_to_level(counts[rowSums(counts) > 0, , drop = FALSE], tax_map, level)
missing <- setdiff(taxa, rownames(counts))
if (length(missing)) stop(length(missing), " DA taxa not found after collapsing, e.g. ", missing[1])
pres <- t(counts[taxa, , drop = FALSE] > 0) * 1
.say(sprintf("%d samples x %d taxa (%d tested + %d untestable); %.1f taxa present per sample on average",
             nrow(pres), ncol(pres), length(tested), length(untest), mean(rowSums(pres))))

# ---- PERMANOVA, binary Jaccard, Type II ----
d <- vegdist(pres, method = "jaccard", binary = TRUE)
f1 <- as.formula(paste("d ~", fac[1], "*", fac[2]))
f2 <- as.formula(paste("d ~", fac[2], "*", fac[1]))
set.seed(42); ad1 <- adonis2(f1, data = meta, permutations = NPERM, by = "terms")
set.seed(42); ad2 <- adonis2(f2, data = meta, permutations = NPERM, by = "terms")
ia <- paste(fac, collapse = ":")
ad <- rbind(ad2[fac[1], ], ad1[fac[2], ], ad1[ia, ])
rownames(ad) <- c(fac, ia)
out <- data.frame(analysis = "permanova_binary_jaccard", term = rownames(ad), n = nrow(pres),
                  taxa = ncol(pres), statistic = ad$F, effect = ad$R2, p = ad$`Pr(>F)`)

# ---- PERMDISP across cells ----
cell <- interaction(meta[[fac[1]]], meta[[fac[2]]], drop = TRUE)
bd <- betadisper(d, cell)
set.seed(42); pt <- permutest(bd, permutations = NPERM)
out <- rbind(out, data.frame(analysis = "permdisp_binary_jaccard", term = "cell", n = nrow(pres),
                             taxa = ncol(pres), statistic = pt$tab[1, "F"], effect = NA_real_,
                             p = pt$tab[1, "Pr(>F)"]))
# Mean number of these taxa present per sample, by cell: the descriptive the
# report quotes beside the test.
bycell <- aggregate(list(taxa_present = rowSums(pres)), meta[fac], mean)
bycell$n <- aggregate(list(n = rowSums(pres)), meta[fac], length)$n
write.table(bycell, file.path(TAB, paste0(stem, "_presence_by_cell.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
path <- file.path(TAB, paste0(stem, "_presence_interaction.tsv"))
write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
.say("Wrote ", path)
print(out, row.names = FALSE)

# ---- per-taxon Firth logistic regression ----
if (firth) {
  suppressPackageStartupMessages(library(logistf))
  fx <- as.formula(paste("y ~", fac[1], "*", fac[2]))
  f0 <- as.formula(paste("y ~", fac[1], "+", fac[2]))
  rows <- lapply(taxa, function(t) {
    dd <- data.frame(y = pres[, t], meta[fac])
    full <- logistf(fx, data = dd)
    # Penalised likelihood-ratio test of the interaction, all its columns at once.
    p <- tryCatch(anova(full, logistf(f0, data = dd))$pval, error = function(e) NA_real_)
    ic <- grep(":", names(coef(full)), value = TRUE)
    data.frame(taxon = t, untestable = t %in% untest, n_present = sum(dd$y),
               interaction_coef = paste(ic, collapse = ";"),
               interaction_logOR = paste(sprintf("%.3f", coef(full)[ic]), collapse = ";"),
               p_interaction = p, check.names = FALSE)
  })
  fr <- do.call(rbind, rows)
  fr$q_interaction <- p.adjust(fr$p_interaction, method = "BH")
  fr <- fr[order(fr$p_interaction), ]
  fpath <- file.path(TAB, paste0(stem, "_presence_firth.tsv"))
  write.table(fr, fpath, sep = "\t", quote = FALSE, row.names = FALSE)
  .say(sprintf("Firth: %d taxa, %d with raw p < 0.05, %d with q < 0.05, smallest q %.3g; wrote %s",
               nrow(fr), sum(fr$p_interaction < 0.05, na.rm = TRUE), sum(fr$q_interaction < 0.05, na.rm = TRUE),
               min(fr$q_interaction, na.rm = TRUE), fpath))
}
