#!/usr/bin/env Rscript
# Do the two arms differ in which taxa appear or vanish between Day 17 and Day 25?
#
# The joint ANCOM-BC2 model's Timepoint x Arm term (the forest plot) can only use
# taxa seen in all four arm x timepoint cells. This asks the same question of
# presence/absence, over the joint model's full prevalence-filtered taxon set:
# the taxa it tests plus the ones it had to leave out, fixed before looking at
# where anything is missing.
#
# The design is paired, so the unit is the mouse: for each mouse, the binary
# Jaccard distance between its Day 17 and Day 25 samples (how much of its
# detected taxon set turns over), plus the counts of taxa gained and lost. The
# arms are compared by the difference in mean turnover, with an exact
# permutation test over all 70 ways to label four of the eight mice IFA. With
# 4 mice per arm the smallest two-sided p is 2/70 = 0.029.
#
# Usage, inside containers/ancombc2-env.sif:
#   Rscript presence_turnover_day17_day25.R

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-day17-day25")
TAB       <- file.path(OUT_DIR, "tables")
LOG_PREFIX <- "[turnover]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages({ library(vegan); library(biomformat) })

stem <- "ANCOMBC2_day17day25_joint_Greengenes2"
meta <- read.delim(file.path(OUT_DIR, "metadata_joint_cells.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
meta$Arm <- sub(" Day\\d+$", "", meta$Group)
meta$Timepoint <- sub("^.* ", "", meta$Group)
tested <- read.delim(file.path(TAB, paste0(stem, "_armrefIFA.tsv")), check.names = FALSE, quote = "")$taxon
untest <- read.delim(file.path(TAB, paste0(stem, "_untestable_absent_from_a_cell.tsv")), check.names = FALSE, quote = "")$taxon
taxa <- c(tested, untest)

counts <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))[, rownames(meta), drop = FALSE]
tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", "Greengenes2_taxonomy.tsv"))
counts <- collapse_to_level(counts[rowSums(counts) > 0, , drop = FALSE], tax_map, COLLAPSE_LEVEL[["Greengenes2"]])
missing <- setdiff(taxa, rownames(counts))
if (length(missing)) stop(length(missing), " joint-model taxa not found after collapsing, e.g. ", missing[1])
pres <- counts[taxa, , drop = FALSE] > 0

mice <- unique(meta$MouseID)
per <- do.call(rbind, lapply(mice, function(m) {
  s17 <- rownames(meta)[meta$MouseID == m & meta$Timepoint == "Day17"]
  s25 <- rownames(meta)[meta$MouseID == m & meta$Timepoint == "Day25"]
  if (length(s17) != 1 || length(s25) != 1) stop("mouse ", m, " does not have exactly one sample per timepoint")
  a <- pres[, s17]; b <- pres[, s25]
  data.frame(mouse = m, arm = unique(meta$Arm[meta$MouseID == m]), day17 = sum(a), day25 = sum(b),
             gained = sum(!a & b), lost = sum(a & !b), shared = sum(a & b),
             turnover = as.numeric(vegdist(rbind(a, b) * 1, method = "jaccard", binary = TRUE)))
}))
per <- per[order(per$arm, per$mouse), ]
arms <- sort(unique(per$arm))
if (length(arms) != 2 || any(table(per$arm) != 4)) stop("expected two arms of four mice")

stat <- function(lab) mean(per$turnover[lab == arms[2]]) - mean(per$turnover[lab == arms[1]])
obs <- stat(per$arm)
splits <- combn(nrow(per), 4)
null <- apply(splits, 2, function(i) { lab <- rep(arms[1], nrow(per)); lab[i] <- arms[2]; stat(lab) })
p <- mean(abs(null) >= abs(obs) - 1e-12)

write.table(per, file.path(TAB, paste0(stem, "_presence_turnover_by_mouse.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
summ <- data.frame(taxa = length(taxa), tested = length(tested), untestable = length(untest),
                   arm_a = arms[1], arm_b = arms[2],
                   turnover_a = mean(per$turnover[per$arm == arms[1]]), turnover_b = mean(per$turnover[per$arm == arms[2]]),
                   gained_a = mean(per$gained[per$arm == arms[1]]), gained_b = mean(per$gained[per$arm == arms[2]]),
                   lost_a = mean(per$lost[per$arm == arms[1]]), lost_b = mean(per$lost[per$arm == arms[2]]),
                   diff_b_minus_a = obs, n_splits = ncol(splits), p_exact = p)
write.table(summ, file.path(TAB, paste0(stem, "_presence_turnover.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(per, row.names = FALSE)
print(summ, row.names = FALSE)
