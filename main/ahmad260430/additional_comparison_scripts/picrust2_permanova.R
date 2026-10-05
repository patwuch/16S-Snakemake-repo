#!/usr/bin/env Rscript
# Whole-profile tests of PICRUSt2's predictions (picrust2/, run_picrust2.sh), one
# per level (MetaCyc pathways, KO, EC), on Bray-Curtis distances between samples'
# relative profiles, 9,999 permutations, seed 42; the same designs as each
# report's beta-diversity section.
#   factorial    Treatment * Radiation over the 24 Day 17 mice, Type II (each main
#                effect entering last); PERMDISP across the six groups. The
#                per-feature tests (ggpicrust2 LinDA) cannot model the interaction;
#                this does.
#   day17-day25  Timepoint within each arm's eight samples, unrestricted
#                permutations (as the taxon PERMANOVA: the 16 within-mouse swaps
#                cannot reach 0.05); PERMDISP between the timepoints per arm.
#
# Writes <report dir>/tables/PICRUSt2_permanova.tsv.
#
# Usage, in containers/ancombc2-env.sif (vegan):
#   Rscript picrust2_permanova.R factorial|day17-day25

suppressPackageStartupMessages(library(vegan))
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
PIC       <- file.path(STUDY_DIR, "picrust2", "pipeline")
NPERM <- 9999L
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
meta_all <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), row.names = 1, check.names = FALSE,
                       stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = character(0))
bray <- function(lv, samples) {
  ab <- as.matrix(read.delim(gzfile(file.path(PIC, FILES[[lv]])), row.names = 1, check.names = FALSE)[, samples])
  ab <- ab[rowSums(ab) > 0, , drop = FALSE]
  list(n = nrow(ab), d = vegdist(t(sweep(ab, 2, colSums(ab), "/")), method = "bray"))
}

REPORTS <- list(
  factorial = list(dir = "comparison-Treatment-x-Radiation-factorial", test = function(lv) {
    m <- meta_all[meta_all$Timepoint == "Day17", ]
    m <- data.frame(Treatment = factor(m$Treatment, c("Control", "IFA", "Inosine")),
                    Radiation = factor(m$Radiation, c("None", "IR")), row.names = rownames(m))
    b <- bray(lv, rownames(m))
    set.seed(42); a1 <- adonis2(b$d ~ Treatment * Radiation, data = m, permutations = NPERM, by = "terms")
    set.seed(42); a2 <- adonis2(b$d ~ Radiation * Treatment, data = m, permutations = NPERM, by = "terms")
    ad <- rbind(a2["Treatment", ], a1["Radiation", ], a1["Treatment:Radiation", ])
    set.seed(42); pt <- permutest(betadisper(b$d, interaction(m$Treatment, m$Radiation, drop = TRUE)), permutations = NPERM)
    data.frame(level = lv, n_features = b$n, arm = "all",
               analysis = c(rep("permanova_bray_curtis", 3), "permdisp_bray_curtis"),
               term = c("Treatment", "Radiation", "Treatment:Radiation", "group"),
               statistic = c(ad$F, pt$tab[1, "F"]), effect = c(ad$R2, NA), p = c(ad$`Pr(>F)`, pt$tab[1, "Pr(>F)"]))
  }),
  "day17-day25" = list(dir = "comparison-day17-day25", test = function(lv) {
    cells <- read.delim(file.path(STUDY_DIR, "comparison-day17-day25", "metadata_joint_cells.tsv"),
                        row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
    do.call(rbind, lapply(c("IFA", "IFA+IR"), function(a) {
      m <- cells[sub(" Day.*", "", cells$Group) == a, , drop = FALSE]
      m <- data.frame(Timepoint = factor(sub(".* ", "", m$Group), c("Day17", "Day25")), row.names = rownames(m))
      b <- bray(lv, rownames(m))
      set.seed(42); ad <- adonis2(b$d ~ Timepoint, data = m, permutations = NPERM)
      set.seed(42); pt <- permutest(betadisper(b$d, m$Timepoint), permutations = NPERM)
      data.frame(level = lv, n_features = b$n, arm = a,
                 analysis = c("permanova_bray_curtis", "permdisp_bray_curtis"), term = c("Timepoint", "Timepoint"),
                 statistic = c(ad$F[1], pt$tab[1, "F"]), effect = c(ad$R2[1], NA), p = c(ad$`Pr(>F)`[1], pt$tab[1, "Pr(>F)"]))
    }))
  }))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1 || !args[1] %in% names(REPORTS)) stop("argument: ", paste(names(REPORTS), collapse = " or "))
report <- REPORTS[[args[1]]]
out <- do.call(rbind, lapply(names(FILES), report$test))
write.table(out, file.path(STUDY_DIR, report$dir, "tables", "PICRUSt2_permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
