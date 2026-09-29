#!/usr/bin/env Rscript
# Community-level Diet x Radiation test of PICRUSt2's predicted profiles, one per
# level (MetaCyc pathways, KO, EC): PERMANOVA on Bray-Curtis distances between
# samples' relative profiles, Diet * Radiation, Type II (each main effect entering
# last), 9,999 permutations, seed 42; PERMDISP across the four groups. The
# per-feature tests (ggpicrust2 LinDA, run_ggpicrust2_linda.R) cannot model the
# interaction; this does.
#
# Same computation as the PERMANOVA block of run_picrust2_factorial.R, which ran
# MetaCyc only and is kept for its ANCOM-BC2 record.
#
# Writes tables/PICRUSt2_permanova.tsv.
#
# Usage, in containers/ancombc2-env.sif (vegan):
#   Rscript picrust2_permanova.R

suppressPackageStartupMessages(library(vegan))
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
PIC      <- file.path(OUT_DIR, "picrust2", "pipeline")
NPERM <- 9999L
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")

meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "", na.strings = character(0))
m <- data.frame(Diet = factor(meta$Diet, c("Control", "SCFAs")), Radiation = factor(meta$Radiation, c("None", "IR")),
                row.names = rownames(meta))
out <- do.call(rbind, lapply(names(FILES), function(lv) {
  ab <- as.matrix(read.delim(gzfile(file.path(PIC, FILES[[lv]])), row.names = 1, check.names = FALSE)[, rownames(meta)])
  d <- vegdist(t(sweep(ab, 2, colSums(ab), "/")), method = "bray")
  set.seed(42); a1 <- adonis2(d ~ Diet * Radiation, data = m, permutations = NPERM, by = "terms")
  set.seed(42); a2 <- adonis2(d ~ Radiation * Diet, data = m, permutations = NPERM, by = "terms")
  ad <- rbind(a2["Diet", ], a1["Radiation", ], a1["Diet:Radiation", ])
  set.seed(42); pt <- permutest(betadisper(d, interaction(m$Diet, m$Radiation, drop = TRUE)), permutations = NPERM)
  data.frame(level = lv, n_features = nrow(ab),
             analysis = c(rep("permanova_bray_curtis", 3), "permdisp_bray_curtis"),
             term = c("Diet", "Radiation", "Diet:Radiation", "group"),
             statistic = c(ad$F, pt$tab[1, "F"]), effect = c(ad$R2, NA), p = c(ad$`Pr(>F)`, pt$tab[1, "Pr(>F)"]))
}))
write.table(out, file.path(OUT_DIR, "tables", "PICRUSt2_permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(out, row.names = FALSE)
