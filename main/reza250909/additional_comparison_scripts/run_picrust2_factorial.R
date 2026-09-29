#!/usr/bin/env Rscript
# Diet x Radiation tests on PICRUSt2's predicted MetaCyc pathway abundances
# (run_picrust2.sh), with the same design, contrasts and ANCOM-BC2 settings as
# the taxon analysis (run_ancombc2_factorial.R).
#
#  * PERMANOVA on Bray-Curtis distances between samples' relative pathway
#    profiles, Diet * Radiation, Type II (each main effect entering last),
#    9,999 permutations, seed 42; PERMDISP across the four groups.
#  * ANCOM-BC2 per pathway, Diet * Radiation, group = NULL, BH across pathways,
#    pseudo-count sensitivity check. Pathway abundances are rounded to integers
#    first (PICRUSt2 reports predicted counts with decimals). A pathway is tested
#    only if it is predicted in at least one sample of every group and in at
#    least 4 of the 32, the taxon rule. Each report row is one coefficient,
#    read from the fit whose reference matches it:
#      SCFAs vs Control        DietSCFAs,   reference Control diet, unirradiated
#      SCFAs+IR vs IR          DietSCFAs,   reference Control diet, irradiated
#      IR vs Control           RadiationIR, reference Control diet
#      SCFAs+IR vs SCFAs       RadiationIR, reference SCFAs diet
#      interaction             DietSCFAs:RadiationIR (identical in every fit)
#
# Writes to comparison-Diet-x-Radiation-factorial/tables/:
#   PICRUSt2_pathways_permanova.tsv, PICRUSt2_pathways_ancombc2.tsv (one row per
#   pathway x contrast), PICRUSt2_pathways_ancombc2_<fit>.tsv (raw fits).
#
# Usage, in containers/ancombc2-env.sif:
#   Rscript run_picrust2_factorial.R

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
PIC_DIR   <- file.path(OUT_DIR, "picrust2")
LOG_PREFIX <- "[picrust2]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
suppressPackageStartupMessages({ library(ANCOMBC); library(vegan) })
require_ancombc("2.14.0")

ALPHA <- 0.05; NPERM <- 9999L; MIN_PRESENCE <- 4L
DIET_LEVELS <- c("Control", "SCFAs"); RADIATION_LEVELS <- c("None", "IR")

meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE, quote = "", comment.char = "")
ab <- read.delim(gzfile(file.path(PIC_DIR, "pipeline", "pathways_out", "path_abun_unstrat.tsv.gz")),
                 row.names = 1, check.names = FALSE)
desc <- read.delim(gzfile(file.path(PIC_DIR, "described", "pathway_abun_unstrat_described.tsv.gz")),
                   check.names = FALSE, quote = "")
desc <- setNames(desc$description, desc$pathway)
missing <- setdiff(rownames(meta), colnames(ab))
if (length(missing)) stop("samples missing from the pathway table: ", paste(missing, collapse = ", "))
ab <- as.matrix(ab[, rownames(meta)])
.say(sprintf("%d MetaCyc pathways x %d samples", nrow(ab), ncol(ab)))

# ---- PERMANOVA on relative pathway profiles ----
rel <- sweep(ab, 2, colSums(ab), "/")
d <- vegdist(t(rel), method = "bray")
m <- data.frame(Diet = factor(meta$Diet, DIET_LEVELS), Radiation = factor(meta$Radiation, RADIATION_LEVELS),
                row.names = rownames(meta))
set.seed(42); a1 <- adonis2(d ~ Diet * Radiation, data = m, permutations = NPERM, by = "terms")
set.seed(42); a2 <- adonis2(d ~ Radiation * Diet, data = m, permutations = NPERM, by = "terms")
ad <- rbind(a2["Diet", ], a1["Radiation", ], a1["Diet:Radiation", ])
bd <- betadisper(d, interaction(m$Diet, m$Radiation, drop = TRUE))
set.seed(42); pt <- permutest(bd, permutations = NPERM)
perm <- data.frame(analysis = c(rep("permanova_bray_curtis", 3), "permdisp_bray_curtis"),
                   term = c("Diet", "Radiation", "Diet:Radiation", "group"),
                   statistic = c(ad$F, pt$tab[1, "F"]), effect = c(ad$R2, NA), p = c(ad$`Pr(>F)`, pt$tab[1, "Pr(>F)"]))
write.table(perm, file.path(TAB_DIR, "PICRUSt2_pathways_permanova.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(perm, row.names = FALSE)

# ---- ANCOM-BC2 per pathway ----
# Main table: the taxon rule (predicted in every group and in >= 4 of 32 mice).
# Sweep: the same fits among pathways predicted in >= 8, 12 and 16 mice, as the
# taxon table's prevalence column does. ANCOM-BC2 treats a zero as a missing
# log value, so a pathway predicted in a few mice is estimated from those few
# alone; the sweep shows which calls rest on them.
counts_all <- round(ab)
cell <- interaction(m$Diet, m$Radiation, drop = TRUE)
in_every <- sapply(levels(cell), function(k) rowSums(counts_all[, cell == k, drop = FALSE] > 0) > 0)
n_pred <- rowSums(counts_all > 0)

fit <- function(counts, diet_ref, rad_ref) {
  mm <- m
  mm$Diet <- relevel(mm$Diet, diet_ref); mm$Radiation <- relevel(mm$Radiation, rad_ref)
  r <- withCallingHandlers(
    ancombc2(data = counts, taxa_are_rows = TRUE, meta_data = mm, fix_formula = "Diet * Radiation",
             group = NULL, prv_cut = 0.05, lib_cut = 0, p_adj_method = "BH", alpha = ALPHA,
             pseudo_sens = TRUE, struc_zero = FALSE, neg_lb = FALSE, global = FALSE, pairwise = FALSE,
             dunnet = FALSE, trend = FALSE, verbose = FALSE),
    warning = function(w) if (grepl("multiple of", conditionMessage(w))) stop(conditionMessage(w)))
  as.data.frame(r$res)
}
rows <- list(c("SCFAs vs Control", "default", "DietSCFAs"),
             c("SCFAs+IR vs IR", "radrefIR", "DietSCFAs"),
             c("IR vs Control", "default", "RadiationIR"),
             c("SCFAs+IR vs SCFAs", "refSCFAs", "RadiationIR"),
             c("interaction", "default", "DietSCFAs:RadiationIR"))

run_at <- function(min_presence) {
  keep <- rowSums(!in_every) == 0 & n_pred >= min_presence
  .say(sprintf("min presence %d: %d of %d pathways tested", min_presence, sum(keep), nrow(counts_all)))
  counts <- counts_all[keep, , drop = FALSE]
  fits <- list(default = fit(counts, "Control", "None"), radrefIR = fit(counts, "Control", "IR"),
               refSCFAs = fit(counts, "SCFAs", "None"))
  long <- do.call(rbind, lapply(rows, function(r) {
    f <- fits[[r[2]]]; co <- r[3]
    data.frame(contrast = r[1], pathway = f$taxon, description = unname(desc[f$taxon]),
               n_predicted = unname(n_pred[f$taxon]),
               log2fc = f[[paste0("lfc_", co)]] / log(2), se_log2 = f[[paste0("se_", co)]] / log(2),
               p = f[[paste0("p_", co)]], q = f[[paste0("q_", co)]],
               passed_ss = f[[paste0("passed_ss_", co)]] %in% c(TRUE, "TRUE"), check.names = FALSE)
  }))
  long$robust <- long$q < ALPHA & long$passed_ss
  list(fits = fits, long = long, n = sum(keep))
}

main <- run_at(MIN_PRESENCE)
for (k in names(main$fits))
  write.table(main$fits[[k]], file.path(TAB_DIR, paste0("PICRUSt2_pathways_ancombc2_", k, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
write.table(main$long, file.path(TAB_DIR, "PICRUSt2_pathways_ancombc2.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

sweep_long <- list()
sweep <- do.call(rbind, lapply(c(MIN_PRESENCE, 8L, 12L, 16L), function(mp) {
  res <- if (mp == MIN_PRESENCE) main else run_at(mp)
  sweep_long[[as.character(mp)]] <<- cbind(min_presence = mp, res$long)
  do.call(rbind, lapply(rows, function(r) {
    x <- res$long[res$long$contrast == r[1], ]
    data.frame(min_presence = mp, n_tested = res$n, contrast = r[1], robust = sum(x$robust),
               q_below = sum(x$q < ALPHA), p_below = sum(x$p < ALPHA))
  }))
}))
write.table(sweep, file.path(TAB_DIR, "PICRUSt2_pathways_ancombc2_prevalence_sweep.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(sweep, row.names = FALSE)
write.table(do.call(rbind, sweep_long), file.path(TAB_DIR, "PICRUSt2_pathways_ancombc2_prevalence_sweep_long.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
