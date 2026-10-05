#!/usr/bin/env Rscript
# Per-feature LinDA tests on PICRUSt2's predictions (picrust2/, run_picrust2.sh)
# for both ahmad260430 reports, at three levels:
#   MetaCyc  pathways_out/path_abun_unstrat.tsv.gz            (pathways)
#   KO       KO_metagenome_out/pred_metagenome_unstrat.tsv.gz (KEGG orthologs)
#   EC       EC_metagenome_out/pred_metagenome_unstrat.tsv.gz (enzymes)
#
# factorial   (24 Day 17 mice, Treatment x Radiation, 4 per group)
#   ggpicrust2 (CRAN 2.5.17) pathway_daa(daa_method = "LinDA"), which fits
#   `~ Group` over the six groups and compares a reference group with each other
#   group; it has no formula, covariate or interaction option (the interaction is
#   left to the community PERMANOVA in picrust2_permanova.R). Each contrast is
#   read from the call whose reference is its second-named group.
# day17-day25 (the IFA and IFA+IR arms, the same 4 mice per arm at both days)
#   ggpicrust2 cannot pair samples, so MicrobiomeStat::linda, the function
#   ggpicrust2 wraps, is called directly with ggpicrust2's settings and a random
#   intercept per mouse: `~ Timepoint + (1 | MouseID)` within each arm, and
#   `~ Arm * Timepoint + (1 | MouseID)` on both arms for whether the shift differs.
#
# LinDA settings throughout are ggpicrust2's: count input, pseudo-count 0.5 for
# zeros, adaptive zero handling, no prevalence filter, features absent from every
# sample of a fit dropped. q is Benjamini-Hochberg across features within each
# contrast and level (ggpicrust2 would pool the contrasts of one call).
# Descriptions come from ggpicrust2's pathway_annotation, completed from
# PICRUSt2's add_descriptions.py maps where ggpicrust2's lists lack an ID.
#
# Writes <report dir>/tables/PICRUSt2_<level>_linda[_minpresN][_exclID].tsv,
# one row per feature x contrast.
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript run_picrust2_linda.R factorial|day17-day25 [--min-presence N] [--exclude ID[,ID]] [MetaCyc KO EC]

suppressPackageStartupMessages({ library(ggpicrust2); library(readr); library(dplyr) })
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
PIC       <- file.path(STUDY_DIR, "picrust2")
message("[linda] ggpicrust2 ", packageVersion("ggpicrust2"), ", MicrobiomeStat ", packageVersion("MicrobiomeStat"))

LEVELS <- list(
  MetaCyc = list(file = "pipeline/pathways_out/path_abun_unstrat.tsv.gz",
                 described = "described/pathway_abun_unstrat_described.tsv.gz"),
  KO      = list(file = "pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 described = "described/KO_metagenome_unstrat_described.tsv.gz"),
  EC      = list(file = "pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 described = "described/EC_metagenome_unstrat_described.tsv.gz"))

study_meta <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), check.names = FALSE,
                         stringsAsFactors = FALSE, na.strings = character(0))

# ggpicrust2 with one reference per call; `want` maps each reference to the groups read against it.
fit_ggpicrust2 <- function(ab, meta, want) {
  do.call(rbind, lapply(names(want), function(ref) {
    d <- pathway_daa(abundance = ab, metadata = meta, group = "Group", daa_method = "LinDA",
                     reference = ref, p_adjust_method = "BH")
    d <- d[d$group1 == ref & d$group2 %in% want[[ref]], , drop = FALSE]
    data.frame(contrast = paste(d$group2, "vs", d$group1), feature = d$feature,
               log2fc = d$log2_fold_change, p = d$p_values, method = "LinDA (ggpicrust2)")
  }))
}

# MicrobiomeStat::linda with ggpicrust2's settings; `terms` maps each contrast name to its coefficient.
fit_linda <- function(ab, meta, formula, terms) {
  ab <- ab[rowSums(ab) > 0, , drop = FALSE]
  o <- MicrobiomeStat::linda(feature.dat = as.matrix(ab), meta.dat = meta, formula = formula,
                             feature.dat.type = "count", prev.filter = 0, mean.abund.filter = 0,
                             adaptive = TRUE, zero.handling = "pseudo-count", pseudo.cnt = 0.5,
                             p.adj.method = "BH", alpha = 0.05, n.cores = 1, verbose = FALSE)$output
  do.call(rbind, lapply(names(terms), function(k) {
    x <- o[[terms[[k]]]]
    data.frame(contrast = k, feature = rownames(x), log2fc = x$log2FoldChange, p = x$pvalue,
               method = paste("LinDA, mixed model", formula))
  }))
}

REPORTS <- list(
  factorial = list(
    dir = "comparison-Treatment-x-Radiation-factorial",
    meta = function() {
      m <- study_meta[study_meta$Timepoint == "Day17", ]
      data.frame(sample = m$sampleid, Group = m$Group, stringsAsFactors = FALSE)
    },
    fit = function(ab, meta) fit_ggpicrust2(ab, meta, list(
      Control = c("IFA", "Inosine", "IR"), Inosine = c("IFA", "Inosine+IR"), IR = c("IFA+IR", "Inosine+IR"),
      "Inosine+IR" = "IFA+IR", IFA = "IFA+IR"))),
  "day17-day25" = list(
    dir = "comparison-day17-day25",
    meta = function() {
      m <- read.delim(file.path(STUDY_DIR, "comparison-day17-day25", "metadata_joint_cells.tsv"),
                      check.names = FALSE, stringsAsFactors = FALSE, na.strings = character(0))
      data.frame(sample = m$sampleid, MouseID = m$MouseID,
                 Arm = factor(sub(" Day.*", "", m$Group), c("IFA", "IFA+IR")),
                 Timepoint = factor(sub(".* ", "", m$Group), c("Day17", "Day25")),
                 row.names = m$sampleid, stringsAsFactors = FALSE)
    },
    fit = function(ab, meta) {
      arm <- function(a) {
        s <- meta$Arm == a
        fit_linda(ab[, s, drop = FALSE], meta[s, ], "~ Timepoint + (1 | MouseID)",
                  setNames(list("TimepointDay25"), paste(a, "Day25 vs", a, "Day17")))
      }
      rbind(arm("IFA"), arm("IFA+IR"),
            fit_linda(ab, meta, "~ Arm * Timepoint + (1 | MouseID)",
                      list("Shift, IFA+IR vs IFA" = "ArmIFA+IR:TimepointDay25")))
    }))

args <- commandArgs(trailingOnly = TRUE)
if (!length(args) || !args[1] %in% names(REPORTS)) stop("first argument: ", paste(names(REPORTS), collapse = " or "))
REPORT <- args[1]; report <- REPORTS[[REPORT]]; args <- args[-1]
TAB_DIR <- file.path(STUDY_DIR, report$dir, "tables")
# --min-presence N: refit on the features predicted in at least N of the report's samples.
i <- match("--min-presence", args)
MIN_PRESENCE <- if (is.na(i)) 0L else as.integer(args[i + 1L])
if (!is.na(i)) args <- args[-c(i, i + 1L)]
suffix <- if (MIN_PRESENCE > 0) paste0("_minpres", MIN_PRESENCE) else ""
# --exclude ID[,ID]: refit without these samples (a one-mouse sensitivity check).
j <- match("--exclude", args)
EXCLUDE <- if (is.na(j)) character(0) else strsplit(args[j + 1L], ",")[[1]]
if (!is.na(j)) args <- args[-c(j, j + 1L)]
if (length(EXCLUDE)) suffix <- paste0(suffix, "_excl", paste(EXCLUDE, collapse = "-"))
which_levels <- if (length(args)) args else names(LEVELS)

meta <- report$meta()
if (length(EXCLUDE)) {
  if (!all(EXCLUDE %in% meta$sample)) stop("--exclude: unknown sample(s)")
  meta <- meta[!meta$sample %in% EXCLUDE, ]
}
decode <- function(t) vapply(t, function(x) if (is.na(x)) NA_character_ else
  xml2::xml_text(xml2::read_html(paste0("<p>", x, "</p>"))), "", USE.NAMES = FALSE)

for (lv in which_levels) {
  spec <- LEVELS[[lv]]
  ab <- as.data.frame(read_delim(file.path(PIC, spec$file), delim = "\t", show_col_types = FALSE))
  rownames(ab) <- ab[[1]]
  ab <- ab[, meta$sample]
  if (MIN_PRESENCE > 0) ab <- ab[rowSums(ab > 0) >= MIN_PRESENCE, , drop = FALSE]
  message(sprintf("[linda] %s %s: %d features x %d samples (min presence %d)",
                  REPORT, lv, nrow(ab), ncol(ab), MIN_PRESENCE))

  res <- report$fit(ab, meta)
  res <- res %>% group_by(contrast) %>% mutate(q = p.adjust(p, method = "BH")) %>% ungroup() %>% as.data.frame()

  ann <- pathway_annotation(pathway = lv, ko_to_kegg = FALSE,
                            daa_results_df = data.frame(feature = unique(res$feature), p_adjust = 1))
  res$description <- ann$description[match(res$feature, ann$feature)]
  p2 <- read.delim(gzfile(file.path(PIC, spec$described)), check.names = FALSE, quote = "")
  miss <- is.na(res$description) | res$description %in% c("", "not_found")
  message(sprintf("[linda] %s: %d of %d features lack a ggpicrust2 description; filled from PICRUSt2",
                  lv, length(unique(res$feature[miss])), length(unique(res$feature))))
  res$description[miss] <- p2$description[match(res$feature[miss], p2[[1]])]
  res$description <- decode(res$description)
  res$n_predicted <- unname(rowSums(ab > 0)[res$feature])
  res$level <- lv

  out <- res[, c("level", "contrast", "feature", "description", "n_predicted", "log2fc", "p", "q", "method")]
  out <- out[order(out$contrast, out$p), ]
  write.table(out, file.path(TAB_DIR, paste0("PICRUSt2_", lv, "_linda", suffix, ".tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("[linda] %s per contrast: q < 0.05 | raw p < 0.05 | features tested", lv))
  for (k in unique(out$contrast)) {
    x <- out[out$contrast == k, ]
    message(sprintf("    %-28s %5d | %5d | %5d", k, sum(x$q < 0.05, na.rm = TRUE), sum(x$p < 0.05, na.rm = TRUE), nrow(x)))
  }
}
