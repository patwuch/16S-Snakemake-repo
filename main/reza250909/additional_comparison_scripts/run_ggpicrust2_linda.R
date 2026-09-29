#!/usr/bin/env Rscript
# Per-feature tests on PICRUSt2's predictions with ggpicrust2 (CRAN 2.5.17):
# pathway_daa(daa_method = "LinDA") at three levels, each annotated by
# ggpicrust2's pathway_annotation:
#   MetaCyc  pathways_out/path_abun_unstrat.tsv.gz          (pathways)
#   KO       KO_metagenome_out/pred_metagenome_unstrat.tsv.gz (KEGG orthologs)
#   EC       EC_metagenome_out/pred_metagenome_unstrat.tsv.gz (enzymes)
#
# ggpicrust2 fits `~ Group` over the four Diet x Radiation groups and compares
# a reference group with each other group; it has no formula, covariate or
# interaction option (the interaction is left to the community PERMANOVA in
# picrust2_permanova.R). The report's four contrasts come from three calls:
#   reference Control   -> SCFAs vs Control, IR vs Control
#   reference IR        -> SCFAs+IR vs IR
#   reference SCFAs     -> SCFAs+IR vs SCFAs
# LinDA runs with ggpicrust2's settings (count input, pseudo-count 0.5 for zeros,
# adaptive zero handling, no prevalence filter) on PICRUSt2's unrounded output.
#
# q: ggpicrust2 applies one BH pass across all rows of a call, i.e. across the
# contrasts that call happens to produce, so its q would change with the
# reference chosen. Here q is BH across features within each contrast and
# level, as LinDA itself reports it and as the taxon analysis does;
# ggpicrust2's pooled value is kept as q_ggpicrust2_pooled.
#
# Descriptions: ggpicrust2's reference lists lack some PICRUSt2 2.6 IDs; those
# are filled from PICRUSt2's own add_descriptions.py output.
#
# Writes tables/PICRUSt2_<level>_linda[_minpresN].tsv (one row per feature x contrast).
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript run_ggpicrust2_linda.R [--min-presence N] [--exclude ID[,ID]] [MetaCyc KO EC]

suppressPackageStartupMessages({ library(ggpicrust2); library(readr); library(dplyr) })
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR   <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
PIC       <- file.path(OUT_DIR, "picrust2")
message("[ggpicrust2] ggpicrust2 ", packageVersion("ggpicrust2"), ", MicrobiomeStat ", packageVersion("MicrobiomeStat"))

LEVELS <- list(
  MetaCyc = list(file = "pipeline/pathways_out/path_abun_unstrat.tsv.gz",
                 described = "described/pathway_abun_unstrat_described.tsv.gz"),
  KO      = list(file = "pipeline/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 described = "described/KO_metagenome_unstrat_described.tsv.gz"),
  EC      = list(file = "pipeline/EC_metagenome_out/pred_metagenome_unstrat.tsv.gz",
                 described = "described/EC_metagenome_unstrat_described.tsv.gz"))
args <- commandArgs(trailingOnly = TRUE)
# --min-presence N: refit on the features predicted in at least N of the 32 mice
# (ggpicrust2's LinDA wrapper has no prevalence filter); outputs get _minpresN.
i <- match("--min-presence", args)
MIN_PRESENCE <- if (is.na(i)) 0L else as.integer(args[i + 1L])
if (!is.na(i)) args <- args[-c(i, i + 1L)]
suffix <- if (MIN_PRESENCE > 0) paste0("_minpres", MIN_PRESENCE) else ""
# --exclude ID[,ID]: refit without these mice (a one-mouse sensitivity check); outputs get _exclID.
j <- match("--exclude", args)
EXCLUDE <- if (is.na(j)) character(0) else strsplit(args[j + 1L], ",")[[1]]
if (!is.na(j)) args <- args[-c(j, j + 1L)]
if (length(EXCLUDE)) suffix <- paste0(suffix, "_excl", paste(EXCLUDE, collapse = "-"))
which_levels <- if (length(args)) args else names(LEVELS)

meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE, na.strings = character(0))
meta <- data.frame(sample = meta$sampleid, Group = meta$Group, stringsAsFactors = FALSE)
if (length(EXCLUDE)) {
  if (!all(EXCLUDE %in% meta$sample)) stop("--exclude: unknown sample(s)")
  meta <- meta[!meta$sample %in% EXCLUDE, ]
}
want <- list(Control = c("SCFAs", "IR"), IR = "SCFAs+IR", SCFAs = "SCFAs+IR")
decode <- function(t) vapply(t, function(x) if (is.na(x)) NA_character_ else
  xml2::xml_text(xml2::read_html(paste0("<p>", x, "</p>"))), "", USE.NAMES = FALSE)

for (lv in which_levels) {
  spec <- LEVELS[[lv]]
  ab <- as.data.frame(read_delim(file.path(PIC, spec$file), delim = "\t", show_col_types = FALSE))
  rownames(ab) <- ab[[1]]
  ab <- ab[, meta$sample]
  if (MIN_PRESENCE > 0) ab <- ab[rowSums(ab > 0) >= MIN_PRESENCE, , drop = FALSE]
  message(sprintf("[ggpicrust2] %s: %d features x %d samples (min presence %d)", lv, nrow(ab), ncol(ab), MIN_PRESENCE))

  res <- do.call(rbind, lapply(names(want), function(ref) {
    d <- pathway_daa(abundance = ab, metadata = meta, group = "Group", daa_method = "LinDA",
                     reference = ref, p_adjust_method = "BH")
    d <- d[d$group1 == ref & d$group2 %in% want[[ref]], , drop = FALSE]
    d$q_ggpicrust2_pooled <- d$p_adjust
    d
  }))
  res$contrast <- paste(res$group2, "vs", res$group1)
  res <- res %>% group_by(contrast) %>% mutate(q = p.adjust(p_values, method = "BH")) %>% ungroup() %>% as.data.frame()

  ann <- pathway_annotation(pathway = lv, daa_results_df = res[!duplicated(res$feature), ], ko_to_kegg = FALSE)
  res$description <- ann$description[match(res$feature, ann$feature)]
  p2 <- read.delim(gzfile(file.path(PIC, spec$described)), check.names = FALSE, quote = "")
  miss <- is.na(res$description) | res$description %in% c("", "not_found")
  message(sprintf("[ggpicrust2] %s: %d of %d features lack a ggpicrust2 description; filled from PICRUSt2",
                  lv, length(unique(res$feature[miss])), length(unique(res$feature))))
  res$description[miss] <- p2$description[match(res$feature[miss], p2[[1]])]
  res$description <- decode(res$description)
  n_pred <- rowSums(ab > 0)
  res$n_predicted <- unname(n_pred[res$feature])
  res$level <- lv

  out <- res[, c("level", "contrast", "feature", "description", "n_predicted", "log2_fold_change", "p_values", "q",
                 "q_ggpicrust2_pooled", "method", "group1", "group2")]
  names(out)[names(out) == "log2_fold_change"] <- "log2fc"
  names(out)[names(out) == "p_values"] <- "p"
  out <- out[order(out$contrast, out$p), ]
  write.table(out, file.path(TAB_DIR, paste0("PICRUSt2_", lv, "_linda", suffix, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("[ggpicrust2] %s per contrast: q < 0.05 | raw p < 0.05 | features tested", lv))
  for (k in unique(out$contrast)) {
    x <- out[out$contrast == k, ]
    message(sprintf("    %-20s %5d | %5d | %5d", k, sum(x$q < 0.05, na.rm = TRUE), sum(x$p < 0.05, na.rm = TRUE), nrow(x)))
  }
}
