#!/usr/bin/env Rscript
# ggpicrust2 figures for the PICRUSt2 section of reza250909's factorial report.
#   picrust2_<level>_pca      - ggpicrust2::pathway_pca per level (MetaCyc, KO, EC)
#   picrust2_<level>_heatmap  - ggpicrust2::pathway_heatmap of the TOP_N features with the
#                               smallest abundance p in each MaAsLin3 row, pooled
#   picrust2_<level>_volcano  - ggpicrust2::pathway_volcano per MaAsLin3 row (abundance
#                               model): coefficient against raw p, robust calls named
# Results come from run_picrust2_maaslin3.R (2026-10-04; ggpicrust2's LinDA was dropped after
# its mode check). The LinDA version is kept in the scratchpad as reza_plot_ggpicrust2_linda.R. Group colours come from the study metadata, as in every other figure.
# Ported back from tisha251224/additional_comparison_scripts/plot_ggpicrust2.R.
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript plot_ggpicrust2.R pca|results

suppressPackageStartupMessages({ library(ggpicrust2); library(ggplot2); library(patchwork) })
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
FIG_DIR  <- file.path(OUT_DIR, "figures")
WHAT <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(WHAT) || !WHAT %in% c("pca", "results")) stop("argument: pca or results")

meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE, na.strings = character(0))
spec <- unique(meta[order(meta$PrimaryOrder), c("Group", "PrimaryColor")])
cols <- setNames(spec$PrimaryColor, spec$Group)
md <- data.frame(sample = meta$sampleid, Group = factor(meta$Group, levels = spec$Group))
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
read_ab <- function(lv) as.matrix(read.delim(gzfile(file.path(OUT_DIR, "picrust2", "pipeline", FILES[[lv]])),
                                              row.names = 1, check.names = FALSE)[, md$sample])
save <- function(p, stem, w, h) {
  for (ext in c("svg", "png")) ggsave(file.path(FIG_DIR, paste0(stem, ".", ext)), p, width = w, height = h, dpi = 300)
  message("  wrote ", stem)
}
LEVEL_NAME <- c(MetaCyc = "MetaCyc pathways", KO = "KEGG orthologs (KO)", EC = "enzymes (EC)")

if (WHAT == "pca") {
  # PCA per level (pathway_pca centres and scales each feature, so features with no
  # variance across the 32 mice are dropped first).
  for (lv in names(LEVEL_NAME)) {
    ab <- read_ab(lv)
    ab <- ab[apply(ab, 1, sd) > 0, , drop = FALSE]
    pca <- pathway_pca(abundance = ab, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                       show_marginal = FALSE) +
      ggtitle(sprintf("Predicted %s: PCA of %s features", LEVEL_NAME[[lv]], format(nrow(ab), big.mark = ",")))
    save(pca, paste0("picrust2_", lv, "_pca"), 7, 6)
  }
  quit(save = "no")
}

# MaAsLin3 rows as the report names them: label, the groups it compares (first, second) for colour.
ROWS <- data.frame(
  row = c("SCFAs vs Control", "SCFAs+IR vs IR", "IR vs Control", "SCFAs+IR vs SCFAs", "Interaction: SCFAs vs Control"),
  stratum = "",
  label = c("SCFAs vs Control", "SCFAs+IR vs IR", "IR vs Control", "SCFAs+IR vs SCFAs", "Interaction: SCFAs vs Control"),
  first = c("SCFAs", "SCFAs+IR", "IR", "SCFAs+IR", NA),
  second = c("Control", "IR", "Control", "SCFAs", NA),
  stringsAsFactors = FALSE)
ROWS$tag <- LETTERS[seq_len(nrow(ROWS))]
TOP_N <- 5L
read_rows <- function(lv) {
  d <- read.delim(file.path(OUT_DIR, "tables", paste0("PICRUSt2_", lv, "_maaslin3_rows.tsv")),
                  check.names = FALSE, quote = "", stringsAsFactors = FALSE, na.strings = c("NA", ""))
  d <- d[d$model == "abundance" & !is.na(d$p), ]
  d$stratum[is.na(d$stratum)] <- ""   # the irradiation and interaction rows have no stratum
  # PICRUSt2's descriptions carry HTML entities ("L-N&delta;-acetylornithine").
  d$description <- vapply(d$description, function(x) if (is.na(x)) NA_character_ else
    xml2::xml_text(xml2::read_html(paste0("<p>", x, "</p>"))), "", USE.NAMES = FALSE)
  d$label <- ROWS$label[match(paste(d$row, d$stratum), paste(ROWS$row, ROWS$stratum))]
  stopifnot(!anyNA(d$label))
  d
}
short_label <- function(lv, id, desc) {
  desc[is.na(desc)] <- "not_found"
  d <- if (lv == "KO") sub("^([^;]*); *(.*?)( \\[EC:.*\\])?$", "\\1: \\2", desc) else desc
  if (lv == "KO") d <- ifelse(desc == "not_found", paste0(sub("^ko:", "", id), ": no name in the KO map"), d)
  if (lv == "EC") {
    moved <- grepl("^Transferred entry: ", d)
    d <- ifelse(moved, paste0("EC ", id, " (now ", sub("^Transferred entry: *", "", d), ")"), paste0(d, " (", id, ")"))
  }
  d <- ifelse(nchar(d) > 60, paste0(substr(d, 1, 59), "…"), d)
  ifelse(duplicated(d) | duplicated(d, fromLast = TRUE), paste0(d, " (", id, ")"), d)
}

sel_tab <- list()
for (lv in names(LEVEL_NAME)) {
  res <- read_rows(lv)
  pick <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
    d <- res[res$label == ROWS$label[k], ]; d <- d[order(d$p), ][seq_len(min(TOP_N, nrow(d))), ]
    data.frame(level = lv, row = ROWS$label[k], tag = ROWS$tag[k], id = d$id, description = d$description,
               coef = d$coef, p = d$p, q = d$q, robust = d$robust, stringsAsFactors = FALSE)
  }))
  sel_tab[[lv]] <- pick
  feats <- unique(pick$id)
  ab <- read_ab(lv)
  rel <- sweep(ab, 2, colSums(ab), "/") * 100
  h <- rel[feats, , drop = FALSE]
  tags <- vapply(feats, function(ft) paste(sort(unique(pick$tag[pick$id == ft])), collapse = ","), "")
  rownames(h) <- paste0(short_label(lv, feats, pick$description[match(feats, pick$id)]), "  [", tags, "]")
  hm <- pathway_heatmap(abundance = h, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                        font_size = 9, cluster_rows = TRUE, clustering_method = "average",
                        clustering_distance = "correlation", colorbar_title = "z-score") +
    # ggpicrust2 passes breaks = colorbar_breaks (NULL by default) to its fill scale, and
    # breaks = NULL removes the colour bar. Same colours, default breaks.
    suppressMessages(ggplot2::scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
                                                   name = "z-score")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(theme = ggplot2::theme(
      legend.key.width = grid::unit(0.5, "cm"), legend.key.height = grid::unit(5, "cm")))) +
    ggplot2::theme(legend.position = "right") +
    ggtitle(sprintf("Predicted %s: the %d strongest per row (%d features)", LEVEL_NAME[[lv]], TOP_N, length(feats)))
  save(hm, paste0("picrust2_", lv, "_heatmap"), 13, 0.24 * length(feats) + 3)
}
write.table(do.call(rbind, sel_tab), file.path(OUT_DIR, "tables", "PICRUSt2_heatmap_features.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

for (lv in names(LEVEL_NAME)) {
  res <- read_rows(lv)
  res$name <- if (lv == "KO") sub(";.*", "", res$description) else
    ifelse(nchar(res$description) > 45, paste0(substr(res$description, 1, 44), "…"), res$description)
  panels <- lapply(seq_len(nrow(ROWS)), function(k) {
    d <- res[res$label == ROWS$label[k], ]
    # Name the robust calls only (at most five, smallest p first).
    d$lab <- ""
    r <- which(d$robust %in% TRUE); r <- r[order(d$p[r])][seq_len(min(5, length(r)))]
    d$lab[r] <- d$name[r]
    up <- if (is.na(ROWS$first[k])) "#b2182b" else unname(cols[ROWS$first[k]])
    dn <- if (is.na(ROWS$second[k])) "#2166ac" else unname(cols[ROWS$second[k]])
    pathway_volcano(d, fc_col = "coef", p_col = "p", label_col = "lab", fc_threshold = 0,
                    p_threshold = 0.05, label_top_n = 5, point_size = if (lv == "MetaCyc") 1.4 else 0.8,
                    colors = c(Down = dn, "Not Significant" = "grey70", Up = up),
                    title = sprintf("%s [%s]: %d at q < 0.05, %d robust", ROWS$label[k], ROWS$tag[k],
                                    sum(d$call %in% TRUE), sum(d$robust %in% TRUE)),
                    x_lab = if (is.na(ROWS$first[k])) "interaction coefficient (log2)" else
                      sprintf("log2 change in share (right = higher in %s)", ROWS$first[k]),
                    y_lab = "-log10(raw p)") +
      # ggpicrust2 subtitles each panel with raw-p counts; the title gives the q-based counts.
      labs(subtitle = NULL) +
      theme(plot.title = element_text(size = 9), legend.position = "none")
  })
  vol <- wrap_plots(panels, ncol = 3) +
    plot_annotation(title = paste("Predicted", LEVEL_NAME[[lv]], "· MaAsLin3, change in share of the predicted total"),
                    caption = "Dashed line at raw p = 0.05; colour marks raw p < 0.05; names mark robust calls (q < 0.05 in all 32 leave-one-mouse-out refits)")
  set.seed(42)   # ggrepel places the labels when the plot is drawn
  save(vol, paste0("picrust2_", lv, "_volcano"), 16, 11)
}
