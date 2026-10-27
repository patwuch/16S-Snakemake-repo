#!/usr/bin/env Rscript
# ggpicrust2 figures for the PICRUSt2 section of tisha251224's factorial report.
#   picrust2_<level>_pca      - ggpicrust2::pathway_pca per level (MetaCyc, KO, EC)
#   picrust2_<level>_heatmap  - ggpicrust2::pathway_heatmap of the TOP_N features with the
#                               smallest abundance p in each MaAsLin3 row, pooled
#   picrust2_<level>_volcano  - ggpicrust2::pathway_volcano per MaAsLin3 row (abundance
#                               model): coefficient against raw p, robust calls named
# Results come from run_picrust2_maaslin3.R (2026-10-02; ggpicrust2's LinDA was dropped after
# its mode check). Group colours come from the study metadata, as in every other figure.
# Adapted from reza250909/additional_comparison_scripts/plot_ggpicrust2.R.
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript plot_ggpicrust2.R pca|results

suppressPackageStartupMessages({ library(ggpicrust2); library(ggplot2); library(patchwork) })
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
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

# pathway_pca's 95% stat_ellipse needs at least four samples per group, so it skips every group of
# two or three here. Replace it with the beta-diversity PCoA's ellipses (plot_diversity_factorial.py):
# one standard deviation along each principal axis of the group's points, for every group with at
# least three samples. Visual grouping only, not a confidence region.
sd_ellipses <- function(p) {
  p$layers <- Filter(function(l) !inherits(l$stat, "StatEllipse"), p$layers)
  # pathway_pca matches colours to groups by position, in the order the groups first appear in its
  # data, not in the order given; set both scales by name and the legend in design order.
  p$data$Group <- factor(p$data$Group, levels = levels(md$Group))
  t <- seq(0, 2 * pi, length.out = 100)
  el <- do.call(rbind, lapply(split(p$data, p$data$Group), function(d) {
    if (nrow(d) < 3) return(NULL)
    e <- eigen(cov(d[, c("PC1", "PC2")]), symmetric = TRUE)
    xy <- t(e$vectors %*% (sqrt(pmax(e$values, 0)) * rbind(cos(t), sin(t))))
    data.frame(PC1 = mean(d$PC1) + xy[, 1], PC2 = mean(d$PC2) + xy[, 2], Group = d$Group[1])
  }))
  p$layers <- c(geom_polygon(data = el, aes(PC1, PC2, color = Group, fill = Group), alpha = 0.10,
                             linetype = "dashed", linewidth = 0.5, show.legend = FALSE), p$layers)
  p + scale_colour_manual(values = cols, breaks = levels(md$Group)) + scale_fill_manual(values = cols, guide = "none")
}

if (WHAT == "pca") {
  # PCA per level (pathway_pca centres and scales each feature, so features with no
  # variance across the 19 samples are dropped first).
  for (lv in names(LEVEL_NAME)) {
    ab <- read_ab(lv)
    ab <- ab[apply(ab, 1, sd) > 0, , drop = FALSE]
    pca <- pathway_pca(abundance = ab, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                       show_marginal = FALSE)
    pca <- sd_ellipses(pca) +
      ggtitle(sprintf("Predicted %s: PCA of %s features", LEVEL_NAME[[lv]], format(nrow(ab), big.mark = ",")))
    save(pca, paste0("picrust2_", lv, "_pca"), 7, 6)
  }
  quit(save = "no")
}

# MaAsLin3 rows as the report names them: label, the groups it compares (first, second) for colour.
ROWS <- data.frame(
  row = c(rep(c("Inulin vs Cellulose", "RS vs Cellulose", "RS vs Inulin"), 2),
          "IR within Cellulose", "IR within Inulin", "IR within RS",
          "Interaction: Inulin vs Cellulose", "Interaction: RS vs Cellulose", "Interaction: RS vs Inulin"),
  stratum = c(rep("Control", 3), rep("IR", 3), rep("", 6)),
  label = c("Inulin vs Cellulose", "RS vs Cellulose", "RS vs Inulin",
            "Inulin+IR vs Cellulose+IR", "RS+IR vs Cellulose+IR", "RS+IR vs Inulin+IR",
            "Cellulose+IR vs Cellulose", "Inulin+IR vs Inulin", "RS+IR vs RS",
            "Interaction: Inulin vs Cellulose", "Interaction: RS vs Cellulose", "Interaction: RS vs Inulin"),
  first = c("Inulin", "RS", "RS", "Inulin+IR", "RS+IR", "RS+IR", "Cellulose+IR", "Inulin+IR", "RS+IR", NA, NA, NA),
  second = c("Cellulose", "Cellulose", "Inulin", "Cellulose+IR", "Cellulose+IR", "Inulin+IR", "Cellulose", "Inulin", "RS", NA, NA, NA),
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
                    caption = "Dashed line at raw p = 0.05; colour marks raw p < 0.05; names mark robust calls (q < 0.05 in all 19 leave-one-out refits)")
  set.seed(42)   # ggrepel places the labels when the plot is drawn
  save(vol, paste0("picrust2_", lv, "_volcano"), 16, 19)
}
