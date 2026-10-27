#!/usr/bin/env Rscript
# Figures for the PICRUSt2 sections of both ahmad260430 reports, from picrust2/ (run_picrust2.sh)
# and the MaAsLin3 tables (run_picrust2_maaslin3.R; LinDA was dropped after its mode check,
# 2026-10-04):
#   picrust2_nsti             - weighted NSTI per sample, by group
#   picrust2_<level>_pca      - ggpicrust2::pathway_pca per level (MetaCyc, KO, EC)
#   picrust2_<level>_heatmap  - ggpicrust2::pathway_heatmap of the TOP_N features with the smallest
#                               abundance p in each MaAsLin3 row, pooled
#   picrust2_<level>_volcano  - ggpicrust2::pathway_volcano per MaAsLin3 row, robust calls named
# Factorial groups take the study's colours. Day17/Day25 keeps the arm colours (IFA light blue, IFA+IR
# blue) for both days and marks the day by shape: Day 17 circles, Day 25 triangles (user, 2026-10-04).
# Heatmap column strips stay in the arm colours; their labels name the day.
# Adapted from tisha251224/additional_comparison_scripts/plot_ggpicrust2.R.
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript plot_picrust2.R factorial|day17-day25

suppressPackageStartupMessages({ library(ggpicrust2); library(ggplot2); library(patchwork) })
THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
PIC       <- file.path(STUDY_DIR, "picrust2", "pipeline")
REPORT <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(REPORT) || !REPORT %in% c("factorial", "day17-day25")) stop("argument: factorial or day17-day25")

if (REPORT == "factorial") {
  OUT <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
  m <- read.delim(file.path(OUT, "tables", "metadata_day17.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  spec <- unique(m[order(m$PrimaryOrder), c("Group", "PrimaryColor")])
  cols <- setNames(spec$PrimaryColor, spec$Group)
  md <- data.frame(sample = m$sampleid, Group = factor(m$Group, levels = spec$Group))
  NREF <- "24 mice"
  ROWS <- data.frame(
    row = c(rep(c("IFA vs Control", "Inosine vs Control", "IFA vs Inosine"), 2), "IR within Control", "IR within IFA",
            "IR within Inosine", "Interaction: IFA vs Control", "Interaction: Inosine vs Control", "Interaction: IFA vs Inosine"),
    stratum = c(rep("None", 3), rep("IR", 3), rep("", 6)),
    # Baseline on the left (user, 2026-10-23); the estimate is the right-hand group's change.
    label = c("Control vs IFA", "Control vs Inosine", "Inosine vs IFA", "IR vs IFA+IR", "IR vs Inosine+IR",
              "Inosine+IR vs IFA+IR", "Control vs IR", "IFA vs IFA+IR", "Inosine vs Inosine+IR",
              "Interaction: Control vs IFA", "Interaction: Control vs Inosine", "Interaction: Inosine vs IFA"),
    first = c("IFA", "Inosine", "IFA", "IFA+IR", "Inosine+IR", "IFA+IR", "IR", "IFA+IR", "Inosine+IR", NA, NA, NA),
    second = c("Control", "Control", "Inosine", "IR", "IR", "Inosine+IR", "Control", "IFA", "Inosine", NA, NA, NA),
    stringsAsFactors = FALSE)
  NCOL <- 3; VOL_H <- 19
  vcols <- cols
  ROWS$xname <- ROWS$first
} else {
  OUT <- file.path(STUDY_DIR, "comparison-day17-day25")
  m <- read.delim(file.path(OUT, "metadata_joint_cells.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
  lv_g <- c("IFA Day17", "IFA Day25", "IFA+IR Day17", "IFA+IR Day25")
  arm <- c(IFA = "#90BFF9", "IFA+IR" = "#0000FF")
  cols <- setNames(arm[c("IFA", "IFA", "IFA+IR", "IFA+IR")], lv_g)
  vcols <- c(cols, IFA_ARM = arm[["IFA"]], IFAIR_ARM = arm[["IFA+IR"]])
  DAY_SHAPE <- c("Day 17" = 21, "Day 25" = 24)
  md <- data.frame(sample = m$sampleid, Group = factor(m$Group, levels = lv_g))
  NREF <- "16 samples"
  ROWS <- data.frame(
    row = c("Day 17 to 25 in IFA", "Day 17 to 25 in IFA+IR", "IFA+IR vs IFA at Day 17", "IFA+IR vs IFA at Day 25",
            "Interaction (IFA+IR shift minus IFA shift)", "Average Day 17 to 25 shift", "Average IFA+IR vs IFA"),
    stratum = "",
    label = c("Day 17 to 25, IFA", "Day 17 to 25, IFA+IR", "IFA+IR vs IFA, Day 17", "IFA+IR vs IFA, Day 25",
              "Interaction (shift difference)", "Average Day 17 to 25 shift", "Average IFA+IR vs IFA"),
    first = c("IFA Day25", "IFA+IR Day25", "IFAIR_ARM", "IFAIR_ARM", NA, NA, NA),
    second = c("IFA Day17", "IFA+IR Day17", "IFA_ARM", "IFA_ARM", NA, NA, NA),
    xname = c("Day 25", "Day 25", "IFA+IR", "IFA+IR", NA, NA, NA),
    stringsAsFactors = FALSE)
  NCOL <- 3; VOL_H <- 14
}
ROWS$tag <- LETTERS[seq_len(nrow(ROWS))]
FIG_DIR <- file.path(OUT, "figures")
FILES <- c(MetaCyc = "pathways_out/path_abun_unstrat.tsv.gz",
           KO = "KO_metagenome_out/pred_metagenome_unstrat.tsv.gz",
           EC = "EC_metagenome_out/pred_metagenome_unstrat.tsv.gz")
read_ab <- function(lv) as.matrix(read.delim(gzfile(file.path(PIC, FILES[[lv]])), row.names = 1, check.names = FALSE)[, md$sample])
save <- function(p, stem, w, h) {
  for (ext in c("svg", "png")) ggsave(file.path(FIG_DIR, paste0(stem, ".", ext)), p, width = w, height = h, dpi = 300)
  message("  wrote ", stem)
}
LEVEL_NAME <- c(MetaCyc = "MetaCyc pathways", KO = "KEGG orthologs (KO)", EC = "enzymes (EC)")
TOP_N <- 5L
decode <- function(t) vapply(t, function(x) if (is.na(x)) NA_character_ else
  xml2::xml_text(xml2::read_html(paste0("<p>", x, "</p>"))), "", USE.NAMES = FALSE)

# NSTI per sample.
wn <- read.delim(gzfile(file.path(PIC, "KO_metagenome_out", "weighted_nsti.tsv.gz")), check.names = FALSE)
d <- merge(md, wn, by.x = "sample", by.y = 1)
stopifnot(nrow(d) == nrow(md))
d$Day <- ifelse(grepl("Day17", d$Group), "Day 17", "Day 25")
p <- ggplot(d, aes(Group, weighted_NSTI, fill = Group)) + geom_boxplot(outlier.shape = NA, width = 0.6) +
  (if (REPORT == "day17-day25") geom_jitter(aes(shape = Day), width = 0.12, height = 0, fill = "white", size = 2.2, show.legend = FALSE)
   else geom_jitter(width = 0.12, height = 0, shape = 21, fill = "white", size = 2)) +
  (if (REPORT == "day17-day25") scale_shape_manual(values = DAY_SHAPE) else NULL) +
  scale_fill_manual(values = cols, guide = "none") + ylim(0, max(0.2, max(d$weighted_NSTI) * 1.15)) +
  labs(x = NULL, y = "weighted NSTI", title = "Distance of each sample's ASVs to their nearest reference genome") +
  theme_classic(base_size = 11) + theme(plot.title = element_text(size = 10))
set.seed(42)   # geom_jitter's positions are drawn when the plot is built
save(p, "picrust2_nsti", 7, 4)

pca_two <- function(ab, lv) {
  # Same PCA as ggpicrust2::pathway_pca (each feature centred and scaled): colour gives the arm, shape
  # the day (Day 17 circles, Day 25 triangles); 95% ellipses per arm x day, solid at Day 25.
  pc <- prcomp(t(ab), center = TRUE, scale. = TRUE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  d <- data.frame(md, PC1 = pc$x[, 1], PC2 = pc$x[, 2])
  d$Day <- ifelse(grepl("Day17", d$Group), "Day 17", "Day 25"); d$Arm <- sub(" Day.*", "", d$Group)
  ggplot(d, aes(PC1, PC2)) +
    stat_ellipse(aes(colour = Arm, linetype = Day, group = Group), level = 0.95, show.legend = FALSE) +
    geom_point(aes(fill = Arm, shape = Day), colour = "grey20", size = 3.6, stroke = 0.6) +
    scale_fill_manual(values = arm) + scale_colour_manual(values = arm) +
    scale_shape_manual(values = DAY_SHAPE) +
    scale_linetype_manual(values = c("Day 17" = "dashed", "Day 25" = "solid")) +
    guides(fill = guide_legend(override.aes = list(shape = 21)), shape = guide_legend(override.aes = list(fill = "white"))) +
    labs(x = sprintf("PC1 (%.1f%%)", ve[1]), y = sprintf("PC2 (%.1f%%)", ve[2]), fill = "Arm", shape = "Day",
         title = sprintf("Predicted %s: PCA of %s features", LEVEL_NAME[[lv]], format(nrow(ab), big.mark = ","))) +
    theme_classic(base_size = 13)
}
for (lv in names(LEVEL_NAME)) {
  ab <- read_ab(lv); ab <- ab[apply(ab, 1, sd) > 0, , drop = FALSE]
  if (REPORT == "day17-day25") { save(pca_two(ab, lv), paste0("picrust2_", lv, "_pca"), 7, 6); next }
  pca <- pathway_pca(abundance = ab, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                     show_marginal = FALSE) +
    ggtitle(sprintf("Predicted %s: PCA of %s features", LEVEL_NAME[[lv]], format(nrow(ab), big.mark = ",")))
  save(pca, paste0("picrust2_", lv, "_pca"), 7, 6)
}

read_rows <- function(lv) {
  d <- read.delim(file.path(OUT, "tables", paste0("PICRUSt2_", lv, "_maaslin3_rows.tsv")),
                  check.names = FALSE, quote = "", stringsAsFactors = FALSE, na.strings = "NA")
  d <- d[d$model == "abundance" & !is.na(d$p), ]
  d$stratum[is.na(d$stratum)] <- ""
  d$description <- decode(d$description)
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
    if (!nrow(d)) return(NULL)
    data.frame(level = lv, row = ROWS$label[k], tag = ROWS$tag[k], id = d$id, description = d$description,
               coef = d$coef, p = d$p, q = d$q, robust = d$robust, stringsAsFactors = FALSE)
  }))
  sel_tab[[lv]] <- pick
  feats <- unique(pick$id)
  ab <- read_ab(lv)
  h <- (sweep(ab, 2, colSums(ab), "/") * 100)[feats, , drop = FALSE]
  tags <- vapply(feats, function(ft) paste(sort(unique(pick$tag[pick$id == ft])), collapse = ","), "")
  rownames(h) <- paste0(short_label(lv, feats, pick$description[match(feats, pick$id)]), "  [", tags, "]")
  hm <- pathway_heatmap(abundance = h, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                        font_size = 9, cluster_rows = TRUE, clustering_method = "average",
                        clustering_distance = "correlation", colorbar_title = "z-score") +
    suppressMessages(ggplot2::scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0, name = "z-score")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(theme = ggplot2::theme(
      legend.key.width = grid::unit(0.5, "cm"), legend.key.height = grid::unit(5, "cm")))) +
    ggplot2::theme(legend.position = "right") +
    ggtitle(sprintf("Predicted %s: the %d strongest per row (%d features)", LEVEL_NAME[[lv]], TOP_N, length(feats)))
  save(hm, paste0("picrust2_", lv, "_heatmap"), 13, 0.24 * length(feats) + 3)
}
write.table(do.call(rbind, sel_tab), file.path(OUT, "tables", "PICRUSt2_heatmap_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

for (lv in names(LEVEL_NAME)) {
  res <- read_rows(lv)
  res$name <- if (lv == "KO") sub(";.*", "", res$description) else
    ifelse(nchar(res$description) > 45, paste0(substr(res$description, 1, 44), "…"), res$description)
  panels <- lapply(seq_len(nrow(ROWS)), function(k) {
    d <- res[res$label == ROWS$label[k], ]
    d$lab <- ""
    r <- which(d$robust %in% TRUE); r <- r[order(d$p[r])][seq_len(min(5, length(r)))]
    d$lab[r] <- d$name[r]
    up <- if (is.na(ROWS$first[k])) "#b2182b" else unname(vcols[ROWS$first[k]])
    dn <- if (is.na(ROWS$second[k])) "#2166ac" else unname(vcols[ROWS$second[k]])
    pathway_volcano(d, fc_col = "coef", p_col = "p", label_col = "lab", fc_threshold = 0,
                    p_threshold = 0.05, label_top_n = 5, point_size = if (lv == "MetaCyc") 1.4 else 0.8,
                    colors = c(Down = dn, "Not Significant" = "grey70", Up = up),
                    title = sprintf("%s [%s]: %d at q < 0.05, %d robust", ROWS$label[k], ROWS$tag[k],
                                    sum(d$call %in% TRUE), sum(d$robust %in% TRUE)),
                    x_lab = if (is.na(ROWS$first[k])) "coefficient (log2)" else
                      sprintf("log2 change in share (right = higher in %s)", ROWS$xname[k]),
                    y_lab = "-log10(raw p)") +
      labs(subtitle = NULL) + theme(plot.title = element_text(size = 9), legend.position = "none")
  })
  vol <- wrap_plots(panels, ncol = NCOL) +
    plot_annotation(title = paste("Predicted", LEVEL_NAME[[lv]], "· MaAsLin3, change in share of the predicted total"),
                    caption = sprintf("Dashed line at raw p = 0.05; colour marks raw p < 0.05; names mark robust calls (q < 0.05 in every leave-one-mouse-out refit)"))
  set.seed(42)   # ggrepel places the labels when the plot is drawn
  save(vol, paste0("picrust2_", lv, "_volcano"), 16, VOL_H)
}
