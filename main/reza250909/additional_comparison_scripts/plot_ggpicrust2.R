#!/usr/bin/env Rscript
# ggpicrust2 figures for the PICRUSt2 section of reza250909's factorial report.
#   picrust2_<level>_pca       - ggpicrust2::pathway_pca per level (MetaCyc, KO, EC)
#   picrust2_<level>_heatmap   - ggpicrust2::pathway_heatmap of the 10 strongest
#                                features per contrast (>= 16-mice refit), pooled
#   picrust2_<level>_volcano   - ggpicrust2::pathway_volcano per contrast on the
#                                all-feature LinDA results (run_ggpicrust2_linda.R),
#                                against raw p, for MetaCyc, KO and EC
#   picrust2_linda_filter_sensitivity - q < 0.05 counts in the two SCFAs+IR
#                                contrasts as the prevalence filter changes
# Group colours come from the study metadata, as in every other figure.
#
# Usage, in containers/ggpicrust2-env.sif:
#   Rscript plot_ggpicrust2.R

suppressPackageStartupMessages({ library(ggpicrust2); library(ggplot2); library(patchwork) })
THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
OUT_DIR  <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Radiation-factorial")
FIG_DIR  <- file.path(OUT_DIR, "figures")

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

CONTRASTS <- c("SCFAs vs Control", "SCFAs+IR vs IR", "IR vs Control", "SCFAs+IR vs SCFAs")
LEVEL_NAME <- c(MetaCyc = "MetaCyc pathways", KO = "KEGG orthologs (KO)", EC = "enzymes (EC)")
read_linda <- function(lv, suffix = "")
  read.delim(file.path(OUT_DIR, "tables", paste0("PICRUSt2_", lv, "_linda", suffix, ".tsv")),
             check.names = FALSE, quote = "")

# PCA per level (ggpicrust2::pathway_pca centres and scales each feature, so
# features with no variance across the 32 mice are dropped first).
for (lv in names(LEVEL_NAME)) {
  ab <- read_ab(lv)
  ab <- ab[apply(ab, 1, sd) > 0, , drop = FALSE]
  pca <- pathway_pca(abundance = ab, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                     show_marginal = FALSE) +
    ggtitle(sprintf("Predicted %s: PCA of %s features", LEVEL_NAME[[lv]], format(nrow(ab), big.mark = ",")))
  save(pca, paste0("picrust2_", lv, "_pca"), 7, 6)
}

# Heatmap per level: the 10 features with the smallest raw p in each contrast, from
# the refit among features predicted in >= 16 mice, pooled (duplicates once).
# Rows are each mouse's relative abundance, z-scored per feature by
# ggpicrust2::pathway_heatmap; rows clustered, columns grouped. Rows are labelled
# by description (KO: gene symbol and name), not by ID.
short_label <- function(lv, id, desc) {
  d <- if (lv == "KO") sub("^([^;]*); *(.*?)( \\[EC:.*\\])?$", "\\1: \\2", desc) else desc
  # Neither ggpicrust2's nor PICRUSt2's KO map names some IDs ("not_found").
  if (lv == "KO") d <- ifelse(desc == "not_found", paste0(sub("^ko:", "", id), ": no name in either KO map"), d)
  if (lv == "EC") {
    # Reassigned EC numbers keep only "Transferred entry: <new>" as their name.
    moved <- grepl("^Transferred entry: ", d)
    d <- ifelse(moved, paste0("EC ", id, " (now ", sub("^Transferred entry: *", "", d), ")"), paste0(d, " (", id, ")"))
  }
  d <- ifelse(nchar(d) > 60, paste0(substr(d, 1, 59), "…"), d)
  ifelse(duplicated(d) | duplicated(d, fromLast = TRUE), paste0(d, " (", id, ")"), d)
}
sel_tab <- list()
for (lv in names(LEVEL_NAME)) {
  res <- read_linda(lv, "_minpres16")
  pick <- do.call(rbind, lapply(CONTRASTS, function(k) {
    d <- res[res$contrast == k, ]; d <- d[order(d$p), ][1:10, ]
    data.frame(level = lv, contrast = k, feature = d$feature, description = d$description, p = d$p, q = d$q, log2fc = d$log2fc)
  }))
  sel_tab[[lv]] <- pick
  feats <- unique(pick$feature)
  ab <- read_ab(lv)
  rel <- sweep(ab, 2, colSums(ab), "/") * 100
  h <- rel[feats, , drop = FALSE]
  # Tag each row with the contrast(s) that picked it: A-D, in CONTRASTS order.
  tags <- vapply(feats, function(ft) paste(LETTERS[match(unique(pick$contrast[pick$feature == ft]), CONTRASTS)] |>
                                             sort(), collapse = ","), "")
  rownames(h) <- paste0(short_label(lv, feats, pick$description[match(feats, pick$feature)]), "  [", tags, "]")
  hm <- pathway_heatmap(abundance = h, metadata = md, group = "Group", colors = unname(cols[levels(md$Group)]),
                        font_size = 9, cluster_rows = TRUE, clustering_method = "average",
                        clustering_distance = "correlation", colorbar_title = "z-score") +
    # ggpicrust2 passes breaks = colorbar_breaks (NULL by default) to its fill scale,
    # and breaks = NULL removes the colour bar. Same colours, default breaks.
    suppressMessages(ggplot2::scale_fill_gradient2(low = "#0571b0", mid = "white", high = "#ca0020", midpoint = 0,
                                                   name = "z-score")) +
    ggplot2::guides(fill = ggplot2::guide_colourbar(theme = ggplot2::theme(
      legend.key.width = grid::unit(0.5, "cm"), legend.key.height = grid::unit(5, "cm")))) +
    ggplot2::theme(legend.position = "right") +
    ggtitle(sprintf("Predicted %s: the 10 strongest per contrast (%d features)", LEVEL_NAME[[lv]], length(feats)))
  save(hm, paste0("picrust2_", lv, "_heatmap"), 13, 0.24 * length(feats) + 3)
}
write.table(do.call(rbind, sel_tab), file.path(OUT_DIR, "tables", "PICRUSt2_heatmap_features.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# One volcano figure per level, from the all-feature fits (ggpicrust2's default).
for (lv in names(LEVEL_NAME)) {
  res <- read_linda(lv)
  # KO descriptions are long "SYMBOL; name [EC]" strings: label with the symbol.
  res$label <- if (lv == "KO") sub(";.*", "", res$description) else
    ifelse(nchar(res$description) > 45, paste0(substr(res$description, 1, 44), "…"), res$description)
  panels <- lapply(CONTRASTS, function(k) {
    d <- res[res$contrast == k, ]
    g <- strsplit(k, " vs ")[[1]]
    n <- sum(d$p < 0.05)
    pathway_volcano(d, fc_col = "log2fc", p_col = "p", label_col = "label", fc_threshold = 0,
                    p_threshold = 0.05, label_top_n = 5, point_size = if (lv == "MetaCyc") 1.6 else 0.9,
                    colors = c(Down = unname(cols[g[2]]), "Not Significant" = "grey70", Up = unname(cols[g[1]])),
                    title = sprintf("%s: %d of %d at raw p < 0.05, %d at q < 0.05", k, n, nrow(d), sum(d$q < 0.05)),
                    x_lab = sprintf("log2 fold change (right = higher in %s)", g[1]), y_lab = "-log10(raw p)") +
      theme(plot.title = element_text(size = 10), legend.position = "none")
  })
  vol <- wrap_plots(panels, ncol = 2) +
    plot_annotation(title = paste("Predicted", LEVEL_NAME[[lv]]),
                    caption = "LinDA via ggpicrust2::pathway_daa, all features; dashed line at raw p = 0.05")
  save(vol, paste0("picrust2_", lv, "_volcano"), 13, 10)
}

# How the q < 0.05 counts move with a prevalence filter, per level.
sens <- do.call(rbind, lapply(names(LEVEL_NAME), function(lv)
  do.call(rbind, lapply(c(0L, 8L, 12L, 16L), function(mp) {
    d <- read_linda(lv, if (mp) paste0("_minpres", mp) else "")
    do.call(rbind, lapply(c("SCFAs+IR vs IR", "SCFAs+IR vs SCFAs"), function(k)
      data.frame(level = LEVEL_NAME[[lv]], min_presence = mp, contrast = k,
                 n_q = sum(d$q[d$contrast == k] < 0.05), n_tested = sum(d$contrast == k))))
  }))))
sens$level <- factor(sens$level, levels = LEVEL_NAME)
lab <- sens[sens$min_presence == 16, ]
p <- ggplot(sens, aes(min_presence, n_q, linetype = contrast, colour = contrast)) +
  geom_line(linewidth = 0.8) + geom_point(size = 2.4) +
  geom_text(data = lab, aes(label = contrast), hjust = -0.08, size = 3, show.legend = FALSE) +
  facet_wrap(~ level, scales = "free_y") +
  scale_colour_manual(values = c("SCFAs+IR vs IR" = "#14181d", "SCFAs+IR vs SCFAs" = "#6b7280")) +
  scale_linetype_manual(values = c("SCFAs+IR vs IR" = "solid", "SCFAs+IR vs SCFAs" = "dashed")) +
  scale_x_continuous(breaks = c(0, 8, 12, 16), labels = c("all", "8", "12", "16"), expand = expansion(mult = c(0.05, 0.9))) +
  labs(x = "features kept: predicted in at least this many of 32 mice", y = "features at q < 0.05",
       title = "LinDA q < 0.05 counts depend on which features are fitted",
       subtitle = "SCFAs vs Control and IR vs Control have none at any filter") +
  theme_minimal(base_size = 11) + theme(legend.position = "none", panel.grid.minor = element_blank())
save(p, "picrust2_linda_filter_sensitivity", 12, 4.2)
write.table(sens, file.path(OUT_DIR, "tables", "PICRUSt2_linda_filter_sensitivity.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
