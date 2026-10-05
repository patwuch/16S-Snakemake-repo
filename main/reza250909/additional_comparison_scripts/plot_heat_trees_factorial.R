#!/usr/bin/env Rscript
# Heat trees (metacoder::heat_tree, drawn with ggplot2) for the four ANCOM-BC2
# contrasts of reza250909's Diet x Radiation factorial, in place of the barplot
# and cladogram pairs, and one for the interaction (see draw_interaction).
# Greengenes2 only.
#
# One tree per contrast and view. The tree is the taxonomy of the taxa the
# panel covers: the 104 taxa ANCOM-BC2 tests, plus that contrast's structural
# zeros (tables/..._structural_zeros.tsv).
#  * Tips of tested taxa are coloured by their ANCOM-BC2 log2FC for the contrast,
#    read from the fit whose reference matches it, on a scale running from the
#    second-named group's colour through white to the first-named group's.
#  * Structural zeros have no log2FC; they are drawn in the colour of the group
#    that has them and labelled with their presence counts and their share of
#    reads in the mice that carry them (median and range, from structural_zeros.py).
#  * Internal nodes are a pale neutral and drawn small: ANCOM-BC2 was fitted to the tips only, and a
#    genus or family is never refitted, so an internal node carries no estimate.
#  * A feature classified no further than genus (or family) gets its own tip
#    under that node (see lineage_tokens), so its colour never reads as an
#    estimate for the whole genus.
# Labels: structural zeros, and tested taxa with q < 0.05 (none of which
# passes ANCOM-BC2's sensitivity check here).
#
# Usage, inside containers/metacoder-env.sif:
#   Rscript plot_heat_trees_factorial.R [--genus-names-only]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")

suppressPackageStartupMessages({
  library(metacoder)
  library(ggplot2)
})
if (packageVersion("metacoder") != "0.3.9")
  stop("metacoder ", packageVersion("metacoder"), " is installed; this script is pinned to 0.3.9. ",
       "Run it in containers/metacoder-env.sif.")

GENUS_ONLY <- "--genus-names-only" %in% commandArgs(trailingOnly = TRUE)
DB    <- "Greengenes2"
STEM  <- paste0("ANCOMBC2_factorial_", DB)
ALPHA <- 0.05

GROUP_COLOUR <- c(Control = "#4D4D4D", IR = "#377EB8", SCFAs = "#E41A1C", "SCFAs+IR" = "#4DAF4A")
INTERNAL_INK <- "#E6E3DC"   # pale warm neutral: the skeleton, kept apart from the grey Control side
MID_INK      <- "#F4F4F2"
RANKS <- c("d", "p", "c", "o", "f", "g", "s")
RANK_WORD <- c(d = "domain", p = "phylum", c = "class", o = "order", f = "family", g = "genus", s = "species")

# key -> fit suffix, coefficient, first-named (positive) group, second-named group, title
PANELS <- list(
  DietSCFAs          = list(fit = "",          coef = "DietSCFAs",   pos = "SCFAs",    neg = "Control",
                            title = "SCFAs vs Control, unirradiated"),
  radrefIR_DietSCFAs = list(fit = "_radrefIR", coef = "DietSCFAs",   pos = "SCFAs+IR", neg = "IR",
                            title = "SCFAs+IR vs IR"),
  RadiationIR        = list(fit = "",          coef = "RadiationIR", pos = "IR",       neg = "Control",
                            title = "IR vs Control, control diet"),
  refSCFAs_RadiationIR = list(fit = "_refSCFAs", coef = "RadiationIR", pos = "SCFAs+IR", neg = "SCFAs",
                            title = "SCFAs+IR vs SCFAs"))

sz_all <- read.delim(file.path(TAB_DIR, paste0(STEM, "_structural_zeros.tsv")), check.names = FALSE, quote = "")

# Lineage tokens for one taxon, ending in a tip of its own. A feature named to
# species ends at its species. Otherwise it gets a tip under its deepest named
# rank, so its colour never reads as an estimate for that whole genus or family.
# Two cases stay apart, because collapsing keeps them as different features:
# an empty token ("g__X;s__", the reference names no species) is labelled
# "g__X (no species call)"; a lineage that stops ("g__X", the classifier's
# confidence ran out) is labelled "g__X (classifier stopped at genus)".
# Returns list(tokens, label).
lineage_tokens <- function(taxon) {
  raw  <- trimws(strsplit(taxon, ";")[[1]])
  raw  <- raw[raw != ""]
  toks <- raw[!grepl("^.__$", raw)]
  last <- toks[length(toks)]
  rank <- substr(last, 1, 1)
  if (rank == "s") return(list(tokens = toks, label = last))
  if (length(raw) > length(toks)) {
    missing <- RANKS[seq(match(rank, RANKS) + 1, length(RANKS))]
    words <- RANK_WORD[missing]
    phrase <- if (length(words) == 1) words else paste(paste(words[-length(words)], collapse = ", "), "or", words[length(words)])
    label <- sprintf("%s (no %s call)", last, phrase)
  } else {
    label <- sprintf("%s (classifier stopped at %s)", last, RANK_WORD[[rank]])
  }
  next_rank <- RANKS[match(rank, RANKS) + 1]
  list(tokens = c(toks, paste0(next_rank, "__", sub("^.__", "", label))), label = label)
}

genus_or_above <- function(taxon) substr(lineage_tokens(taxon)$label, 1, 3) != "s__"

# ---------- shared drawing ----------

# Taxmap over the rows of d (taxon, ...), each ending in its own tip.
# Returns list(obj, ids, row): row[k] is the row of d drawn at node ids[k], NA
# for an internal node.
build_tree <- function(d, key) {
  lt <- lapply(d$taxon, lineage_tokens)
  d$lineage <- vapply(lt, function(x) paste(x$tokens, collapse = ";"), character(1))
  d$label   <- vapply(lt, `[[`, character(1), "label")
  if (anyDuplicated(d$lineage)) stop("duplicate lineages in ", key)
  obj <- parse_tax_data(d, class_cols = "lineage", class_sep = ";",
                        class_regex = "^(.)__(.+)$", class_key = c(rank = "taxon_rank", name = "taxon_name"))
  ids <- taxon_ids(obj)
  tip_of <- obj$data$tax_data$taxon_id                  # one row per taxon, its deepest class
  if (any(tip_of %in% unlist(supertaxa(obj)))) stop("a taxon is an internal node in ", key)
  list(obj = obj, ids = ids, row = match(ids, tip_of), label = d$label)
}

# heat_tree() evaluates its arguments inside the taxmap, where tax_data
# columns (label, lineage, ...) would shadow local vectors of the same name,
# hence the ht_ prefix. ht_colour mixes numbers (scaled on ramp over +-lim)
# and hex codes (drawn as given).
draw_tree <- function(tr, ht_label, ht_colour, ramp, lim, tip_size = NULL, internal = 0.45,
                      node_range = c(0.006, 0.014)) {
  # Tips carry the data, so they are drawn larger than the internal skeleton;
  # sizing internal nodes by how many taxa they hold made the hubs dominate.
  # tip_size (per node, tips only, <= 1) scales tips, e.g. by mice carrying the taxon.
  # internal is the value given to internal nodes; node_range maps the smallest and
  # largest value to diameters (a range proportional to the values keeps tip diameter
  # proportional to tip_size).
  ht_size <- ifelse(is.na(tr$row), internal, if (is.null(tip_size)) 1 else tip_size)
  ht_label_size <- rep(1, length(tr$ids))
  set.seed(1)
  suppressWarnings(heat_tree(tr$obj,
                 node_label = ht_label, node_label_size = ht_label_size,
                 node_size = ht_size, node_color = ht_colour,
                 node_color_trans = "linear", node_color_interval = c(-lim, lim), node_color_range = ramp,
                 edge_color = ht_colour, edge_color_trans = "linear", edge_color_interval = c(-lim, lim),
                 edge_color_range = ramp,
                 node_size_range = node_range, edge_size_range = c(0.0015, 0.004),
                 node_label_size_range = c(0.0115, 0.0115), node_label_max = 200,
                 make_node_legend = FALSE, make_edge_legend = FALSE,
                 layout = "davidson-harel", initial_layout = "reingold-tilford",
                 overlap_avoidance = 2, margin_size = c(0.01, 0.01, 0.01, 0.01)))
}

# Rings round chosen tips. heat_tree keeps each node's circle in its first
# layer, grouped "<taxon id>_node". rings: named list of taxon-id vectors;
# styles: linetype, radius factor and legend text per name.
add_rings <- function(p, rings, styles, fixed_r = NULL) {
  shapes <- p$layers[[1]]$data
  one <- function(tids, name) {
    if (!length(tids)) return(NULL)
    do.call(rbind, lapply(tids, function(tid) {
      pts <- shapes[shapes$group == paste0(tid, "_node"), ]
      cx <- mean(range(pts$x)); cy <- mean(range(pts$y))
      # fixed_r: every ring sized for one reference tip, so a smaller tip sits loose in it.
      r <- (if (is.null(fixed_r)) diff(range(pts$x)) / 2 else fixed_r) * styles[[name]]$radius
      a <- seq(0, 2 * pi, length.out = 61)
      data.frame(x = cx + r * cos(a), y = cy + r * sin(a), group = paste(name, tid), lty = name)
    }))
  }
  df <- do.call(rbind, Map(one, rings, names(rings)))
  if (is.null(df)) return(p)
  p + geom_path(data = df, aes(x = x, y = y, group = group, linetype = lty),
                colour = "#2b2b28", linewidth = 0.45, inherit.aes = FALSE) +
    scale_linetype_manual(values = vapply(styles, `[[`, character(1), "lty"),
                          labels = vapply(styles, `[[`, character(1), "text"),
                          breaks = names(styles), drop = FALSE, name = NULL,
                          guide = guide_legend(ncol = 1, theme = theme(legend.key.width = unit(1, "cm"))))
}

# One colour bar for the tip scale, from an invisible layer; heat_tree draws
# its own fills directly, so this scale touches nothing else.
add_colourbar <- function(p, ramp, lim, name) {
  bb <- ggplot_build(p)$layout$panel_params[[1]]
  legend_df <- data.frame(x = rep(mean(bb$x.range), 2), y = rep(mean(bb$y.range), 2), v = c(-lim, lim))
  p + geom_point(data = legend_df, aes(x = x, y = y, colour = v), alpha = 0, inherit.aes = FALSE) +
    scale_colour_gradientn(colours = ramp, limits = c(-lim, lim), name = name,
                           guide = guide_colourbar(theme = theme(legend.title.position = "top",
                                                                 legend.title = element_text(hjust = 0.5),
                                                                 legend.key.width = unit(9, "cm"),
                                                                 legend.key.height = unit(0.35, "cm"))))
}

# Centre and radius of each node's circle (heat_tree's first layer, grouped
# "<taxon id>_node").
node_geometry <- function(p, tids) {
  shapes <- p$layers[[1]]$data
  do.call(rbind, lapply(tids, function(tid) {
    pts <- shapes[shapes$group == paste0(tid, "_node"), ]
    data.frame(tid = tid, x = mean(range(pts$x)), y = mean(range(pts$y)), r = diff(range(pts$x)) / 2)
  }))
}

# Reference dots for a size-scaled tree, drawn in the tree's own coordinates from the
# radii heat_tree produced, so a key dot is exactly the size of a tip with that count.
add_size_key <- function(p, geom, size_of, counts, title) {
  fit <- lm(r ~ v, data = geom)
  bb <- ggplot_build(p)$layout$panel_params[[1]]
  w <- diff(bb$x.range); h <- diff(bb$y.range)
  r <- predict(fit, data.frame(v = size_of(counts)))
  x <- bb$x.range[2] - w * 0.22 + cumsum(c(0, head(r, -1) + r[-1])) + seq_along(r) * w * 0.02
  y <- bb$y.range[1] + h * 0.06
  a <- seq(0, 2 * pi, length.out = 61)
  circ <- do.call(rbind, lapply(seq_along(r), function(i)
    data.frame(x = x[i] + r[i] * cos(a), y = y + r[i] * sin(a), g = i)))
  p + geom_polygon(data = circ, aes(x = x, y = y, group = g), fill = "#9a9a97", colour = NA, inherit.aes = FALSE) +
    annotate("text", x = x, y = y - max(r) - h * 0.012, label = counts, size = 2.8, colour = "#6b6b66", vjust = 1) +
    annotate("text", x = mean(range(x)), y = y + max(r) + h * 0.012, label = title, size = 2.8, colour = "#6b6b66", vjust = 0)
}

letters_key <- function(n) {
  one <- LETTERS
  c(one, as.vector(outer(one, one, paste0)))[seq_len(n)]
}

# Keys for the named taxa, so full names go beside the tree instead of on it.
# sections: list of list(tids, name (ring style), header, text (one per tid),
# keys = "number" | "letter"). Within a section, keys run clockwise from
# twelve o'clock round the tree's centre, so a reader can sweep for them.
assign_keys <- function(p, sections) {
  all_geo <- node_geometry(p, unique(unlist(lapply(sections, `[[`, "tids"))))
  shapes <- p$layers[[1]]$data
  cx0 <- mean(range(shapes$x)); cy0 <- mean(range(shapes$y))
  out <- lapply(sections, function(sec) {
    if (!length(sec$tids)) return(NULL)
    g <- all_geo[match(sec$tids, all_geo$tid), ]
    ang <- (atan2(g$x - cx0, g$y - cy0) + 2 * pi) %% (2 * pi)
    o <- order(ang)
    k <- if (sec$keys == "number") as.character(seq_along(o)) else letters_key(length(o))
    data.frame(key = k, section = sec$name, header = sec$header, text = sec$text[o],
               tid = g$tid[o], x = g$x[o], y = g$y[o], r = g$r[o], stringsAsFactors = FALSE)
  })
  do.call(rbind, out)
}

# Each key in a small white box just outside its node's ring, up and to the right.
add_key_labels <- function(p, keys) {
  if (is.null(keys) || !nrow(keys)) return(p)
  lab <- data.frame(x = keys$x + keys$r * 1.55, y = keys$y + keys$r * 1.35, key = keys$key)
  # Room past the rim, so a key beside an outermost node is not clipped.
  # heat_tree fixes both the coordinate and the scale limits, so widen both.
  xl <- p$coordinates$limits$x; yl <- p$coordinates$limits$y
  pad <- 0.05 * diff(xl)
  xl <- xl + c(-pad, pad); yl <- yl + c(0, pad)
  p <- suppressMessages(p + scale_x_continuous(expand = c(0, 0), limits = xl) +
                          scale_y_continuous(expand = c(0, 0), limits = yl) +
                          coord_fixed(xlim = xl, ylim = yl))
  p + geom_label(data = lab, aes(x = x, y = y, label = key), inherit.aes = FALSE,
                 size = 2.5, fontface = "bold", colour = "#2b2b28", fill = "white",
                 label.padding = unit(0.06, "lines"), label.r = unit(0.05, "lines"),
                 border.colour = "#9a9a97", linewidth = 0.15, hjust = 0, vjust = 0)
}

# The key list, beside the tree: a header per section, then "key  full name  value".
# KEY_ROWS rows fit one column at a fixed pitch; KEY_TEXT is the text size (mm).
KEY_ROWS <- 48
KEY_TEXT <- 3.1
key_panel <- function(keys) {
  if (is.null(keys) || !nrow(keys))
    return(ggplot() + theme_void())
  rows <- list()
  for (sec in unique(keys$section)) {
    k <- keys[keys$section == sec, ]
    rows[[length(rows) + 1]] <- data.frame(key = "", text = k$header[1], head = TRUE)
    rows[[length(rows) + 1]] <- data.frame(key = k$key, text = k$text, head = FALSE)
    rows[[length(rows) + 1]] <- data.frame(key = "", text = "", head = FALSE)
  }
  df <- do.call(rbind, rows)
  df <- df[-nrow(df), ]
  per_col <- if (nrow(df) > KEY_ROWS) ceiling(nrow(df) / 2) else nrow(df)
  df$col <- (seq_len(nrow(df)) - 1) %/% per_col
  df$y <- -((seq_len(nrow(df)) - 1) %% per_col)
  ggplot(df) +
    geom_text(data = df[!df$head, ], aes(x = col * 1.0, y = y, label = key),
              hjust = 1, size = KEY_TEXT, fontface = "bold", colour = "#2b2b28") +
    geom_text(data = df[!df$head, ], aes(x = col * 1.0 + 0.035, y = y, label = text),
              hjust = 0, size = KEY_TEXT, colour = "#2b2b28") +
    geom_text(data = df[df$head, ], aes(x = col * 1.0 - 0.06, y = y, label = text),
              hjust = 0, size = KEY_TEXT, fontface = "bold", colour = "#6b6b66") +
    scale_x_continuous(limits = c(-0.08, max(df$col) + 1.0), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-max(per_col, KEY_ROWS) - 0.5, 0.8), expand = c(0, 0)) +   # fixed row height
    theme_void()
}

save_tree <- function(p, keys, key, title, subtitle, sidecar, key_extra = 0) {
  p <- add_key_labels(p, keys) +
    theme(legend.position = "bottom", legend.box = "vertical", legend.title = element_text(size = 9, colour = "#2b2b28"),
          legend.text = element_text(size = 8, colour = "#6b6b66"),
          plot.margin = margin(4, 4, 4, 4))
  panel <- key_panel(keys)
  two_cols <- !is.null(keys) && nrow(keys) + 2 * length(unique(keys$section)) - 1 > KEY_ROWS
  # key_extra (inches) widens a one-column key: the interaction tree's carrier
  # counts are longer than any other key line.
  extra <- if (two_cols) 0 else key_extra
  body <- cowplot::plot_grid(p, panel, nrow = 1,
                             rel_widths = if (two_cols) c(1, 1.3) else c(1, 0.75 + extra / 9))
  head <- cowplot::ggdraw() +
    cowplot::draw_label(title, x = 0.5, y = 0.72, size = 13, colour = "#2b2b28") +
    cowplot::draw_label(subtitle, x = 0.5, y = 0.25, size = 7.5, colour = "#6b6b66")
  fig <- cowplot::plot_grid(head, body, ncol = 1, rel_heights = c(0.07, 1))
  width <- if (two_cols) 20 else 15.7 + extra   # the tree keeps its 9 in either way
  out <- file.path(FIG_DIR, paste0(STEM, if (GENUS_ONLY) "_gplus" else "", "_heattree_", key))
  ggsave(paste0(out, ".png"), fig, width = width, height = 9.6, dpi = 200, bg = "white")
  ggsave(paste0(out, ".svg"), fig, width = width, height = 9.6, bg = "white")
  # Sidecars: the counts the caption quotes, and the key list, so the report
  # checker can compare them.
  write.table(sidecar, paste0(out, ".tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  if (!is.null(keys))
    write.table(keys[, c("key", "section", "text")], paste0(out, ".keys.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  basename(out)
}

# ---------- one tree per contrast ----------

draw_panel <- function(key, sp) {
  fit <- read.delim(file.path(TAB_DIR, paste0(STEM, sp$fit, ".tsv")), check.names = FALSE, quote = "")
  tested <- data.frame(taxon = fit$taxon, lfc = fit[[paste0("lfc_", sp$coef)]] / log(2),
                       q = fit[[paste0("q_", sp$coef)]], sz = FALSE, present_in = NA_character_,
                       counts = NA_character_, share = NA_character_, stringsAsFactors = FALSE)
  contrast <- paste(sp$pos, "vs", sp$neg)
  z <- sz_all[sz_all$contrast == contrast, ]
  if (nrow(z) == 0) stop("no structural-zero rows for ", contrast)
  zdf <- data.frame(taxon = z$taxon, lfc = NA_real_, q = NA_real_, sz = TRUE, present_in = z$present_in,
                    counts = ifelse(z$present_in == sp$pos,
                                    sprintf("%d/%d vs 0/%d", z$n_present, z$n_group, z$n_group),
                                    sprintf("0/%d vs %d/%d", z$n_group, z$n_present, z$n_group)),
                    share = z$carrier_share,
                    stringsAsFactors = FALSE)
  if (any(zdf$taxon %in% tested$taxon)) stop("a structural zero is also a tested taxon: ", contrast)
  d <- rbind(tested, zdf)
  # Colour limit from every tested taxon, so a colour means the same log2FC in
  # the species view and the genus view of one contrast.
  lim <- max(abs(d$lfc), na.rm = TRUE)
  if (GENUS_ONLY) d <- d[vapply(d$taxon, genus_or_above, logical(1)), ]
  tr <- build_tree(d, key); row <- tr$row

  ht_colour <- rep(INTERNAL_INK, length(tr$ids))
  is_tested_tip <- !is.na(row) & !d$sz[row]
  is_sz_tip <- !is.na(row) & d$sz[row]
  ht_colour[is_tested_tip] <- format(round(d$lfc[row[is_tested_tip]], 4), scientific = FALSE, trim = TRUE)
  ht_colour[is_sz_tip] <- GROUP_COLOUR[d$present_in[row[is_sz_tip]]]
  ht_colour[is_tested_tip & is.na(d$lfc[row])] <- INTERNAL_INK
  hit <- is_tested_tip & !is.na(d$q[row]) & d$q[row] < ALPHA
  ht_label <- rep(NA_character_, length(tr$ids))       # names go in the key panel

  ramp <- c(GROUP_COLOUR[[sp$neg]], MID_INK, GROUP_COLOUR[[sp$pos]])
  p <- draw_tree(tr, ht_label, ht_colour, ramp, lim)
  # Dashed round a structural zero (no log2FC, so its fill is the group colour
  # at full strength, not an effect size), solid round a tested taxon at q < 0.05.
  p <- add_rings(p, list(hit = tr$ids[hit], sz = tr$ids[is_sz_tip]),
                 list(hit = list(lty = "solid", radius = 1.45, text = sprintf("q < %.2f, fails the sensitivity check", ALPHA)),
                      sz = list(lty = "22", radius = 1.45, text = "structural zero: in one group only, no log2FC")))
  p <- add_colourbar(p, ramp, lim, sprintf("log2FC (ANCOM-BC2)   %s higher ←  → %s higher", sp$neg, sp$pos))
  keys <- assign_keys(p, list(
    list(tids = tr$ids[hit], name = "hit", keys = "number",
         header = sprintf("q < %.2f, fails the sensitivity check (solid ring) · log2FC", ALPHA),
         text = sprintf("%s   %s", tr$label[row[hit]], sub("-", "\u2212", sprintf("%+.2f", d$lfc[row[hit]])))),
    list(tids = tr$ids[is_sz_tip], name = "sz", keys = "letter",
         header = sprintf("Structural zeros (dashed ring) · present %s vs %s · share of reads where present", sp$pos, sp$neg),
         text = sprintf("%s   %s · %s", tr$label[row[is_sz_tip]], d$counts[row[is_sz_tip]], d$share[row[is_sz_tip]]))))
  n_hit <- sum(hit); n_sz <- sum(is_sz_tip); n_tips <- sum(!is.na(row))
  view <- if (GENUS_ONLY) "genus-level names only · " else ""
  name <- save_tree(p, keys, key, sp$title,
    sprintf("%s%s · %d tips · pale nodes: ranks above the tested taxa, no estimate · keyed: %d taxa at q < %.2f (none passes the sensitivity check), %d structural zeros (%s vs %s)",
            view, DB, n_tips, n_hit, ALPHA, n_sz, sp$pos, sp$neg),
    data.frame(pos = sp$pos, neg = sp$neg, tips = n_tips, labelled_hits = n_hit, structural_zeros = n_sz,
               sz_in_pos = sum(d$present_in[row[is_sz_tip]] == sp$pos),
               sz_in_neg = sum(d$present_in[row[is_sz_tip]] == sp$neg), lfc_limit = lim))
  message(sprintf("[heattree] %s: %d tips, %d q < %.2f labelled, %d structural zeros", name, n_tips, n_hit, ALPHA, n_sz))
}

# ---------- the interaction ----------
# Tips: all 135 prevalence-filtered taxa. The 104 in every group are coloured
# by the ANCOM-BC2 interaction coefficient (default fit, log2): right of zero,
# irradiation raises the taxon more, or lowers it less, on SCFAs than on the
# control diet. No group owns either side, so the scale is a neutral
# purple-orange diverging one, not two group colours. The 31 absent from a
# whole group have no abundance interaction and are a flat warm grey. Rings:
# solid for ANCOM-BC2 q < 0.05; dotted on every structural zero, keyed with the
# mice carrying it in each group. Presence is not tested (removed 2026-10-01):
# the dotted rings mark taxa to inspect, not calls.
INTERACTION_RAMP <- c("#5E3C99", "#F7F7F7", "#E66101")
NO_ESTIMATE_INK  <- "#A9A395"

# sized: tip size = mice carrying the taxon, rings fixed at a 32-mouse tip (user, 2026-10-01).
draw_interaction <- function(sized = TRUE) {
  key <- "interaction"; coef <- "DietSCFAs:RadiationIR"
  fit <- read.delim(file.path(TAB_DIR, paste0(STEM, ".tsv")), check.names = FALSE, quote = "")
  un <- read.delim(file.path(TAB_DIR, paste0(STEM, "_untestable_absent_from_a_cell.tsv")),
                   check.names = FALSE, quote = "")
  cells <- grep(" \\(n=", names(un), value = TRUE)
  un$carriers <- apply(un[cells], 1, function(v)
    paste(sprintf("%s %s/%s", sub(" \\(n=.*", "", cells), v, sub(".*n=(\\d+)\\)", "\\1", cells)), collapse = " · "))
  d <- data.frame(taxon = c(fit$taxon, un$taxon), stringsAsFactors = FALSE)
  d$carriers <- un$carriers[match(d$taxon, un$taxon)]
  # Mice carrying each taxon, of 32: tested taxa from their zero counts, structural zeros from
  # their sample counts.
  ix <- read.delim(file.path(TAB_DIR, "DA_interaction_intersection_Greengenes2.tsv"), check.names = FALSE, quote = "")
  n_all <- 32
  d$mice <- ifelse(d$taxon %in% un$taxon, un$n_samples[match(d$taxon, un$taxon)],
                   n_all - ix$n_zero_mice[match(d$taxon, ix$taxon)])
  if (anyNA(d$mice)) stop("mice carrying a taxon missing for ", sum(is.na(d$mice)), " taxa")
  d$lfc <- fit[[paste0("lfc_", coef)]][match(d$taxon, fit$taxon)] / log(2)
  d$q   <- fit[[paste0("q_", coef)]][match(d$taxon, fit$taxon)]
  d$tested <- d$taxon %in% fit$taxon
  lim <- max(abs(d$lfc), na.rm = TRUE)
  if (GENUS_ONLY) d <- d[vapply(d$taxon, genus_or_above, logical(1)), ]
  tr <- build_tree(d, key); row <- tr$row

  tip <- !is.na(row)
  est <- tip & d$tested[row] & !is.na(d$lfc[row])
  ht_colour <- rep(INTERNAL_INK, length(tr$ids))
  ht_colour[est] <- format(round(d$lfc[row[est]], 4), scientific = FALSE, trim = TRUE)
  ht_colour[tip & !d$tested[row]] <- NO_ESTIMATE_INK
  hit <- est & !is.na(d$q[row]) & d$q[row] < ALPHA
  sz <- tip & !d$tested[row]
  ht_label <- rep(NA_character_, length(tr$ids))       # names go in the key panel

  tip_size <- NULL; fixed_r <- NULL
  # Tip diameter proportional to mice carrying the taxon, 45% (fewest) to 100% (all 32);
  # internal nodes drop to 30% so the smallest tips stay distinct from them.
  size_of <- function(n) 0.45 + 0.55 * (n - min(d$mice)) / (n_all - min(d$mice))
  if (sized) {
    tip_size <- rep(NA_real_, length(tr$ids))
    tip_size[tip] <- size_of(d$mice[row[tip]])
    p <- draw_tree(tr, ht_label, ht_colour, INTERACTION_RAMP, lim, tip_size, internal = 0.3,
                   node_range = c(0.014 * 0.3, 0.014))
    geom <- node_geometry(p, tr$ids[tip]); geom$v <- tip_size[tip]
    fixed_r <- max(geom$r)   # the ring around a tip carried by all 32 mice
  } else {
    p <- draw_tree(tr, ht_label, ht_colour, INTERACTION_RAMP, lim)
  }
  p <- add_rings(p, list(hit = tr$ids[hit], sz = tr$ids[sz]), fixed_r = fixed_r, styles =
                 list(hit = list(lty = "solid", radius = 1.45, text = sprintf("abundance interaction q < %.2f (ANCOM-BC2), fails the sensitivity check", ALPHA)),
                      sz = list(lty = "11", radius = 1.75, text = "structural zero: absent from a whole group, not tested")))
  if (sized)
    p <- add_size_key(p, geom, size_of, c(min(d$mice), 16, n_all),
                      sprintf("tip size: mice carrying the taxon, of %d", n_all))
  p <- add_colourbar(p, INTERACTION_RAMP, lim,
                     "Interaction (ANCOM-BC2, log2)   irradiation raises it more on Control ←  → more on SCFAs")
  keys <- assign_keys(p, list(
    list(tids = tr$ids[hit], name = "hit", keys = "number",
         header = sprintf("Abundance interaction q < %.2f, fails the sensitivity check (solid ring) · log2", ALPHA),
         text = sprintf("%s   %s", tr$label[row[hit]], sub("-", "\u2212", sprintf("%+.2f", d$lfc[row[hit]])))),
    list(tids = tr$ids[sz], name = "sz", keys = "letter",
         header = "Structural zero, not tested (dotted ring) · mice carrying it per group",
         text = sprintf("%s   %s", tr$label[row[sz]], d$carriers[row[sz]]))))
  n_tips <- sum(tip); n_est <- sum(est); n_none <- sum(tip & !d$tested[row])
  view <- if (GENUS_ONLY) "genus-level names only · " else ""
  name <- save_tree(p, keys, key, "Does radiation's effect differ between diets?",
    sprintf("%s%s · %d tips: %d coloured by the ANCOM-BC2 interaction, %d grey-brown (absent from a whole group, no abundance estimate) · pale nodes: ranks above the tips\nkeyed: %d at q < %.2f (none passes the sensitivity check), %d structural zeros (dotted, untested)",
            view, DB, n_tips, n_est, n_none, sum(hit), ALPHA, sum(sz)),
    data.frame(tips = n_tips, estimated = n_est, no_estimate = n_none, q_hits = sum(hit), lfc_limit = lim,
               min_mice = min(d$mice[row[tip]]), max_mice = max(d$mice[row[tip]]),
               hits_in_all_mice = sum(d$mice[row[hit]] == n_all)),
    key_extra = if (sum(sz)) max(0, (max(nchar(d$carriers[row[sz]]) + nchar(tr$label[row[sz]])) - 80) * 0.075) else 0)
  message(sprintf("[heattree] %s: %d tips (%d estimated), %d q < %.2f", name, n_tips, n_est, sum(hit), ALPHA))
}

for (key in names(PANELS)) draw_panel(key, PANELS[[key]])
draw_interaction()
