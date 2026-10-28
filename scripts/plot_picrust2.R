if (exists("snakemake")) {
  # Send output and messages (incl. errors) to the Snakemake log.
  log_con <- file(snakemake@log[[1]], open = "wt")
  sink(log_con)
  sink(log_con, type = "message")
}

library(tidyverse)

# ===========================================================================
# PICRUSt2 functional heatmaps (MetaCyc pathway / KO / EC)
#
# Output:
#   Three heatmaps, each written as BOTH a vector .svg and a 300-dpi .png:
#       picrust2_MetaCyc.svg / .png
#       picrust2_KO.svg      / .png
#       picrust2_EC.svg      / .png
#   (six files total; no bundled PDF)
#
# Normalization (prep_scores):
#   transform = "clr"  -> centred log-ratio; compositional-safe (default)
#   transform = "cpm"  -> log10 relative abundance (depth-corrected only)
#
# Axis ordering:
#   x (samples)  -> controlled by `sample_order`:
#        "alpha"            natural alphanumeric sort of sampleid (default;
#                           T2 < T10, not string-sorted)
#        "cluster"          hierarchical clustering by functional profile
#        <metadata column>  order by that metadata field, e.g. "Order"
#                           (numeric columns sort numerically; final tiebreak is
#                           a natural alphanumeric sort of sampleid)
#   y (features) -> always hierarchically clustered per plot, so co-varying
#                   functions sit together (no natural order for descriptions).
# ===========================================================================


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

pivot_picrust_long <- function(data_table) {
  id_col   <- colnames(data_table)[1]
  name_col <- colnames(data_table)[2]
  out <- data_table %>%
    pivot_longer(cols = -c(1:2), names_to = "sampleid", values_to = "Abundance")
  names(out)[match(c(id_col, name_col), names(out))] <- c("feature_id", "feature_name")
  out
}

# Coerce a column to numeric ONLY if every non-NA value parses cleanly
# (so an integer column read as text, e.g. "Order", sorts 2 < 10 not "10" < "2");
# genuine text columns are returned unchanged.
maybe_num <- function(x) {
  if (is.character(x)) {
    xn <- suppressWarnings(as.numeric(x))
    if (!any(is.na(xn) & !is.na(x))) return(xn)
  }
  x
}

# Per-sample-normalized `score`. Computed on the FULL feature set (required for
# CLR's per-sample geometric mean), before any top-N filtering.
prep_scores <- function(data_table, transform = c("clr", "cpm"), pseudocount = 1) {
  transform <- match.arg(transform)
  long <- pivot_picrust_long(data_table)
  if (transform == "cpm") {
    long %>%
      group_by(sampleid) %>%
      mutate(score = if (sum(Abundance) > 0)
                       log10(Abundance / sum(Abundance) * 1e6 + 1) else 0) %>%
      ungroup()
  } else {  # clr
    long %>%
      group_by(sampleid) %>%
      mutate(.l = log(Abundance + pseudocount), score = .l - mean(.l)) %>%
      ungroup() %>% select(-.l)
  }
}

save_plot <- function(plot, svg_path, png_path, width = 12, height = 10, dpi = 300) {
  svg_device <- if (requireNamespace("svglite", quietly = TRUE)) "svg" else grDevices::svg
  ggsave(svg_path, plot, device = svg_device, width = width, height = height,
         bg = "white")
  ggsave(png_path, plot, device = "png",      width = width, height = height, dpi = dpi,
         bg = "white")
  message(sprintf("  wrote %s  +  %s (%d dpi)", svg_path, png_path, dpi))
}

# ---------------------------------------------------------------------------
# Global sample (x-axis) ordering, shared across all three plots
# ---------------------------------------------------------------------------
get_global_sample_order <- function(file_path, metadata_path, sample_order = "alpha",
                                    transform = "clr", pseudocount = 1) {
  data_table <- read_tsv(file_path, col_types = cols())
  samples <- colnames(data_table)[-(1:2)]

  # 1. Cluster by functional profile (original behaviour)
  if (identical(sample_order, "cluster")) {
    mat <- data_table %>%
      prep_scores(transform, pseudocount) %>%
      select(feature_id, sampleid, score) %>%
      pivot_wider(names_from = "sampleid", values_from = "score", values_fill = 0) %>%
      column_to_rownames("feature_id") %>% as.matrix()
    return(colnames(mat)[hclust(dist(t(mat)))$order])
  }

  # 2. Natural alphanumeric sort of the sample IDs
  if (identical(sample_order, "alpha")) {
    return(str_sort(samples, numeric = TRUE))
  }

  # 3. Otherwise treat `sample_order` as one or more metadata column names.
  #    Numeric-looking columns (e.g. "Order") sort numerically; the final
  #    tiebreak is a natural alphanumeric sort of sampleid.
  md <- read_tsv(metadata_path, col_types = cols()) %>% filter(sampleid %in% samples)
  missing <- setdiff(sample_order, colnames(md))
  if (length(missing) > 0)
    stop("sample_order column(s) not found in metadata: ", paste(missing, collapse = ", "))

  md <- md %>%
    mutate(across(all_of(sample_order), maybe_num),
           .natkey = match(sampleid, str_sort(sampleid, numeric = TRUE)))
  ord <- md %>% arrange(across(all_of(sample_order)), .natkey) %>% pull(sampleid)
  c(ord, str_sort(setdiff(samples, ord), numeric = TRUE))  # any un-annotated samples last
}


# ---------------------------------------------------------------------------
# Build one heatmap
# ---------------------------------------------------------------------------
generate_plot <- function(file_path, category_label, metadata_path,
                          global_sample_order, top_n,
                          transform = "clr", pseudocount = 1) {

  message("Reading file: ", file_path)
  data_table <- read_tsv(file_path, col_types = cols())
  metadata   <- read_tsv(metadata_path, col_types = cols())
  message(sprintf("Loaded %s: %d features x %d samples [transform=%s]",
                  category_label, nrow(data_table), ncol(data_table) - 2, transform))

  df <- prep_scores(data_table, transform = transform, pseudocount = pseudocount)

  df_merged <- inner_join(df, metadata, by = "sampleid")
  if (nrow(df_merged) == 0) stop("No matching sampleids between metadata and PICRUSt2 output.")

  feature_labels <- df_merged %>%
    distinct(feature_id, feature_name) %>%
    group_by(feature_id) %>% slice(1) %>% ungroup()
  label_lookup <- setNames(feature_labels$feature_name, feature_labels$feature_id)

  # Top-N most variable features, on the transformed score
  top_features <- df_merged %>%
    group_by(feature_id) %>%
    summarize(v = var(score), .groups = "drop") %>%
    slice_max(v, n = top_n, with_ties = FALSE) %>%
    pull(feature_id)
  plot_data <- df_merged %>% filter(feature_id %in% top_features)

  # Cluster features (y-axis), keyed on the unique id
  wide_mat <- plot_data %>%
    select(feature_id, sampleid, score) %>%
    pivot_wider(names_from = "sampleid", values_from = "score", values_fill = 0) %>%
    column_to_rownames("feature_id") %>% as.matrix()
  feature_order <- rownames(wide_mat)[hclust(dist(wide_mat))$order]

  # Per-feature z-score for display
  plot_data <- plot_data %>%
    group_by(feature_id) %>%
    mutate(z = { s <- sd(score)
                 if (is.na(s) || s == 0) rep(0, n()) else (score - mean(score)) / s }) %>%
    ungroup()

  plot_data$sampleid   <- factor(plot_data$sampleid,   levels = global_sample_order)
  plot_data$feature_id <- factor(plot_data$feature_id, levels = feature_order)

  ggplot(plot_data, aes(x = sampleid, y = feature_id, fill = z)) +
    geom_tile() +
    scale_y_discrete(labels = label_lookup) +
    scale_fill_gradient2(low = "#3B4CC0", mid = "grey95", high = "#B40426",
                         midpoint = 0, name = sprintf("Row z-score\n(%s)", transform)) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8),
      axis.text.y = element_text(size = 7),
      plot.title  = element_text(hjust = 0.5)
    ) +
    labs(title = sprintf("Top %d %s", top_n, category_label),
         x = "sampleid", y = "Description")
}


# ---------------------------------------------------------------------------
# Main execution (Snakemake-driven; guarded so functions are testable alone)
#
# Expected Snakemake wiring (six named outputs, two per category):
#
#   rule picrust2_heatmaps:
#       input:
#           pathway  = "...path_abun_unstrat_descrip.tsv",   # MetaCyc pathways
#           ko       = "...ko_pred_metagenome_unstrat_descrip.tsv",
#           ec       = "...ec_pred_metagenome_unstrat_descrip.tsv",
#           metadata = "metadata.tsv"
#       output:
#           metacyc_svg = "results/picrust2_MetaCyc.svg",
#           metacyc_png = "results/picrust2_MetaCyc.png",
#           ko_svg      = "results/picrust2_KO.svg",
#           ko_png      = "results/picrust2_KO.png",
#           ec_svg      = "results/picrust2_EC.svg",
#           ec_png      = "results/picrust2_EC.png"
#       params:
#           top_n = 30, transform = "clr", sample_order = "Order",
#           width = 12, height = 10, dpi = 300
#       script: "picrust2_heatmaps.R"
# ---------------------------------------------------------------------------
if (exists("snakemake")) {

  top_n        <- snakemake@params[["top_n"]];        if (is.null(top_n))        top_n        <- 30
  transform    <- snakemake@params[["transform"]];    if (is.null(transform))    transform    <- "clr"
  sample_order <- snakemake@params[["sample_order"]]; if (is.null(sample_order)) sample_order <- "Order"
  width        <- snakemake@params[["width"]];        if (is.null(width))        width        <- 12
  height       <- snakemake@params[["height"]];       if (is.null(height))       height       <- 10
  dpi          <- snakemake@params[["dpi"]];          if (is.null(dpi))          dpi          <- 300

  global_sample_order <- get_global_sample_order(
    snakemake@input[["pathway"]], snakemake@input[["metadata"]],
    sample_order = sample_order, transform = transform)

  # category -> (input table, plot title, svg output, png output)
  plots <- list(
    MetaCyc = list(input = snakemake@input[["pathway"]], label = "MetaCyc Pathways",
                   svg = snakemake@output[["metacyc_svg"]], png = snakemake@output[["metacyc_png"]]),
    KO      = list(input = snakemake@input[["ko"]],       label = "KOs",
                   svg = snakemake@output[["ko_svg"]],      png = snakemake@output[["ko_png"]]),
    EC      = list(input = snakemake@input[["ec"]],       label = "EC Numbers",
                   svg = snakemake@output[["ec_svg"]],      png = snakemake@output[["ec_png"]])
  )

  for (nm in names(plots)) {
    p   <- plots[[nm]]
    plt <- generate_plot(p$input, p$label, snakemake@input[["metadata"]],
                         global_sample_order, top_n, transform)
    save_plot(plt, p$svg, p$png, width = width, height = height, dpi = dpi)
  }
  message("Success: 6 files written (svg + png for MetaCyc, KO, EC)")
}