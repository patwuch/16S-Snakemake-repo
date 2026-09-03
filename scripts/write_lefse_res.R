# Writer for LEfSe `.res` files, sourced by run_lefse_MicrobiomeMarker.R.
#
# A `.res` file is a prefix-closed table over the taxonomy tree, one TAB-separated
# row per node:
#
#   lineage.dot.separated \t log10(max(max(class_means),1)) \t class \t |LDA| \t pvalue
#
# Non-significant nodes carry an empty class, an empty LDA and "-" for pvalue.
# Prefix closure (every ancestor of every row is itself a row) is what draws the
# cladogram's background tree; markers only color it. This writer guarantees
# closure structurally, by generating rows from the cumulative prefixes of every
# feature lineage rather than from the marker list.
#
# The reference semantics are save_res() in scripts/lefse/plugin_lefse_run.py.

# microbiomeMarker names a deeper rank's node by prepending its parent's name
# ("g__NK4A214_group_s__muris", or "g__NK4A214_group_s__" when the species is
# unnamed). Reduce it to the plain rank token so all three methods key tree
# nodes by the same vocabulary.
unwrap_composed_rank <- function(tok) sub("^[a-z]__.*_([a-z]__)", "\\1", tok)

# The two renderers parse rows with a bare Python .split(), so a space anywhere
# in a field silently shifts every column after it, and "." is the lineage
# separator. Both are stripped here rather than guarded against downstream.
sanitize_rank_token <- function(tok) {
  prefix <- ifelse(grepl("^[a-z]__", tok), substr(tok, 1, 3), "")
  name   <- sub("^[a-z]__", "", tok)
  # make.unique() suffixes ("...|s__.177") are an artifact of de-duplicating
  # identical lineages, not taxonomy; dropping them merges those rows back onto
  # the one tree node they describe.
  name   <- sub("\\.[0-9]+$", "", name)
  name   <- gsub("[^A-Za-z0-9_]+", "_", name)
  name   <- gsub("_+", "_", name)
  paste0(prefix, gsub("^_|_$", "", name))
}

# A bare rank token like "s__" (species not named beyond genus) has an empty
# name. Since tree nodes are keyed by their exact lineage string, every
# unresolved species under every genus would otherwise collapse onto the *same*
# node. Backfill the blank from the nearest named ancestor so each lineage keeps
# a distinct, readable tip.
fill_unresolved_ranks <- function(parts) {
  last_named <- NA_character_
  for (i in seq_along(parts)) {
    tok    <- parts[i]
    prefix <- sub("^([a-z]__).*", "\\1", tok)
    name   <- sub("^[a-z]__", "", tok)
    if (nzchar(name)) {
      last_named <- name
    } else if (!is.na(last_named)) {
      tag <- if (identical(prefix, "s__")) "sp" else "unresolved"
      tok <- paste0(prefix, last_named, "_", tag)
    }
    parts[i] <- tok
  }
  parts
}

# One delimited lineage string -> the canonical rank-token vector this pipeline
# keys tree nodes by. Run every lineage through here, markers included, so the
# two sides match by construction instead of by coincidence.
canonical_lineage <- function(lineage, delim = "|") {
  parts <- strsplit(lineage, delim, fixed = TRUE)[[1]]
  parts <- parts[nzchar(parts)]
  # Classifiers disagree on the top rank's prefix: GTDB-derived ones emit "d__"
  # (domain) while classic LEfSe/MetaPhlAn lineages use "k__" (kingdom). Settle
  # on one so a lineage keys the same tree node whichever produced it.
  parts <- sub("^k__", "d__", parts)
  fill_unresolved_ranks(sanitize_rank_token(unwrap_composed_rank(parts)))
}

canonical_key <- function(lineage, delim = "|") {
  paste(canonical_lineage(lineage, delim), collapse = ".")
}

# formatC() pads a numeric vector to a common width, and one space anywhere in
# a row shifts every column after it once a renderer splits on whitespace.
num_field <- function(x) trimws(formatC(as.numeric(x), format = "g", digits = 15))

# counts       numeric matrix, features (tree leaves) x samples, raw counts
# group_vec    class of each column of `counts`, in column order
# lineage_map  list of rank-token vectors, one per row of `counts` (row order)
# markers      data.frame(lineage, class, lda, pvalue) keyed by dot-joined
#              lineage; NULL or 0 rows when nothing was significant
# path         output file
write_lefse_res <- function(counts, group_vec, lineage_map, markers, path) {
  stopifnot(ncol(counts) == length(group_vec),
            nrow(counts) == length(lineage_map))

  depths <- colSums(counts, na.rm = TRUE)
  if (any(depths <= 0))
    stop("Sample(s) with zero total counts: ",
         paste(colnames(counts)[depths <= 0], collapse = ", "))
  # LEfSe's own input formatter normalizes to counts-per-million (`-o 1000000`)
  # before it sums the hierarchy, and column 2 is thresholded at log10(1), so the
  # scale is not cosmetic.
  cpm <- sweep(counts, 2, depths, "/") * 1e6

  nodes <- new.env(hash = TRUE, parent = emptyenv())
  for (i in seq_len(nrow(cpm))) {
    parts <- lineage_map[[i]]
    row   <- cpm[i, ]
    for (d in seq_along(parts)) {
      key  <- paste(parts[seq_len(d)], collapse = ".")
      prev <- nodes[[key]]
      nodes[[key]] <- if (is.null(prev)) row else prev + row
    }
  }
  if (length(ls(nodes)) == 0) stop("No tree nodes to write: every lineage was empty.")

  keys <- sort(ls(nodes), method = "radix")
  abund <- do.call(rbind, mget(keys, envir = nodes))

  classes <- if (is.factor(group_vec)) levels(group_vec) else sort(unique(group_vec))
  classes <- classes[classes %in% as.character(group_vec)]
  means <- vapply(classes, function(cl)
    rowMeans(abund[, as.character(group_vec) == cl, drop = FALSE]),
    numeric(nrow(abund)))
  value <- log10(pmax(apply(means, 1, max), 1))

  cls <- rep("", length(keys))
  lda <- rep("", length(keys))
  pval <- rep("-", length(keys))

  if (!is.null(markers) && nrow(markers) > 0) {
    markers <- markers[order(-abs(as.numeric(markers$lda))), , drop = FALSE]
    # Distinct input features can canonicalize onto one node, since species left
    # unresolved under one genus all backfill to the same name. The strongest
    # effect wins that node.
    markers <- markers[!duplicated(markers$lineage), , drop = FALSE]

    idx <- match(markers$lineage, keys)
    orphans <- markers$lineage[is.na(idx)]
    if (length(orphans) > 0)
      warning("Dropping ", length(orphans), " marker(s) with no matching tree node: ",
              paste(utils::head(orphans, 5), collapse = ", "))

    keep <- !is.na(idx)
    cls[idx[keep]]  <- as.character(markers$class[keep])
    lda[idx[keep]]  <- num_field(abs(as.numeric(markers$lda[keep])))
    pval[idx[keep]] <- ifelse(is.na(markers$pvalue[keep]), "-",
                              num_field(markers$pvalue[keep]))
  }

  fields <- list(keys, num_field(value), cls, lda, pval)
  if (any(vapply(fields, function(f) any(grepl("[[:space:]]", f)), logical(1))))
    stop("Whitespace in a .res field would shift every column after it once the ",
         "renderers split the row.")

  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  writeLines(do.call(paste, c(fields, sep = "\t")), path)
  message("[lefse_res] Wrote ", length(keys), " tree node(s), ",
          sum(nzchar(cls)), " significant, to ", path)
}
