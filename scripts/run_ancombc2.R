suppressPackageStartupMessages({
  library(biomformat)
  library(phyloseq)
  library(ANCOMBC)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}

natural_level_order <- function(x) {
  pad_nums <- function(s) {
    parts <- strsplit(s, "(?<=\\D)(?=\\d)|(?<=\\d)(?=\\D)", perl = TRUE)[[1]]
    paste(ifelse(grepl("^\\d+$", parts),
                 formatC(as.integer(parts), width = 10, flag = "0"),
                 parts), collapse = "")
  }
  x[order(vapply(as.character(x), pad_nums, character(1)))]
}

extract_taxonomy_matrix <- function(taxonomy_df, taxa_ids) {
  tax_col <- if ("Taxon" %in% colnames(taxonomy_df)) {
    "Taxon"
  } else if ("taxonomy" %in% colnames(taxonomy_df)) {
    "taxonomy"
  } else {
    NULL
  }

  if (is.null(tax_col)) {
    mat <- matrix("", nrow = length(taxa_ids), ncol = 1)
    rownames(mat) <- taxa_ids
    colnames(mat) <- "Taxon"
    return(mat)
  }

  tax_vec <- taxonomy_df[[tax_col]]
  names(tax_vec) <- rownames(taxonomy_df)
  tax_vec <- as.character(tax_vec[taxa_ids])
  tax_vec[is.na(tax_vec)] <- ""

  split_tax <- strsplit(tax_vec, ";\\s*")
  max_rank <- max(vapply(split_tax, length, integer(1)), 1L)

  mat <- matrix("", nrow = length(taxa_ids), ncol = max_rank)
  rownames(mat) <- taxa_ids
  colnames(mat) <- paste0("Rank", seq_len(max_rank))

  for (i in seq_along(split_tax)) {
    vec <- split_tax[[i]]
    if (length(vec) > 0) {
      mat[i, seq_len(length(vec))] <- vec
    }
  }

  mat
}

# ---------- LEfSe-style marker selection ----------
# ANCOM-BC2 has no LDA score and no single "enriched group" column; the .res
# format wants both. Derive them from the two-layer FDR control ancombc2()
# already performed: the omnibus test decides WHETHER a taxon differs, the
# mdFDR-corrected pairwise family decides BETWEEN WHICH levels.
#
# Returns data.frame(taxon, class, lda, pvalue); `taxon` is the feature id, which
# the caller maps onto its canonical lineage key.
build_ancombc2_markers <- function(res, group_col, group_levels, ref_level, alpha) {
  empty <- data.frame(taxon = character(0), class = character(0),
                      lda = numeric(0), pvalue = numeric(0),
                      stringsAsFactors = FALSE)

  res_main <- as.data.frame(res$res)
  if (is.null(res_main) || nrow(res_main) == 0) return(empty)
  res_global <- if (is.null(res$res_global)) NULL else as.data.frame(res$res_global)
  res_pair   <- if (is.null(res$res_pair))   NULL else as.data.frame(res$res_pair)

  non_ref <- setdiff(group_levels, ref_level)
  term_of <- setNames(paste0(group_col, group_levels), group_levels)

  # Contrast column suffixes are BUILT from the model's term names, never parsed
  # back out of them: a level name may itself contain "_", so splitting
  # "lfc_GroupC_GroupB" into its two sides would be guesswork.
  # Sign convention (verified empirically against ancombc2): lfc_<A>_<B> > 0
  # means level A is the higher-abundance side; for a vs-reference term
  # lfc_<A> > 0 means A is higher than the reference level.
  contrasts <- list()
  for (lvl in non_ref)
    contrasts[[length(contrasts) + 1L]] <-
      list(suffix = term_of[[lvl]], pos = lvl, neg = ref_level)
  for (a in non_ref) for (b in non_ref) if (a != b)
    contrasts[[length(contrasts) + 1L]] <-
      list(suffix = paste0(term_of[[a]], "_", term_of[[b]]), pos = a, neg = b)

  # Strongest |lfc| among the contrasts that clear alpha in `tbl` for one taxon.
  pick_best <- function(tbl, row_idx, cand) {
    best <- NULL
    for (ct in cand) {
      lfc_col <- paste0("lfc_", ct$suffix)
      q_col   <- paste0("q_",   ct$suffix)
      if (!all(c(lfc_col, q_col) %in% colnames(tbl))) next
      lfc <- tbl[[lfc_col]][row_idx]
      q   <- tbl[[q_col]][row_idx]
      if (is.na(lfc) || is.na(q) || q >= alpha) next
      if (is.null(best) || abs(lfc) > abs(best$lfc))
        best <- list(class = if (lfc > 0) ct$pos else ct$neg, lfc = lfc, q = q)
    }
    best
  }

  # res$res only ever carries vs-reference model terms.
  main_cand <- Filter(function(ct) identical(ct$neg, ref_level), contrasts)

  # The omnibus test is the gate. A two-level design has no omnibus test to run,
  # so ancombc2 returns nothing here and the pairwise q-values are the only FDR
  # control there is — gate on every tested taxon instead of on nothing.
  gate <- if (!is.null(res_global) && nrow(res_global) > 0) {
    res_global$taxon[!is.na(res_global$q_val) & res_global$q_val < alpha]
  } else {
    res_main$taxon
  }

  out <- list()
  for (tx in gate) {
    best <- NULL
    if (!is.null(res_pair)) {
      i <- match(tx, res_pair$taxon)
      if (!is.na(i)) best <- pick_best(res_pair, i, contrasts)
    }
    if (is.null(best)) {
      i <- match(tx, res_main$taxon)
      if (!is.na(i)) best <- pick_best(res_main, i, main_cand)
    }
    # Omnibus-significant but no contrast survives: the mdFDR pairwise family
    # never localized the effect to a specific pair of levels, so there is no
    # single group to color the taxon for. Dropping it is the conservative read.
    if (is.null(best)) next
    # ancombc2 models log(count) with R's log(), so lfc is a NATURAL-log fold
    # change. Divide by log(2) to publish the log2 fold change the effect-size
    # axis names and the field reports. Ranking is unaffected; only the scale is.
    out[[length(out) + 1L]] <- data.frame(
      taxon = tx, class = best$class, lda = abs(best$lfc) / log(2), pvalue = best$q,
      stringsAsFactors = FALSE)
  }
  if (length(out) == 0) return(empty)
  do.call(rbind, out)
}

# ---------- BIOM format helpers ----------
.is_hdf5_biom <- function(path) {
  magic <- tryCatch(readBin(path, raw(), n = 4), error = function(e) raw(0))
  length(magic) >= 4 &&
    identical(magic[1:4], as.raw(c(0x89, 0x48, 0x44, 0x46)))  # \x89HDF
}

.read_hdf5_biom <- function(path) {
  if (!requireNamespace("rhdf5", quietly = TRUE))
    stop("Package 'rhdf5' is required for HDF5 BIOM files.\n",
         "  Install with: BiocManager::install('rhdf5')")
  obs_ids  <- as.character(rhdf5::h5read(path, "observation/ids"))
  samp_ids <- as.character(rhdf5::h5read(path, "sample/ids"))
  data_v   <- rhdf5::h5read(path, "observation/matrix/data")
  indices  <- as.integer(rhdf5::h5read(path, "observation/matrix/indices")) + 1L
  indptr   <- as.integer(rhdf5::h5read(path, "observation/matrix/indptr"))  + 1L
  mat <- matrix(0L, nrow = length(obs_ids), ncol = length(samp_ids),
                dimnames = list(obs_ids, samp_ids))
  for (i in seq_along(obs_ids)) {
    col_idx        <- indices[indptr[i]:(indptr[i + 1L] - 1L)]
    mat[i, col_idx] <- data_v[indptr[i]:(indptr[i + 1L] - 1L)]
  }
  mat
}

.read_biom_auto <- function(path) {
  if (.is_hdf5_biom(path)) {
    message("[ANCOMBC2] Detected HDF5 (BIOM v2) format — using rhdf5 reader")
    .read_hdf5_biom(path)
  } else {
    message("[ANCOMBC2] Detected JSON (BIOM v1) format — using biomformat reader")
    as.matrix(biomformat::biom_data(biomformat::read_biom(path)))
  }
}

biom_path            <- snakemake@input[["feature_table"]]
metadata_path        <- snakemake@input[["metadata"]]
taxonomy_path        <- snakemake@input[["taxonomy"]]
out_path             <- snakemake@output[["output_file"]]
res_path             <- snakemake@output[["res_file"]]
group_col            <- snakemake@params[["group_col"]]           %||% ""
db                   <- snakemake@params[["db"]]                  %||% ""
min_sample_presence  <- snakemake@params[["min_sample_presence"]] %||% 3
prv_cut              <- snakemake@params[["prv_cut"]]             %||% 0.28
lib_cut              <- snakemake@params[["lib_cut"]]             %||% 5000
p_adj_method         <- snakemake@params[["p_adj_method"]]        %||% "BH"
alpha                <- snakemake@params[["alpha"]]               %||% 0.05
pseudo_sens          <- snakemake@params[["pseudo_sens"]]         %||% TRUE
struc_zero           <- snakemake@params[["struc_zero"]]          %||% TRUE
neg_lb               <- snakemake@params[["neg_lb"]]              %||% TRUE
precollapsed         <- isTRUE(snakemake@params[["precollapsed"]])

dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

# canonical_lineage()/write_res(): the shared LEfSe .res interface, so the
# LDA bar chart and cladogram can be drawn for ANCOM-BC2 by the same renderer.
snakemake@source("write_res.R")

# Every early exit has to leave BOTH declared outputs behind or snakemake fails
# the rule. An empty .res is a valid "nothing to draw" signal: the renderer
# writes empty figures for it rather than erroring.
write_empty_outputs <- function() {
  write.table(data.frame(FeatureID = character(0), Taxon = character(0)),
              file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  dir.create(dirname(res_path), recursive = TRUE, showWarnings = FALSE)
  file.create(res_path)
}

message("[ANCOMBC2] Reading BIOM table: ", biom_path)
otu <- .read_biom_auto(biom_path)

if (!is.matrix(otu) || nrow(otu) == 0 || ncol(otu) == 0) {
  stop("BIOM table is empty or unreadable.")
}

message("[ANCOMBC2] Reading metadata: ", metadata_path)
meta <- read.delim(
  metadata_path,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

if (!(group_col %in% colnames(meta))) {
  stop("Grouping column '", group_col, "' not found in metadata.")
}

message("[ANCOMBC2] Reading taxonomy: ", taxonomy_path)
tax_df <- read.delim(
  taxonomy_path,
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE,
  stringsAsFactors = FALSE,
  quote = "",
  comment.char = ""
)

common_samples <- intersect(colnames(otu), rownames(meta))
if (length(common_samples) < 2) {
  stop("Too few overlapping samples between BIOM table and metadata.")
}

otu <- otu[, common_samples, drop = FALSE]
meta <- meta[common_samples, , drop = FALSE]
meta[[group_col]] <- factor(meta[[group_col]],
                            levels = natural_level_order(unique(as.character(meta[[group_col]]))))

# ---------- Model-safe column aliases ----------
# ancombc2() resolves variable names by parsing `fix_formula` as a string, so a
# metadata column named "Tumor Volume" cannot be referenced at all -- backticks
# do not help, because it does not strip them:
#   fix_formula = "Tumor Volume"   -> "variables specified are not in the meta data"
#   fix_formula = "`Tumor Volume`" -> the same error
# Nothing else in this pipeline minds such a name: snakemake round-trips it
# through wildcards, and no shell rule ever has a grouping column in its paths.
# So instead of forbidding
# it, give the model a spaceless ALIAS column and read the term names off that.
#
# Group VALUES are deliberately untouched. ancombc2 keeps the space in its term
# names ("lfc_Tumor_VolumeLow tumor volume"), and write_res() wants it in the
# class column. Only the COLUMN NAME is rewritten, and only for the fit.
#
# When a name is already spaceless the alias equals it and nothing is added, so
# this is a no-op for every study that does not need it.
safe_name <- function(x) gsub("[^A-Za-z0-9_]+", "_", x)

model_cols <- group_col[group_col %in% colnames(meta)]

alias_of <- setNames(as.list(safe_name(model_cols)), model_cols)
for (col in names(alias_of)) {
  a <- alias_of[[col]]
  # Sanitizing can land on a name some other column already uses; suffix until
  # it does not, so the alias never shadows real metadata.
  while (!identical(a, col) && a %in% colnames(meta)) a <- paste0(a, "_")
  alias_of[[col]] <- a
  if (!identical(a, col)) {
    meta[[a]] <- meta[[col]]
    message("[ANCOMBC2] Model alias for column '", col, "': '", a, "'")
  }
}
alias_or_self <- function(x) {
  vapply(as.character(x),
         function(v) if (!is.null(alias_of[[v]])) alias_of[[v]] else v,
         character(1), USE.NAMES = FALSE)
}

keep_taxa <- rowSums(otu, na.rm = TRUE) > 0
otu <- otu[keep_taxa, , drop = FALSE]

if (nrow(otu) == 0) {
  stop("No non-zero taxa remain after filtering.")
}

# ── Per-database taxonomic resolution ──────────────────────────────────────
# Same policy, and the same rationale, as scripts/run_lefse_Huttenhower.py:
# resolution is a property of the database, not of the method. Greengenes2's
# full-length classifier resolves real species; Silva138's species field is
# empty for every ASV, so collapsing past genus would only add noise.
#
# ANCOM-BC2 previously ran on the raw ASV table while LEfSe collapsed, which
# made the two methods incomparable and, worse, unrepresentable in the shared
# `.res` renderer: the `.res` is a taxonomy tree, so ASVs sharing a taxonomy
# string collapse onto one node no matter how many were significant.
COLLAPSE_LEVEL <- c(Silva138 = 6, Greengenes2 = 7)

collapse_to_level <- function(counts, taxonomy_df, level) {
  if (is.na(level) || level <= 0) return(counts)
  tax_col <- intersect(c("Taxon", "taxonomy"), colnames(taxonomy_df))[1]
  if (is.na(tax_col)) {
    warning("No taxonomy column found; leaving the table uncollapsed.")
    return(counts)
  }
  tax_vec <- setNames(as.character(taxonomy_df[[tax_col]]), trimws(rownames(taxonomy_df)))
  paths <- vapply(rownames(counts), function(id) {
    s <- unname(tax_vec[trimws(id)])
    if (length(s) == 0 || is.na(s)) return(NA_character_)
    parts <- trimws(strsplit(s, ";")[[1]])
    parts <- parts[nzchar(parts)]
    if (length(parts) == 0) return(NA_character_)
    paste(parts[seq_len(min(level, length(parts)))], collapse = ";")
  }, character(1))
  drop <- is.na(paths)
  if (any(drop))
    message("[ANCOMBC2] ", sum(drop), " feature(s) had no taxonomy and were dropped")
  rowsum(counts[!drop, , drop = FALSE], group = paths[!drop], reorder = FALSE)
}

collapse_level <- if (db %in% names(COLLAPSE_LEVEL)) COLLAPSE_LEVEL[[db]] else NA_integer_
if (precollapsed) {
  message("[ANCOMBC2] Using pre-collapsed genus table for ", db)
  tax_df <- data.frame(Taxon = rownames(otu), row.names = rownames(otu),
                       stringsAsFactors = FALSE, check.names = FALSE)
} else {
  if (is.na(collapse_level))
    warning("No ANCOMBC2 collapse level for db='", db, "'; running at feature level. ",
            "Known: ", paste(names(COLLAPSE_LEVEL), collapse = ", "))
  n_before <- nrow(otu)
  otu <- collapse_to_level(otu, tax_df, collapse_level)
  message("[ANCOMBC2] ", db, ": collapsed ", n_before, " feature(s) to ", nrow(otu),
          " at level ", collapse_level,
          " (features are taxonomy paths, so they map 1:1 onto .res tree nodes)")

  # After collapsing, a feature id IS its taxonomy path. Re-key the taxonomy
  # frame on those ids so every downstream lookup keeps working unchanged.
  if (!is.na(collapse_level) && collapse_level > 0) {
    tax_df <- data.frame(Taxon = rownames(otu), row.names = rownames(otu),
                         stringsAsFactors = FALSE, check.names = FALSE)
  }
}

tax_mat <- extract_taxonomy_matrix(tax_df, rownames(otu))

ps <- phyloseq::phyloseq(
  phyloseq::otu_table(otu, taxa_are_rows = TRUE),
  phyloseq::sample_data(meta),
  phyloseq::tax_table(tax_mat)
)

ps <- phyloseq::filter_taxa(ps, function(x) sum(x > 0) >= min_sample_presence, TRUE)

# Remove taxa with zero variance within any group level — ancombc2 cannot fit these
otu_mat <- as.matrix(phyloseq::otu_table(ps))
if (!phyloseq::taxa_are_rows(ps)) otu_mat <- t(otu_mat)
sample_groups <- meta[colnames(otu_mat), group_col]
zero_var_mask <- apply(otu_mat, 1, function(x) {
  group_vars <- tapply(x, sample_groups, var)
  # Keep a taxon if at least one group has non-zero variance;
  # only drop if ALL groups have zero variance (truly constant everywhere)
  all(group_vars == 0, na.rm = TRUE)
})
if (any(zero_var_mask)) {
  message("[ANCOMBC2] Removing ", sum(zero_var_mask),
          " taxa with zero within-group variance")
  ps <- phyloseq::prune_taxa(!zero_var_mask, ps)
}

if (phyloseq::ntaxa(ps) == 0) {
  message("[ANCOMBC2] No taxa remain after filtering. Writing empty results: ", out_path)
  write_empty_outputs()
  quit(save = "no", status = 0)
}

ref_level <- levels(meta[[group_col]])[1]
message("[ANCOMBC2] Reference (baseline) level for '", group_col, "': ", ref_level)

# The pipeline fits the flattest model there is: the grouping column alone, as a
# flat categorical, one run per DESIGN.FACTORS column. Nothing else joins the
# formula, ever.
#
# Interactions, covariates and random effects are deliberately absent rather than
# missing. Whether a column is a coefficient or a nuisance confounder is a
# property of the study, not of the pipeline, and it changes faster than a shared
# config schema can track. A study needing any of them writes its own script under
# main/<study>/additional_comparison_scripts/, where the choice is visible next to
# the reasoning for it.
#
# Interactions are additionally incompatible with `group=` mode by construction:
# ancombc2() selects coefficients with grepl(group, colnames(model_matrix)), which
# catches "DietIFA:RadiationIR" columns and differences them against main effects,
# yielding contrasts that mean nothing.
group_term  <- alias_or_self(group_col)
fix_formula <- group_term

run_ancombc2 <- function(ps_obj) {
  ANCOMBC::ancombc2(
    data       = ps_obj,
    fix_formula = fix_formula,
    group      = group_term,
    prv_cut    = prv_cut,
    lib_cut    = lib_cut,
    p_adj_method = p_adj_method,
    pseudo_sens = pseudo_sens,
    struc_zero = struc_zero,
    neg_lb     = neg_lb,
    alpha      = alpha,
    global     = TRUE,
    pairwise   = TRUE,
    dunnet     = FALSE,
    trend      = FALSE,
    verbose    = FALSE
  )
}

message("[ANCOMBC2] Running with formula: ", fix_formula)
ps_run <- ps
res <- NULL
repeat {
  res <- tryCatch(
    run_ancombc2(ps_run),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("Zero variances have been detected", msg)) {
        taxa_block <- gsub(".*Zero variances have been detected for the following taxa:\\s*", "", msg)
        taxa_block <- gsub("\\s*Please remove these taxa.*", "", taxa_block)
        bad_taxa   <- trimws(unlist(strsplit(taxa_block, "[,\n\\s]+")))
        bad_taxa   <- bad_taxa[nzchar(bad_taxa)]
        bad_taxa   <- bad_taxa[bad_taxa %in% phyloseq::taxa_names(ps_run)]
        if (length(bad_taxa) > 0) {
          message("[ANCOMBC2] Removing ", length(bad_taxa),
                  " zero-variance taxa and retrying")
          ps_run <<- phyloseq::prune_taxa(
            setdiff(phyloseq::taxa_names(ps_run), bad_taxa), ps_run)
          if (phyloseq::ntaxa(ps_run) == 0) {
            message("[ANCOMBC2] No taxa remain. Writing empty results: ", out_path)
            write_empty_outputs()
            quit(save = "no", status = 0)
          }
          return(NULL)  # signal: retry
        }
      }
      message("[ANCOMBC2] ancombc2 failed: ", msg)
      write_empty_outputs()
      quit(save = "no", status = 0)
    }
  )
  if (!is.null(res)) break
}

# 1. Extract the main results
res_df <- as.data.frame(res$res)

# 2. If pairwise was enabled, merge those columns in.
# ancombc2() returns the pairwise table as `res_pair`, not `pairwise` — the old
# key never existed, so this merge silently never ran and the mdFDR-corrected
# pairwise q-values were missing from the TSV.
if (!is.null(res$res_pair)) {
  pair_df <- as.data.frame(res$res_pair)
  # res_pair repeats every vs-reference contrast under the same column names as
  # res$res (with mdFDR-corrected q-values). Merging those would hand the plotter
  # duplicate ".x"/".y" copies of one contrast, so take only the contrasts that
  # res$res cannot express: the non-reference-vs-non-reference ones.
  pair_new <- setdiff(colnames(pair_df), setdiff(colnames(res_df), "taxon"))
  # Match by taxon to ensure rows align correctly
  res_df <- merge(res_df, pair_df[, pair_new, drop = FALSE], by = "taxon", all = TRUE)
}

# 3. Ensure we have a FeatureID column for the plotter
if ("taxon" %in% colnames(res_df)) {
  res_df$FeatureID <- res_df$taxon
} else {
  res_df$FeatureID <- rownames(res_df)
}

# 4. Map the Taxonomy (Handling the merge safely)
tax_col <- if ("Taxon" %in% colnames(tax_df)) "Taxon" else if ("taxonomy" %in% colnames(tax_df)) "taxonomy" else NULL

if (!is.null(tax_col)) {
  tax_map <- as.character(tax_df[[tax_col]])
  names(tax_map) <- trimws(rownames(tax_df))
  res_df$Taxon <- tax_map[trimws(as.character(res_df$FeatureID))]
} else {
  res_df$Taxon <- NA_character_
}

# 5. Attach baseline metadata so the plot script can label axes correctly
res_df$group_col    <- group_col
res_df$ref_level    <- ref_level

# 6. Dynamic Column Reordering (The Safety Check)
# Only include columns that actually exist in the data frame
first_cols <- c("FeatureID", "Taxon", "group_col", "ref_level")
existing_first <- first_cols[first_cols %in% colnames(res_df)]
other_cols <- setdiff(colnames(res_df), existing_first)

res_df <- res_df[, c(existing_first, other_cols), drop = FALSE]

write.table(
  res_df,
  file = out_path,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

message("[ANCOMBC2] Results written: ", out_path)

# ---------- Second rendering path: a LEfSe-style .res ----------
# The same LDA-bar-chart + cladogram the LEfSe methods emit, so a permissive
# LEfSe hit list and ANCOM-BC2's FDR-controlled one can be read side by side in
# one chart format. The tree is built from the FINAL retained phyloseq object,
# which is exactly the taxon universe res$res covers.
run_taxa   <- phyloseq::taxa_names(ps_run)
counts_run <- as.matrix(phyloseq::otu_table(ps_run))
if (!phyloseq::taxa_are_rows(ps_run)) counts_run <- t(counts_run)

sd_run    <- as(phyloseq::sample_data(ps_run), "data.frame")
# Read the values through the ALIAS column, not the real name: coercing
# sample_data to a data.frame puts every column name through make.names(), so
# "Tumor Volume" comes back as "Tumor.Volume" and a lookup by the real name
# silently yields NA. The alias is spaceless by construction, so make.names()
# leaves it alone -- and when the column needed no alias the two are identical.
group_run <- factor(as.character(sd_run[colnames(counts_run), group_term]),
                    levels = levels(meta[[group_col]]))

# Pruning taxa can leave a sample with nothing in it, which the .res writer
# refuses to normalize. Such a sample contributes no abundance to any node.
depths <- colSums(counts_run, na.rm = TRUE)
if (any(depths <= 0)) {
  message("[ANCOMBC2] Excluding ", sum(depths <= 0),
          " sample(s) with zero counts across retained taxa from the .res")
  counts_run <- counts_run[, depths > 0, drop = FALSE]
  group_run  <- droplevels(group_run[depths > 0])
}

tax_run     <- extract_taxonomy_matrix(tax_df, run_taxa)
lineage_map <- lapply(seq_len(nrow(tax_run)),
                      function(i) canonical_lineage(paste(tax_run[i, ], collapse = "|")))

mk <- build_ancombc2_markers(res, group_term, levels(meta[[group_col]]),
                             ref_level, alpha)
mk <- mk[mk$taxon %in% run_taxa, , drop = FALSE]
message("[ANCOMBC2] ", nrow(mk), " taxa assigned an enrichment class for the .res")

markers <- NULL
if (nrow(mk) > 0) {
  keys <- setNames(vapply(lineage_map, paste, character(1), collapse = "."), run_taxa)
  markers <- data.frame(
    lineage = unname(keys[mk$taxon]),
    class   = mk$class,
    lda     = mk$lda,
    pvalue  = mk$pvalue,
    stringsAsFactors = FALSE
  )
}
write_res(counts_run, group_run, lineage_map, markers, res_path)