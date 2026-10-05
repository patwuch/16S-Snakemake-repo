#!/usr/bin/env Rscript
# Figure inputs from run_maaslin3_rows.R's abundance results (>= 2 samples):
#   * one .res per two-group report row, for the shared barplot/cladogram
#     renderer. The tree covers the features the row's coefficient was fitted
#     for; markers are the robust calls (q < 0.05 in the full fit and in all 19
#     leave-one-out refits), with |log2 change in share| as the effect size.
#     Names follow the ANCOM-BC2 convention (_ref<Diet>, _condref<Condition>,
#     then the coefficient) so the plotting and structural-zero scripts read
#     them unchanged.
#   * <stem>_tested.tsv beside each .res: the fitted features, for
#     write_structural_zero_res.R.
#   * MaAsLin3_factorial_Greengenes2_int[_refInulin].tsv: the interaction
#     coefficients in the wide layout plot_interaction_forest.py reads, in log2.
#
# Base R only. Runs in containers/maaslin3-env.sif:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/maaslin3-env.sif Rscript write_maaslin3_res.R

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
REPO_ROOT <- dirname(dirname(dirname(THIS_DIR)))
source(file.path(THIS_DIR, "factorial_rows.R"))
source(file.path(REPO_ROOT, "scripts", "write_res.R"))
OUT_DIR <- file.path(dirname(THIS_DIR), "comparison-Diet-x-Condition-factorial")
MA_DIR  <- file.path(OUT_DIR, "tables", "maaslin3")
TAB_DIR <- file.path(OUT_DIR, "tables")
PRE <- "MaAsLin3_factorial_Greengenes2"

res <- read.delim(file.path(MA_DIR, "MaAsLin3_factorial_Greengenes2_min2_rows.tsv"), stringsAsFactors = FALSE)
res <- res[res$model == "abundance", ]
counts <- read.delim(file.path(MA_DIR, "counts_Greengenes2.tsv"), row.names = 1, check.names = FALSE)
meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), row.names = 1, stringsAsFactors = FALSE)[rownames(counts), ]
key <- read.delim(file.path(MA_DIR, "feature_key.tsv"), stringsAsFactors = FALSE)
taxon_of <- setNames(key$taxon, key$feature)
lineage_of <- function(t) canonical_lineage(paste(trimws(strsplit(t, ";")[[1]]), collapse = "|"))

group <- function(d, c) if (c == "Control") d else paste0(d, "+IR")
stem_of <- function(k) paste0(PRE, if (ROWS$c_ref[k] == "IR") "_condrefIR" else "",
                              if (ROWS$d_ref[k] != "Cellulose") paste0("_ref", ROWS$d_ref[k]) else "",
                              "_", ROWS$coef[k])

for (k in which(ROWS$family != "interaction")) {
  d <- res[res$row == ROWS$row[k] & res$stratum == ROWS$stratum[k] & !is.na(res$p), ]
  if (grepl("^Diet", ROWS$coef[k])) {
    pos <- group(sub("^Diet", "", ROWS$coef[k]), ROWS$c_ref[k]); neg <- group(ROWS$d_ref[k], ROWS$c_ref[k])
  } else {
    pos <- group(ROWS$d_ref[k], "IR"); neg <- group(ROWS$d_ref[k], "Control")
  }
  keep <- meta$Group %in% c(pos, neg)
  cs <- t(counts[keep, d$feature, drop = FALSE])
  rownames(cs) <- taxon_of[d$feature]
  cs <- cs[, colSums(cs) > 0, drop = FALSE]
  lineage_map <- lapply(rownames(cs), lineage_of)
  hit <- d$robust %in% c(TRUE, "TRUE")
  markers <- if (any(hit)) data.frame(
    lineage = vapply(taxon_of[d$feature[hit]], function(t) paste(lineage_of(t), collapse = "."), character(1)),
    class = ifelse(d$coef[hit] > 0, pos, neg), lda = abs(d$coef[hit]), pvalue = d$q[hit],
    stringsAsFactors = FALSE) else NULL
  stem <- stem_of(k)
  write_res(cs, meta$Group[keep][match(colnames(cs), rownames(meta)[keep])], lineage_map, markers,
            file.path(TAB_DIR, paste0(stem, ".res")))
  write.table(data.frame(taxon = unname(taxon_of[d$feature])), file.path(TAB_DIR, paste0(stem, "_tested.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("[maaslin3-res] %s: %d fitted, %d robust (%s vs %s)", stem, nrow(d), sum(hit), pos, neg))
}

# Interaction coefficients, wide, one table per fit, for the forest plot.
# W is the t statistic, from which the plot recovers each feature's df.
for (d_ref in c("Cellulose", "Inulin")) {
  ks <- which(ROWS$family == "interaction" & ROWS$d_ref == d_ref)
  wide <- Reduce(function(a, b) merge(a, b, by = "taxon", all = TRUE), lapply(ks, function(k) {
    d <- res[res$row == ROWS$row[k] & !is.na(res$p), ]
    cf <- ROWS$coef[k]
    out <- data.frame(taxon = unname(taxon_of[d$feature]), d$coef, d$stderr, d$coef / d$stderr, d$p, d$q,
                      d$robust %in% c(TRUE, "TRUE"), stringsAsFactors = FALSE)
    # passed_ss_ is the forest plot's "counts as a hit" column: here, robust.
    names(out)[-1] <- paste0(c("lfc_", "se_", "W_", "p_", "q_", "passed_ss_"), cf)
    out
  }))
  f <- file.path(TAB_DIR, paste0(PRE, "_int", if (d_ref == "Cellulose") "" else "_refInulin", ".tsv"))
  write.table(wide, f, sep = "\t", quote = FALSE, row.names = FALSE)
  message(sprintf("[maaslin3-res] %s: %d features", basename(f), nrow(wide)))
}
