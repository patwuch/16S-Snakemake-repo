#!/usr/bin/env Rscript
# ANCOM-BC2 for the Day17-vs-Day25 contrast, held within one treatment arm.
#
# Why this is a project-local script and not the pipeline. The project-agnostic
# Snakefile runs only the most flattened ANCOM-BC2: one DESIGN.FACTORS column,
# flat categorical levels, `group = <column>`. Declaring a covariate or an
# interaction there is now a hard config error, because whether a column is a
# coefficient or a nuisance confounder is a property of the study and it shifts
# faster than a shared schema can track. This fit needs a random effect, so it
# belongs here, next to the reasoning for it.
#
# Why a random effect. The two timepoints are the SAME four mice sampled twice,
# so the eight observations per arm are four pairs, not eight independent draws.
# `rand_formula = "(1|MouseID)"` says so. Ignoring it treats within-mouse
# variation as between-mouse variation, which inflates the residual and throws
# away the one thing this design has going for it at n = 4.
#
# Why `group = NULL`. Timepoint has two levels, so the multi-group machinery
# (global omnibus, mdFDR pairwise, structural zeros) has nothing to add over the
# single coefficient, and mdFDR would mislabel the correction on the figure axis.
# Plain Benjamini-Hochberg across features, as the axis says.
#
# Scope. This is a paired within-arm comparison of two timepoints eight days
# apart, with treatment held fixed. The `colon` suffix on the Day25 sample ids is
# a labelling convention, not a body site; both timepoints were collected from
# the same place. See comparison-day17-day25/report.html.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_ancombc2_day17_day25.R [--db Silva138|Greengenes2|both] [--arm IFA|IFA_IR]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
REPO_ROOT <- dirname(dirname(STUDY_DIR))

LOG_PREFIX <- "[day17_day25]"
source(file.path(THIS_DIR, "ancombc2_common.R"))
# canonical_lineage()/write_res(): the shared .res interface the effect-size
# barplot and cladogram renderer consumes, so these coefficients draw in the same
# chart format as the LEfSe results beside them.
source(file.path(REPO_ROOT, "scripts", "write_res.R"))

suppressPackageStartupMessages({
  library(biomformat)
  library(ANCOMBC)
})
require_ancombc("2.14.0")

# --- Settings. Mirror config/config_ahmad260430.yaml's ANCOMBC2_* block. ---
MIN_SAMPLE_PRESENCE <- 2
PRV_CUT             <- 0.05
LIB_CUT             <- 1000
P_ADJ_METHOD        <- "BH"
ALPHA               <- 0.05
PSEUDO_SENS         <- TRUE

# Slug -> the Group value it selects. Only these two arms have both timepoints;
# the other four are Day17 only. Matches TREATMENTS in run_comparison.py.
ARMS <- c(IFA = "IFA", IFA_IR = "IFA+IR")

DAY_LEVELS <- c("Day17", "Day25")   # Day17 is the reference, so a positive
                                    # coefficient means enriched at Day25

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
db_arg <- arg_val("--db", "both")
dbs <- if (db_arg == "both") c("Silva138", "Greengenes2") else db_arg
arm_arg <- arg_val("--arm", NULL)
slugs <- if (is.null(arm_arg)) names(ARMS) else arm_arg
if (!all(slugs %in% names(ARMS)))
  stop("--arm must be one of: ", paste(names(ARMS), collapse = ", "))

# metadata.tsv has no MouseID column, so the pairing is recovered from the
# sample id. `E1_1.00` (Day17) and `1.00colon` (Day25) are both mouse 1_00. The
# September 10 comparison folders carry an explicit MouseID column that this
# reproduces; it was dropped when run_comparison.py took over.
mouse_id <- function(sid) {
  m <- regmatches(sid, regexec("^E\\d+_(\\d+)\\.(\\d+)$", sid))[[1]]
  if (length(m) == 0) m <- regmatches(sid, regexec("^(\\d+)\\.(\\d+)colon$", sid))[[1]]
  if (length(m) == 0) stop("Cannot parse a mouse id out of sample '", sid, "'")
  paste(m[2], m[3], sep = "_")
}

biom_path <- file.path(STUDY_DIR, "tables", "study-seqs.biom")
meta_path <- file.path(STUDY_DIR, "metadata.tsv")

.say("Reading BIOM table: ", biom_path)
otu_all <- read_biom_auto(biom_path)
meta_all <- read.delim(meta_path, header = TRUE, row.names = 1, sep = "\t",
                       check.names = FALSE, stringsAsFactors = FALSE,
                       quote = "", comment.char = "")
for (col in c("Group", "Timepoint"))
  if (!(col %in% colnames(meta_all))) stop("metadata.tsv has no '", col, "' column.")

common <- intersect(colnames(otu_all), rownames(meta_all))
otu_all  <- otu_all[, common, drop = FALSE]
meta_all <- meta_all[common, , drop = FALSE]

for (slug in slugs) {
  arm_name <- ARMS[[slug]]
  out_dir  <- file.path(STUDY_DIR, "comparison-day17-day25")

  meta <- meta_all[meta_all$Group == arm_name, , drop = FALSE]
  if (nrow(meta) == 0) stop("No samples with Group == '", arm_name, "'")
  meta$Timepoint <- factor(meta$Timepoint, levels = DAY_LEVELS)
  if (anyNA(meta$Timepoint))
    stop("Arm '", arm_name, "' has Timepoint values outside ", paste(DAY_LEVELS, collapse = "/"))
  meta$MouseID <- factor(vapply(rownames(meta), mouse_id, character(1)))

  # A random effect only buys anything if the pairing is complete. An arm where
  # some mouse appears once contributes a singleton group, which the LME fits as
  # residual and which silently costs the power the pairing was meant to add.
  pairing <- table(meta$MouseID, meta$Timepoint)
  if (!all(pairing == 1))
    stop("Arm '", arm_name, "' is not cleanly paired; every mouse must appear ",
         "exactly once per timepoint:\n",
         paste(capture.output(print(pairing)), collapse = "\n"))

  .say(slug, " (Group == '", arm_name, "'): ", nrow(meta), " samples, ",
       nlevels(meta$MouseID), " mice x ", nlevels(meta$Timepoint), " timepoints")

  counts_arm <- otu_all[, rownames(meta), drop = FALSE]
  counts_arm <- counts_arm[rowSums(counts_arm, na.rm = TRUE) > 0, , drop = FALSE]

  for (db in dbs) {
    tax_path <- file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(db, "_taxonomy.tsv"))
    if (!file.exists(tax_path)) { .say("Skipping ", db, ": no taxonomy at ", tax_path); next }
    tax_map <- read_taxonomy_map(tax_path)
    level <- if (db %in% names(COLLAPSE_LEVEL)) COLLAPSE_LEVEL[[db]] else 0L

    counts_db <- collapse_to_level(counts_arm, tax_map, level)
    keep_db <- apply(counts_db, 1, function(x) sum(x > 0) >= MIN_SAMPLE_PRESENCE)
    counts_db <- counts_db[keep_db, , drop = FALSE]
    .say(sprintf("%s/%s: collapse level %d -> %d feature(s), %d retained (>= %d samples)",
                 slug, db, level, length(keep_db), nrow(counts_db), MIN_SAMPLE_PRESENCE))
    if (nrow(counts_db) == 0) { .say("nothing left for ", slug, "/", db); next }

    fit <- function(rand) {
      ancombc2(data = counts_db, taxa_are_rows = TRUE, meta_data = meta,
               fix_formula = "Timepoint", rand_formula = rand,
               group = NULL, prv_cut = PRV_CUT, lib_cut = LIB_CUT,
               p_adj_method = P_ADJ_METHOD, alpha = ALPHA,
               pseudo_sens = PSEUDO_SENS, struc_zero = FALSE, neg_lb = FALSE,
               global = FALSE, pairwise = FALSE, dunnet = FALSE, trend = FALSE,
               verbose = FALSE)
    }

    # The paired fit is the result. The unpaired one is written beside it so the
    # report can say what the pairing bought instead of asserting that it helped.
    .say("Fitting paired: Timepoint + (1|MouseID), group = NULL")
    res_paired <- tryCatch(fit("(1|MouseID)"),
      error = function(e) stop("paired ancombc2 failed for ", slug, "/", db, ": ",
                               conditionMessage(e)))
    .say("Fitting unpaired: Timepoint, group = NULL")
    res_unpaired <- tryCatch(fit(NULL), error = function(e) {
      .say("WARNING: unpaired fit failed (", conditionMessage(e), "); paired result stands alone")
      NULL
    })

    tab_dir <- file.path(out_dir, "tables")
    dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
    # One flat folder for the whole comparison, self-describing filenames, the
    # same grammar run_comparison.py's art() builds and the same shape
    # comparison-Treatment-x-Radiation-factorial uses. The arm is a filename token
    # because the two arms are one comparison asked twice.
    artefact <- function(...) file.path(
      tab_dir, paste0(paste(c("ANCOMBC2", "day17day25", slug, db, ...), collapse = "_")))

    write_table_for <- function(res_obj, path) {
      out <- as.data.frame(res_obj$res)
      # Collapsed features ARE taxonomy paths, so the id is already the label.
      out$Taxon <- if (level > 0) as.character(out$taxon) else
                   unname(tax_map[trimws(as.character(out$taxon))])
      out <- out[, c("taxon", "Taxon", setdiff(colnames(out), c("taxon", "Taxon")))]
      write.table(out, path, sep = "\t", quote = FALSE, row.names = FALSE)
      .say("Wrote ", path)
      out
    }

    out_paired <- write_table_for(
      res_paired, paste0(artefact(), ".tsv"))

    # Per-sample log abundance with the estimated sampling fraction removed.
    # Relative abundance is compositional, so one taxon collapsing inflates every
    # other one mechanically; this is the quantity the model actually compares,
    # and the only per-sample one that is not subject to that. The `res` table's
    # lfc is a single number per taxon, so it cannot carry a per-mouse view.
    bc <- as.data.frame(res_paired$bias_correct_log_table)
    bc <- cbind(taxon = rownames(bc), bc)
    bc_path <- paste0(artefact("biascorrected"), ".tsv")
    write.table(bc, bc_path, sep = "\t", quote = FALSE, row.names = FALSE)
    .say("Wrote ", bc_path)
    if (!is.null(res_unpaired))
      write_table_for(res_unpaired, paste0(artefact("unpaired"), ".tsv"))

    # ---- one .res for the single coefficient, for the shared renderer ----
    coef    <- paste0("Timepoint", DAY_LEVELS[2])
    res_df  <- as.data.frame(res_paired$res)
    qc <- paste0("q_", coef); lc <- paste0("lfc_", coef); sc <- paste0("passed_ss_", coef)
    if (!all(c(qc, lc) %in% colnames(res_df)))
      stop("Coefficient '", coef, "' absent from the fit for ", slug, "/", db)

    res_taxa   <- as.character(res_df$taxon)
    counts_res <- counts_db[res_taxa[res_taxa %in% rownames(counts_db)], , drop = FALSE]
    res_taxa   <- rownames(counts_res)
    lineage_map <- lapply(res_taxa, function(t) {
      s <- if (level > 0) t else unname(tax_map[trimws(t)])
      if (length(s) == 0 || is.na(s)) s <- ""
      canonical_lineage(paste(trimws(strsplit(s, ";")[[1]]), collapse = "|"))
    })
    keys   <- vapply(lineage_map, paste, character(1), collapse = ".")
    row_of <- match(res_taxa, as.character(res_df$taxon))

    q  <- res_df[[qc]][row_of]; l <- res_df[[lc]][row_of]
    ss <- if (sc %in% colnames(res_df)) res_df[[sc]][row_of] else TRUE
    # passed_ss is ancombc2's own pseudo-count sensitivity check: FALSE means the
    # call flips with the pseudo-count, i.e. a zero-handling artifact rather than
    # signal. Every count reported here passes both it and q < alpha.
    sig <- !is.na(q) & !is.na(l) & q < ALPHA & (is.na(ss) | ss %in% c(TRUE, "TRUE"))

    # ancombc2 fits lm() on log(count) with R's log(), so lfc is a NATURAL-log
    # fold change. Divide by log(2) to publish the log2 fold change the
    # effect-size axis names. Ranking and q-values are unaffected; only scale.
    markers <- if (any(sig)) data.frame(
      lineage = keys[sig],
      class   = ifelse(l[sig] > 0, DAY_LEVELS[2], DAY_LEVELS[1]),
      lda     = abs(l[sig]) / log(2),
      pvalue  = q[sig], stringsAsFactors = FALSE) else NULL

    cs <- counts_res[, rownames(meta), drop = FALSE]
    cs <- cs[, colSums(cs, na.rm = TRUE) > 0, drop = FALSE]
    res_path <- paste0(artefact(), ".res")
    write_res(cs, as.character(meta[colnames(cs), "Timepoint"]), lineage_map, markers, res_path)

    n_raw <- sum(!is.na(q) & q < ALPHA)
    .say(sprintf("%s/%s: %d marker(s) at q < %g AND passed_ss (%d at raw q < %g), of %d tested",
                 slug, db, if (is.null(markers)) 0L else nrow(markers), ALPHA,
                 n_raw, ALPHA, nrow(res_df)))
    if (!is.null(res_unpaired)) {
      u <- as.data.frame(res_unpaired$res)
      us <- if (sc %in% colnames(u)) u[[sc]] %in% c(TRUE, "TRUE") else TRUE
      .say(sprintf("%s/%s: unpaired fit would give %d marker(s) on the same filter",
                   slug, db, sum(!is.na(u[[qc]]) & u[[qc]] < ALPHA & us)))
    }
  }
}

.say("Done.")
