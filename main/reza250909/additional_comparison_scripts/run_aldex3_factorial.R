#!/usr/bin/env Rscript
# Fit ALDEx3 as a Diet x Radiation factorial for reza250909, the secondary
# check on the ANCOM-BC2 fit (run_ancombc2_factorial.R). Greengenes2 only.
#
# Needs ALDEx3 1.3.1: run it in containers/ancombc2-env.sif.
#
# Same taxa as the ANCOM-BC2 default fit: the counts are read, collapsed and
# filtered by the same rules (every Diet x Radiation cell sees the taxon; >= 4
# of 32 samples), and the script stops if the resulting set differs from the
# taxa in ANCOMBC2_factorial_Greengenes2.tsv.
#
# Model: aldex(Y, ~ Diet * Radiation, method = "lm", test = "t.HC3") with the
# CLR scale model (clr.sm). gamma is the SD (log2) of the random scale shift
# clr.sm adds to every model coefficient in each Monte Carlo draw; gamma = 0 is
# the fixed CLR assumption. Estimates are log2.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_aldex3_factorial.R [--gamma 0.5] [--seed 42] [--nsample 1000]

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Radiation-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
META_PATH <- file.path(OUT_DIR, "metadata_factorial.tsv")   # build_factorial_metadata.py

LOG_PREFIX <- "[aldex3]"
source(file.path(THIS_DIR, "ancombc2_common.R"))

suppressPackageStartupMessages(library(ALDEx3))
if (packageVersion("ALDEx3") != "1.3.1")
  stop("ALDEx3 ", packageVersion("ALDEx3"), " is installed; this script is pinned to 1.3.1. ",
       "Run it in containers/ancombc2-env.sif.")
.say("ALDEx3 ", packageVersion("ALDEx3"))

DB                   <- "Greengenes2"
P_ADJ_METHOD         <- "BH"
DEFAULT_MIN_PRESENCE <- 4L
DIET_LEVELS          <- c("Control", "SCFAs")
RADIATION_LEVELS     <- c("None", "IR")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
GAMMA   <- as.numeric(arg_val("--gamma", "0.5"))
SEED    <- as.integer(arg_val("--seed", "42"))
NSAMPLE <- as.integer(arg_val("--nsample", "1000"))
if (is.na(GAMMA) || GAMMA < 0) stop("--gamma must be a non-negative number")
if (is.na(SEED) || is.na(NSAMPLE) || NSAMPLE < 1) stop("--seed and --nsample must be integers")

# Stability reruns (another seed or nsample) get their own file.
suffix <- paste0("_gamma", format(GAMMA))
if (SEED != 42L) suffix <- paste0(suffix, "_seed", SEED)
if (NSAMPLE != 1000L) suffix <- paste0(suffix, "_ns", NSAMPLE)

# ---------- Counts and taxa, exactly as run_ancombc2_factorial.R ----------

otu_all <- read_biom_auto(file.path(STUDY_DIR, "tables", "study-seqs.biom"))
meta <- read.delim(META_PATH, header = TRUE, row.names = 1, sep = "\t",
                   check.names = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
common <- intersect(colnames(otu_all), rownames(meta))
otu_all <- otu_all[, common, drop = FALSE]
meta    <- meta[common, , drop = FALSE]
meta$Diet      <- factor(meta$Diet, levels = DIET_LEVELS)
meta$Radiation <- factor(meta$Radiation, levels = RADIATION_LEVELS)
if (anyNA(meta$Diet) || anyNA(meta$Radiation)) stop("Diet/Radiation outside the declared levels")

tax_map <- read_taxonomy_map(file.path(STUDY_DIR, "tables", "exported-taxonomy", paste0(DB, "_taxonomy.tsv")))
counts <- collapse_to_level(otu_all[rowSums(otu_all) > 0, , drop = FALSE], tax_map, COLLAPSE_LEVEL)

cell <- interaction(meta$Diet, meta$Radiation, drop = TRUE)
in_every_cell <- sapply(levels(cell), function(k) rowSums(counts[, cell == k, drop = FALSE] > 0) > 0)
keep <- rowSums(!in_every_cell) == 0 & rowSums(counts > 0) >= DEFAULT_MIN_PRESENCE
counts <- counts[keep, , drop = FALSE]

ancombc_taxa <- read.delim(file.path(TAB_DIR, paste0("ANCOMBC2_factorial_", DB, ".tsv")),
                           check.names = FALSE, quote = "")$taxon
if (!setequal(rownames(counts), ancombc_taxa))
  stop("taxon set differs from the ANCOM-BC2 default fit: ",
       length(setdiff(rownames(counts), ancombc_taxa)), " only here, ",
       length(setdiff(ancombc_taxa, rownames(counts))), " only in ANCOM-BC2")
.say(sprintf("%s: %d taxa, identical to the ANCOM-BC2 default fit", DB, nrow(counts)))

# ---------- Fit ----------

.say(sprintf("Fitting ~ Diet * Radiation, clr.sm gamma = %s, nsample = %d, seed = %d",
             format(GAMMA), NSAMPLE, SEED))
set.seed(SEED)
fit <- aldex(counts, ~ Diet * Radiation, data = meta, method = "lm",
             nsample = NSAMPLE, scale = clr.sm, gamma = GAMMA,
             test = "t.HC3", p.adjust.method = P_ADJ_METHOD, n.cores = 1,
             return.pars = c("X", "estimate", "std.error", "p.val", "p.val.adj"))

# estimate/std.error are coefficient x taxon x draw; report the mean over draws.
coefs <- setdiff(rownames(fit$p.val), "(Intercept)")
out <- data.frame(taxon = colnames(fit$p.val), check.names = FALSE)
# ANCOM-BC2 drops a zero as a missing log value; ALDEx3 keeps it as a Dirichlet
# draw. How many mice have a zero is what separates the two methods' answers.
out$n_zero_mice <- rowSums(counts[out$taxon, , drop = FALSE] == 0)
for (co in coefs) {
  out[[paste0("est_", co)]] <- rowMeans(fit$estimate[co, , , drop = TRUE])
  out[[paste0("se_", co)]]  <- rowMeans(fit$std.error[co, , , drop = TRUE])
  out[[paste0("p_", co)]]   <- fit$p.val[co, ]
  out[[paste0("q_", co)]]   <- fit$p.val.adj[co, ]
}
out_path <- file.path(TAB_DIR, paste0("ALDEx3_factorial_", DB, suffix, ".tsv"))
write.table(out, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
.say("Wrote ", out_path)
for (co in coefs)
  .say(sprintf("  %-24s q < 0.05: %d of %d", co, sum(out[[paste0("q_", co)]] < 0.05), nrow(out)))
