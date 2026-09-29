#!/usr/bin/env Rscript
# Alpha and beta diversity for ahmad260430 under the Treatment x Radiation factorial.
#
# Companion to run_ancombc2_factorial.R. Differential abundance asks which taxa
# changed and pays a ~300-taxon multiple-testing tax to do it; these metrics ask
# whether the community changed at all, in one test per metric. That makes them
# far better powered for a diffuse effect spread thinly across many taxa —
# exactly the regime where the DA analysis found almost nothing for Radiation.
#
# The pipeline already computes every one of these, but only over `Group`, the
# flattened 6-level factor. This refits them as the 3x2 factorial.
#
# Day17 only: Day25 exists for the IFA arm alone, so including it aliases the
# timepoint difference into the Treatment terms. Confirmed that both databases'
# rarefaction depths are set by Day17 samples, so subsetting costs no depth.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/qiime2-2025.10-amplicon-core.sif \
#     Rscript run_diversity_factorial.R [--timepoint Day17] [--permutations 9999]

suppressPackageStartupMessages({ library(vegan) })

TREATMENT_LEVELS <- c("Control", "IFA", "Inosine")
RADIATION_LEVELS <- c("None", "IR")

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
CORE      <- file.path(STUDY_DIR, "qiime2_artifacts", "core-metrics-results")
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
# Results in tables/, plots in figures/, prose at the root. Routing by what the
# file is keeps the folder navigable however many parameterizations run.
TAB_DIR   <- file.path(OUT_DIR, "tables")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
timepoint <- arg_val("--timepoint", "Day17")
if (identical(timepoint, "all")) timepoint <- NULL
NPERM <- as.integer(arg_val("--permutations", "9999"))

# Metric registry. `db` NA means the artifact is database-agnostic, and that is
# verified rather than assumed: rules/diversity.smk feeds the non-phylogenetic
# metrics `table-analysis-rarefied.qza`, which carries no {db} wildcard, while
# the phylogenetic ones read `{db}-table-rarefied.qza` plus that database's tree.
#
# So the two databases differ for the phylogenetic metrics only -- but they
# differ a great deal more than "same data, different tree", and the difference
# is NOT symmetric:
#
#   Silva138    2321 ASV features, rarefied to 11008. Fragment insertion keeps
#               essentially everything (99.8% of features).
#   Greengenes2  182 backbone-genome features, rarefied to 5751. `qiime
#               greengenes2 non-v4-16s` MAPS ASVs onto reference genomes rather
#               than inserting them, collapsing 2325 ASVs onto 182 genomes and
#               dropping 40% of reads as unmappable.
#
# Consequences for reading the output: absolute Faith PD / UniFrac values are
# not comparable across databases (different feature space AND ~2x depth), so
# only the within-database test results transfer. Greengenes2 has less headroom
# for richness differences against a 182-unit ceiling, and its unmappable reads
# are not a random sample -- divergent and novel lineages are exactly what fails
# to map -- so its view is biased toward well-characterised taxa. Agreement
# across the two is therefore strong evidence; disagreement is ambiguous rather
# than informative.
alpha_metrics <- rbind(
  data.frame(metric = "shannon",  db = NA,            file = "shannon-vector.qza"),
  data.frame(metric = "simpson",  db = NA,            file = "simpson-vector.qza"),
  data.frame(metric = "evenness", db = NA,            file = "evenness-vector.qza"),
  data.frame(metric = "chao1",    db = NA,            file = "chao1-vector.qza"),
  data.frame(metric = "faith_pd", db = "Silva138",    file = "Silva138-faith-pd-vector.qza"),
  data.frame(metric = "faith_pd", db = "Greengenes2", file = "Greengenes2-faith-pd-vector.qza")
)
beta_metrics <- rbind(
  data.frame(metric = "jaccard",            db = NA,            file = "jaccard-distance-matrix.qza"),
  data.frame(metric = "bray_curtis",        db = NA,            file = "bray-curtis-distance-matrix.qza"),
  data.frame(metric = "unweighted_unifrac", db = "Silva138",    file = "Silva138-unweighted-unifrac-distance-matrix.qza"),
  data.frame(metric = "weighted_unifrac",   db = "Silva138",    file = "Silva138-weighted-unifrac-distance-matrix.qza"),
  data.frame(metric = "unweighted_unifrac", db = "Greengenes2", file = "Greengenes2-unweighted-unifrac-distance-matrix.qza"),
  data.frame(metric = "weighted_unifrac",   db = "Greengenes2", file = "Greengenes2-weighted-unifrac-distance-matrix.qza")
)

export_qza <- function(qza, want) {
  tmp <- file.path(tempdir(), basename(tools::file_path_sans_ext(qza)))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  ok <- system2("qiime", c("tools", "export", "--input-path", shQuote(qza),
                           "--output-path", shQuote(tmp)),
                stdout = FALSE, stderr = FALSE)
  f <- file.path(tmp, want)
  if (ok != 0 || !file.exists(f)) stop("export failed for ", qza)
  f
}

# ---------- metadata ----------
meta <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
if (!is.null(timepoint)) meta <- meta[meta$Timepoint == timepoint, , drop = FALSE]
meta$Treatment <- factor(meta$Treatment, levels = TREATMENT_LEVELS)
meta$Radiation <- factor(meta$Radiation, levels = RADIATION_LEVELS)
if (anyNA(meta$Treatment) || anyNA(meta$Radiation)) stop("Treatment/Radiation outside declared levels.")

message("[diversity] ", nrow(meta), " samples; design:")
print(table(meta$Treatment, meta$Radiation))
balanced <- length(unique(as.vector(table(meta$Treatment, meta$Radiation)))) == 1
message("[diversity] balanced: ", balanced,
        if (balanced) " (sequential and marginal SS coincide, term order is not a choice)" else
                      " -- UNBALANCED: adonis2 term order matters, results are sequential")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
results <- list()
# Tidy per-sample values for the plot script. betadisper distances in
# particular can only come from vegan, so they are emitted here rather than
# reimplemented downstream.
alpha_long <- list(); disper_long <- list()

# ---------- alpha ----------
for (i in seq_len(nrow(alpha_metrics))) {
  m <- alpha_metrics[i, ]
  qza <- file.path(CORE, m$file)
  if (!file.exists(qza)) { message("[alpha] missing, skipped: ", m$file); next }
  tab <- read.delim(export_qza(qza, "alpha-diversity.tsv"), row.names = 1, check.names = FALSE)
  common <- intersect(rownames(meta), rownames(tab))
  d <- meta[common, ]; d$value <- tab[common, 1]

  alpha_long[[length(alpha_long) + 1L]] <- data.frame(
    sampleid = rownames(d), metric = m$metric, db = ifelse(is.na(m$db), "any", m$db),
    Treatment = as.character(d$Treatment), Radiation = as.character(d$Radiation),
    value = d$value, stringsAsFactors = FALSE)

  fit <- lm(value ~ Treatment * Radiation, data = d)
  av  <- anova(fit)
  # Diagnostics: these tests are themselves weak at n=24, so treat them as
  # flags rather than gates.
  shapiro_p  <- tryCatch(shapiro.test(residuals(fit))$p.value, error = function(e) NA_real_)
  bartlett_p <- tryCatch(bartlett.test(value ~ interaction(Treatment, Radiation), data = d)$p.value,
                         error = function(e) NA_real_)
  for (term in c("Treatment", "Radiation", "Treatment:Radiation")) {
    results[[length(results) + 1L]] <- data.frame(
      analysis = "alpha", metric = m$metric, db = ifelse(is.na(m$db), "any", m$db),
      term = term, n = nrow(d),
      statistic = av[term, "F value"], effect = NA_real_,
      p_raw = av[term, "Pr(>F)"],
      resid_normal_p = shapiro_p, var_homog_p = bartlett_p,
      stringsAsFactors = FALSE)
  }
  message(sprintf("[alpha] %-9s %-12s n=%d  Treatment p=%.4f  Radiation p=%.4f  Interaction p=%.4f",
                  m$metric, ifelse(is.na(m$db), "-", m$db), nrow(d),
                  av["Treatment", "Pr(>F)"], av["Radiation", "Pr(>F)"], av["Treatment:Radiation", "Pr(>F)"]))
}

# ---------- beta ----------
for (i in seq_len(nrow(beta_metrics))) {
  m <- beta_metrics[i, ]
  qza <- file.path(CORE, m$file)
  if (!file.exists(qza)) { message("[beta] missing, skipped: ", m$file); next }
  dm <- as.matrix(read.delim(export_qza(qza, "distance-matrix.tsv"), row.names = 1, check.names = FALSE))
  common <- intersect(rownames(meta), rownames(dm))
  d <- meta[common, ]
  dist_sub <- as.dist(dm[common, common])

  set.seed(42)
  ad <- adonis2(dist_sub ~ Treatment * Radiation, data = d, permutations = NPERM, by = "terms")
  for (term in c("Treatment", "Radiation", "Treatment:Radiation")) {
    results[[length(results) + 1L]] <- data.frame(
      analysis = "beta_permanova", metric = m$metric, db = ifelse(is.na(m$db), "any", m$db),
      term = term, n = length(common),
      statistic = ad[term, "F"], effect = ad[term, "R2"], p_raw = ad[term, "Pr(>F)"],
      resid_normal_p = NA_real_, var_homog_p = NA_real_, stringsAsFactors = FALSE)
  }
  # PERMANOVA cannot separate a centroid shift from a spread difference, so a
  # significant term is uninterpretable without this.
  # "Cell" is the 6-way Treatment x Radiation combination: betadisper only takes one
  # grouping, so a PERMANOVA interaction driven by one unusually spread cell
  # would be invisible to the single-factor tests.
  groupings <- list(Treatment = d$Treatment, Radiation = d$Radiation,
                    Cell = droplevels(interaction(d$Treatment, d$Radiation, sep = "/")))
  for (fac in names(groupings)) {
    bd <- betadisper(dist_sub, groupings[[fac]])
    set.seed(42)
    pt <- permutest(bd, permutations = NPERM)
    results[[length(results) + 1L]] <- data.frame(
      analysis = "beta_dispersion", metric = m$metric, db = ifelse(is.na(m$db), "any", m$db),
      term = fac, n = length(common),
      statistic = pt$tab[1, "F"], effect = NA_real_, p_raw = pt$tab[1, "Pr(>F)"],
      resid_normal_p = NA_real_, var_homog_p = NA_real_, stringsAsFactors = FALSE)
    disper_long[[length(disper_long) + 1L]] <- data.frame(
      sampleid = names(bd$distances), metric = m$metric,
      db = ifelse(is.na(m$db), "any", m$db), grouping = fac,
      group = as.character(groupings[[fac]]), distance = unname(bd$distances),
      stringsAsFactors = FALSE)
    if (fac == "Cell") z_cell <- bd$distances
  }
  # The Cell test says whether any cell is unusually spread, not which factor the
  # spread follows. betadisper takes no formula, but PERMDISP is just an ANOVA on
  # distance-to-centroid, so fit the factorial model on the Cell distances and
  # test each term by permuting distances across samples.
  zF <- function(z) {
    zd <- data.frame(z = z, Treatment = d$Treatment, Radiation = d$Radiation)
    anova(lm(z ~ Treatment * Radiation, data = zd))[1:3, "F value"]
  }
  z_obs <- unname(z_cell)
  obs <- zF(z_obs)
  set.seed(42)
  perm <- replicate(NPERM, zF(sample(z_obs)))
  zp <- (rowSums(perm >= obs) + 1) / (NPERM + 1)
  for (j in 1:3) {
    results[[length(results) + 1L]] <- data.frame(
      analysis = "beta_dispersion_factorial", metric = m$metric,
      db = ifelse(is.na(m$db), "any", m$db),
      term = c("Treatment", "Radiation", "Treatment:Radiation")[j], n = length(common),
      statistic = obs[j], effect = NA_real_, p_raw = zp[j],
      resid_normal_p = NA_real_, var_homog_p = NA_real_, stringsAsFactors = FALSE)
  }
  message(sprintf("[beta]  %-18s %-12s  Treatment p=%.4f (R2=%.3f)  Radiation p=%.4f (R2=%.3f)  Interaction p=%.4f",
                  m$metric, ifelse(is.na(m$db), "-", m$db),
                  ad["Treatment", "Pr(>F)"], ad["Treatment", "R2"],
                  ad["Radiation", "Pr(>F)"], ad["Radiation", "R2"],
                  ad["Treatment:Radiation", "Pr(>F)"]))
}

res <- do.call(rbind, results)

# BH within each analysis x term family, since the metrics inside one family are
# alternative views of the same samples rather than independent questions.
# Silva138 is its own family: it is the secondary database, reported apart from
# the main results, and its UniFrac/Faith PD re-ask the Greengenes2 questions on
# another tree, so pooling them would charge the main results for tests they
# never make.
res$p_BH <- NA_real_
fam <- paste(res$analysis, res$term, res$db == "Silva138")
for (f in unique(fam)) {
  k <- fam == f
  res$p_BH[k] <- p.adjust(res$p_raw[k], method = "BH")
}

suffix  <- if (is.null(timepoint)) "_all-timepoints" else paste0("_", timepoint)
outfile <- file.path(TAB_DIR, paste0("diversity_factorial", suffix, ".tsv"))
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("\n[diversity] Wrote ", outfile)

for (nm in c("alpha_values", "dispersion_distances")) {
  dat <- if (nm == "alpha_values") do.call(rbind, alpha_long) else do.call(rbind, disper_long)
  if (is.null(dat)) next
  f <- file.path(TAB_DIR, paste0("diversity_", nm, suffix, ".tsv"))
  write.table(dat, f, sep = "\t", quote = FALSE, row.names = FALSE)
  message("[diversity] Wrote ", f)
}

message("\n[diversity] Significant after BH (q < 0.05):")
sig <- res[!is.na(res$p_BH) & res$p_BH < 0.05, ]
if (nrow(sig) == 0) message("    none") else
  for (j in seq_len(nrow(sig))) message(sprintf("    %-16s %-18s %-12s %-14s p=%.4g  q=%.4g",
    sig$analysis[j], sig$metric[j], sig$db[j], sig$term[j], sig$p_raw[j], sig$p_BH[j]))
