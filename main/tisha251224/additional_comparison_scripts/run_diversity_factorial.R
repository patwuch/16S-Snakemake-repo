#!/usr/bin/env Rscript
# Alpha and beta diversity for tisha251224 under the Diet x Condition factorial.
# Greengenes2 only: Silva138's phylogenetic metrics are not used.
#
# Ported from ahmad260430/additional_comparison_scripts/run_diversity_factorial.R.
# The one change of substance: that design was balanced (4 per cell), this one
# is not (2-4 per cell), so sequential sums of squares would make each main
# effect depend on term order. Every test here is Type II instead: each main
# effect adjusted for the other, the interaction adjusted for both, all against
# the full model's residual. Computed without the car package (not in the
# container) by fitting both term orders and taking the main effect that enters
# second from each: its sequential SS is exactly the Type II SS, and the
# residual is the full model's either way. For PERMANOVA the same holds for
# adonis2(by = "terms"), whose permutations are free under the full model.
#
# Usage (from anywhere):
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/qiime2-2025.10-amplicon-core.sif \
#     Rscript run_diversity_factorial.R [--permutations 9999]

suppressPackageStartupMessages({ library(vegan) })

DIET_LEVELS      <- c("Cellulose", "Inulin", "RS")
CONDITION_LEVELS <- c("Control", "IR")
TERMS <- c("Diet", "Condition", "Diet:Condition")

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "",
              grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
CORE      <- file.path(STUDY_DIR, "qiime2_artifacts", "core-metrics-results")
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Diet-x-Condition-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
dir.create(TAB_DIR, recursive = TRUE, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default = NULL) {
  i <- match(flag, args); if (is.na(i) || i == length(args)) default else args[i + 1L]
}
NPERM <- as.integer(arg_val("--permutations", "9999"))

# Shannon, Simpson, evenness, Chao1, Jaccard and Bray-Curtis come from the
# database-agnostic ASV table (rules/diversity.smk feeds them
# table-analysis-rarefied.qza, rarefied to 11,479 here); Faith PD and UniFrac
# read the Greengenes2 table and tree, rarefied to 5,509.
alpha_metrics <- data.frame(
  metric = c("shannon", "simpson", "evenness", "chao1", "faith_pd"),
  db     = c("any", "any", "any", "any", "Greengenes2"),
  file   = c("shannon-vector.qza", "simpson-vector.qza", "evenness-vector.qza",
             "chao1-vector.qza", "Greengenes2-faith-pd-vector.qza"))
beta_metrics <- data.frame(
  metric = c("jaccard", "bray_curtis", "unweighted_unifrac", "weighted_unifrac"),
  db     = c("any", "any", "Greengenes2", "Greengenes2"),
  file   = c("jaccard-distance-matrix.qza", "bray-curtis-distance-matrix.qza",
             "Greengenes2-unweighted-unifrac-distance-matrix.qza",
             "Greengenes2-weighted-unifrac-distance-matrix.qza"))

export_qza <- function(qza, want) {
  tmp <- file.path(tempdir(), basename(tools::file_path_sans_ext(qza)))
  dir.create(tmp, recursive = TRUE, showWarnings = FALSE)
  ok <- system2("qiime", c("tools", "export", "--input-path", shQuote(qza),
                           "--output-path", shQuote(tmp)), stdout = FALSE, stderr = FALSE)
  f <- file.path(tmp, want)
  if (ok != 0 || !file.exists(f)) stop("export failed for ", qza)
  f
}

# Type II F table for a univariate response: rows Diet, Condition, Diet:Condition.
type2_lm <- function(y, d) {
  d$y <- y
  a1 <- anova(lm(y ~ Diet * Condition, data = d))   # Condition | Diet
  a2 <- anova(lm(y ~ Condition * Diet, data = d))   # Diet | Condition
  rbind(Diet = a2["Diet", ], Condition = a1["Condition", ], "Diet:Condition" = a1["Diet:Condition", ])
}

# ---------- metadata ----------
meta <- read.delim(file.path(OUT_DIR, "metadata_factorial.tsv"), header = TRUE, row.names = 1,
                   sep = "\t", check.names = FALSE, stringsAsFactors = FALSE,
                   quote = "", comment.char = "")
meta$Diet      <- factor(meta$Diet, levels = DIET_LEVELS)
meta$Condition <- factor(meta$Condition, levels = CONDITION_LEVELS)
if (anyNA(meta$Diet) || anyNA(meta$Condition)) stop("Diet/Condition outside declared levels.")
message("[diversity] ", nrow(meta), " samples; design:")
print(table(meta$Diet, meta$Condition))

results <- list(); alpha_long <- list(); disper_long <- list()
row_of <- function(...) data.frame(..., stringsAsFactors = FALSE)
check_n <- function(common, what)
  if (length(common) != nrow(meta)) message("[", what, "] WARNING: ", nrow(meta) - length(common),
                                            " sample(s) absent from the artifact")

# ---------- alpha ----------
for (i in seq_len(nrow(alpha_metrics))) {
  m <- alpha_metrics[i, ]
  tab <- read.delim(export_qza(file.path(CORE, m$file), "alpha-diversity.tsv"), row.names = 1, check.names = FALSE)
  common <- intersect(rownames(meta), rownames(tab)); check_n(common, m$metric)
  d <- meta[common, ]; d$value <- tab[common, 1]
  alpha_long[[i]] <- row_of(sampleid = rownames(d), metric = m$metric, db = m$db,
                            Diet = as.character(d$Diet), Condition = as.character(d$Condition), value = d$value)
  av <- type2_lm(d$value, d)
  # Weak at n = 19; flags, not gates.
  fit <- lm(value ~ Diet * Condition, data = d)
  shapiro_p  <- tryCatch(shapiro.test(residuals(fit))$p.value, error = function(e) NA_real_)
  levene_p   <- tryCatch({ z <- abs(residuals(fit)); anova(lm(z ~ interaction(d$Diet, d$Condition)))[1, "Pr(>F)"] },
                         error = function(e) NA_real_)
  for (term in TERMS)
    results[[length(results) + 1L]] <- row_of(analysis = "alpha", metric = m$metric, db = m$db, term = term,
      n = nrow(d), statistic = av[term, "F value"], effect = NA_real_, p_raw = av[term, "Pr(>F)"],
      resid_normal_p = shapiro_p, var_homog_p = levene_p)
  message(sprintf("[alpha] %-9s n=%d  Diet p=%.4f  Condition p=%.4f  Interaction p=%.4f",
                  m$metric, nrow(d), av["Diet", "Pr(>F)"], av["Condition", "Pr(>F)"], av["Diet:Condition", "Pr(>F)"]))
}

# ---------- beta ----------
for (i in seq_len(nrow(beta_metrics))) {
  m <- beta_metrics[i, ]
  dm <- as.matrix(read.delim(export_qza(file.path(CORE, m$file), "distance-matrix.tsv"), row.names = 1, check.names = FALSE))
  common <- intersect(rownames(meta), rownames(dm)); check_n(common, m$metric)
  d <- meta[common, ]
  dist_sub <- as.dist(dm[common, common])

  set.seed(42); ad1 <- adonis2(dist_sub ~ Diet * Condition, data = d, permutations = NPERM, by = "terms")
  set.seed(42); ad2 <- adonis2(dist_sub ~ Condition * Diet, data = d, permutations = NPERM, by = "terms")
  ad <- rbind(Diet = ad2["Diet", ], Condition = ad1["Condition", ], "Diet:Condition" = ad1["Diet:Condition", ])
  for (term in TERMS)
    results[[length(results) + 1L]] <- row_of(analysis = "beta_permanova", metric = m$metric, db = m$db,
      term = term, n = length(common), statistic = ad[term, "F"], effect = ad[term, "R2"],
      p_raw = ad[term, "Pr(>F)"], resid_normal_p = NA_real_, var_homog_p = NA_real_)
  # Type II R2 do not sum with the residual to 1 in an unbalanced design; the
  # residual share is reported separately for the R2 figure.
  results[[length(results) + 1L]] <- row_of(analysis = "beta_permanova", metric = m$metric, db = m$db,
    term = "Residual", n = length(common), statistic = NA_real_, effect = ad1["Residual", "R2"],
    p_raw = NA_real_, resid_normal_p = NA_real_, var_homog_p = NA_real_)

  # PERMANOVA cannot separate a centroid shift from a spread difference.
  # "Cell" is the 6-way Diet x Condition grouping PERMANOVA's location test
  # assumes equally spread.
  groupings <- list(Diet = d$Diet, Condition = d$Condition,
                    Cell = droplevels(interaction(d$Diet, d$Condition, sep = "/")))
  for (fac in names(groupings)) {
    bd <- betadisper(dist_sub, groupings[[fac]])
    set.seed(42); pt <- permutest(bd, permutations = NPERM)
    results[[length(results) + 1L]] <- row_of(analysis = "beta_dispersion", metric = m$metric, db = m$db,
      term = fac, n = length(common), statistic = pt$tab[1, "F"], effect = NA_real_,
      p_raw = pt$tab[1, "Pr(>F)"], resid_normal_p = NA_real_, var_homog_p = NA_real_)
    disper_long[[length(disper_long) + 1L]] <- row_of(sampleid = names(bd$distances), metric = m$metric,
      db = m$db, grouping = fac, group = as.character(groupings[[fac]]), distance = unname(bd$distances))
    if (fac == "Cell") z_cell <- unname(bd$distances)
  }
  # Split the Cell dispersion by term: PERMDISP is an ANOVA on
  # distance-to-centroid, so fit the factorial (Type II) on the Cell distances
  # and test each term by permuting distances across samples.
  zF <- function(z) type2_lm(z, d)[TERMS, "F value"]
  obs <- zF(z_cell)
  set.seed(42); perm <- replicate(NPERM, zF(sample(z_cell)))
  zp <- (rowSums(perm >= obs) + 1) / (NPERM + 1)
  for (j in 1:3)
    results[[length(results) + 1L]] <- row_of(analysis = "beta_dispersion_factorial", metric = m$metric,
      db = m$db, term = TERMS[j], n = length(common), statistic = obs[j], effect = NA_real_,
      p_raw = zp[j], resid_normal_p = NA_real_, var_homog_p = NA_real_)
  message(sprintf("[beta]  %-18s Diet p=%.4f (R2=%.3f)  Condition p=%.4f (R2=%.3f)  Interaction p=%.4f",
                  m$metric, ad["Diet", "Pr(>F)"], ad["Diet", "R2"], ad["Condition", "Pr(>F)"],
                  ad["Condition", "R2"], ad["Diet:Condition", "Pr(>F)"]))
}

res <- do.call(rbind, results)
# BH within each analysis x term family: the metrics in one family are
# alternative views of the same samples, not independent questions.
res$p_BH <- NA_real_
fam <- paste(res$analysis, res$term)
for (f in unique(fam)) { k <- fam == f & !is.na(res$p_raw); res$p_BH[k] <- p.adjust(res$p_raw[k], method = "BH") }

outfile <- file.path(TAB_DIR, "diversity_factorial.tsv")
write.table(res, outfile, sep = "\t", quote = FALSE, row.names = FALSE)
message("\n[diversity] Wrote ", outfile)
write.table(do.call(rbind, alpha_long), file.path(TAB_DIR, "diversity_alpha_values.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(do.call(rbind, disper_long), file.path(TAB_DIR, "diversity_dispersion_distances.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("\n[diversity] Significant after BH (q < 0.05):")
sig <- res[!is.na(res$p_BH) & res$p_BH < 0.05, ]
if (nrow(sig) == 0) message("    none") else
  for (j in seq_len(nrow(sig))) message(sprintf("    %-26s %-18s %-14s p=%.4g  q=%.4g",
    sig$analysis[j], sig$metric[j], sig$term[j], sig$p_raw[j], sig$p_BH[j]))
