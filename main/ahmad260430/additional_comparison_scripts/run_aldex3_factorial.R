#!/usr/bin/env Rscript
# ALDEx3 on the Treatment x Radiation factorial for ahmad260430: the primary
# abundance method (settings fixed in ANALYSIS_LOG.md before the first fit).
# Fitted to exactly the taxa ANCOM-BC2 and LinDA test:
# run_ancombc2_factorial.R --export-only writes tables/tested_counts_<db><suffix>.tsv.
#
# Model: aldex(Y, ~ Treatment * Radiation, method = "lm", test = "t.HC3",
# scale = clr.sm, gamma). clr.sm adds N(0, gamma^2) scale noise to every model
# coefficient in each Monte Carlo draw, so each report row is read as a direct
# coefficient of a fit whose reference levels make it one: every row then gets
# the same scale prior. Estimates are log2.
#
# Correction follows ALDEx3's own rule (aldex.pvals: BH within each draw and
# side, mean over draws, smaller side), applied to the report's families rather
# than one coefficient at a time. Per-draw p-values are rebuilt from the
# returned per-draw estimate and std.error with fflm's t(N - P); the script
# checks that this reproduces aldex()'s own p.val and p.val.adj.
#
# Usage:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/ancombc2-env.sif \
#     Rscript run_aldex3_factorial.R [--min-presence 6] [--gamma 0.5] [--seed 42] [--nsample 1000]

suppressPackageStartupMessages(library(ALDEx3))
stopifnot(packageVersion("ALDEx3") == "1.3.1")

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
TAB_DIR   <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial", "tables")

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
DB      <- "Greengenes2"
MINPRES <- as.integer(arg_val("--min-presence", "6"))
GAMMA   <- as.numeric(arg_val("--gamma", "0.5"))
SEED    <- as.integer(arg_val("--seed", "42"))
NSAMPLE <- as.integer(arg_val("--nsample", "1000"))
ALPHA   <- 0.05
TAG     <- sprintf("%s_Day17_minpres%d", DB, MINPRES)
suffix  <- paste0("_gamma", format(GAMMA))
if (SEED != 42L) suffix <- paste0(suffix, "_seed", SEED)
if (NSAMPLE != 1000L) suffix <- paste0(suffix, "_ns", NSAMPLE)

counts <- as.matrix(read.delim(file.path(TAB_DIR, paste0("tested_counts_", TAG, ".tsv")),
                               row.names = 1, check.names = FALSE))
meta <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)[colnames(counts), c("Group", "Treatment", "Radiation")]
stopifnot(!anyNA(meta$Treatment), all(table(meta$Treatment, meta$Radiation) == 4))
message(sprintf("[aldex3] %s: %d taxa x %d samples; gamma %s, seed %d, nsample %d",
                TAG, nrow(counts), ncol(counts), format(GAMMA), SEED, NSAMPLE))

# Report rows and the fit (treatment ref, radiation ref) + coefficient that
# estimates each directly; the same table as run_linda_factorial.R.
ROWS <- data.frame(
  row    = c("Inosine vs Control", "IFA vs Control", "IFA vs Inosine",
             "Inosine vs Control", "IFA vs Control", "IFA vs Inosine",
             "Radiation within Control", "Radiation within Inosine", "Radiation within IFA",
             "Interaction: Inosine vs Control", "Interaction: IFA vs Control",
             "Interaction: IFA vs Inosine"),
  stratum = c(rep("unirradiated", 3), rep("irradiated", 3), rep("", 6)),
  family = c(rep("treatment_unirradiated", 3), rep("treatment_irradiated", 3),
             "radiation_Control", "radiation_Inosine", "radiation_IFA", rep("interaction", 3)),
  t_ref = c("Control", "Control", "Inosine", "Control", "Control", "Inosine",
            "Control", "Inosine", "IFA", "Control", "Control", "Inosine"),
  r_ref = c("None", "None", "None", "IR", "IR", "IR", rep("None", 6)),
  coef  = c("TreatmentInosine", "TreatmentIFA", "TreatmentIFA",
            "TreatmentInosine", "TreatmentIFA", "TreatmentIFA",
            rep("RadiationIR", 3),
            "TreatmentInosine:RadiationIR", "TreatmentIFA:RadiationIR", "TreatmentIFA:RadiationIR"),
  stringsAsFactors = FALSE)

fit <- function(t_ref, r_ref) {
  m <- meta
  m$Treatment <- relevel(factor(m$Treatment, levels = c("Control", "Inosine", "IFA")), t_ref)
  m$Radiation <- relevel(factor(m$Radiation, levels = c("None", "IR")), r_ref)
  set.seed(SEED)
  aldex(counts, ~ Treatment * Radiation, data = m, method = "lm", nsample = NSAMPLE,
        scale = clr.sm, gamma = GAMMA, test = "t.HC3", p.adjust.method = "BH", n.cores = 1,
        return.pars = c("X", "estimate", "std.error", "p.val", "p.val.adj"))
}
refs <- unique(ROWS[, c("t_ref", "r_ref")])
fits <- lapply(seq_len(nrow(refs)), function(i) fit(refs$t_ref[i], refs$r_ref[i]))
names(fits) <- paste(refs$t_ref, refs$r_ref)

# Per-draw two-sided p for each side, as aldex.pvals builds them.
side_p <- function(f, coef) {
  dof <- ncol(counts) - dim(f$estimate)[1]                   # f$X is coefficients x samples
  t <- f$estimate[coef, , ] / f$std.error[coef, , ]          # taxon x draw
  up <- pt(t, dof, lower.tail = FALSE)
  list(lower = pmin(2 * (1 - up), 1), upper = pmin(2 * up, 1))
}
# ALDEx3's combination over a set of (row, taxon) tests: BH within each draw
# and side, mean over draws, smaller side.
combine <- function(sides) {
  adj <- lapply(sides, function(s) apply(s, 2, p.adjust, method = "BH"))
  list(p = pmin(rowMeans(sides$lower), rowMeans(sides$upper)),
       q = pmin(rowMeans(adj$lower), rowMeans(adj$upper)))
}

# The rebuild must reproduce aldex()'s own numbers coefficient by coefficient.
for (k in seq_len(nrow(ROWS))) {
  f <- fits[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]
  own <- combine(side_p(f, ROWS$coef[k]))
  stopifnot(isTRUE(all.equal(unname(own$p), unname(f$p.val[ROWS$coef[k], ]), tolerance = 1e-10)),
            isTRUE(all.equal(unname(own$q), unname(f$p.val.adj[ROWS$coef[k], ]), tolerance = 1e-10)))
}
message("[aldex3] per-draw rebuild reproduces aldex()'s p.val and p.val.adj for all 12 rows")

res <- do.call(rbind, lapply(split(seq_len(nrow(ROWS)), factor(ROWS$family, levels = unique(ROWS$family))), function(ks) {
  sides <- lapply(ks, function(k) side_p(fits[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]], ROWS$coef[k]))
  stacked <- list(lower = do.call(rbind, lapply(sides, `[[`, "lower")),
                  upper = do.call(rbind, lapply(sides, `[[`, "upper")))
  fam <- combine(stacked)
  do.call(rbind, lapply(seq_along(ks), function(j) {
    k <- ks[j]; f <- fits[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]
    idx <- (j - 1) * nrow(counts) + seq_len(nrow(counts))
    data.frame(ROWS[k, c("row", "stratum", "family")], taxon = rownames(counts),
               log2FC = rowMeans(f$estimate[ROWS$coef[k], , ]),
               se = rowMeans(f$std.error[ROWS$coef[k], , ]),
               p = fam$p[idx], q = fam$q[idx],
               q_single_coef = unname(f$p.val.adj[ROWS$coef[k], ]),
               n_zero_mice = rowSums(counts == 0), row.names = NULL, stringsAsFactors = FALSE)
  }))
}))
res <- res[order(match(paste(res$row, res$stratum), paste(ROWS$row, ROWS$stratum)), res$p), ]
res$gamma <- GAMMA

# 95% interval: quantiles of the mixture over draws of estimate + se * t(N - P).
# Its tail mass below zero is the mean over draws of the per-draw one-sided p,
# so the interval excludes zero exactly when the mean two-sided p (before
# aldex.pvals caps each draw at 1) is below 0.05.
dof <- ncol(counts) - 6L
mix_q <- function(est, se, prob) {
  f <- function(x) mean(pt((x - est) / se, dof)) - prob
  uniroot(f, c(min(est - 50 * se), max(est + 50 * se)), tol = 1e-8)$root
}
ci <- t(mapply(function(row, stratum, taxon) {
  k <- which(ROWS$row == row & ROWS$stratum == stratum)
  f <- fits[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]
  est <- f$estimate[ROWS$coef[k], taxon, ]; se <- f$std.error[ROWS$coef[k], taxon, ]
  c(mix_q(est, se, 0.025), mix_q(est, se, 0.975))
}, res$row, res$stratum, res$taxon))
res$ci_lo <- ci[, 1]; res$ci_hi <- ci[, 2]
n_disagree <- sum((res$ci_lo > 0 | res$ci_hi < 0) != (res$p < 0.05))
message(sprintf("[aldex3] 95%% interval excludes zero vs p < 0.05: %d of %d disagree", n_disagree, nrow(res)))

# Model-fitted mean of each Treatment x Radiation group (log2, before any
# centring), mean over draws, from the Control/None fit.
base <- fits[["Control None"]]
stopifnot(identical(dimnames(base$estimate)[[1]], c("(Intercept)", "TreatmentInosine", "TreatmentIFA", "RadiationIR",
                                                   "TreatmentInosine:RadiationIR", "TreatmentIFA:RadiationIR")))
CELLS <- rbind("Control" = c(1, 0, 0, 0, 0, 0), "Inosine" = c(1, 1, 0, 0, 0, 0), "IFA" = c(1, 0, 1, 0, 0, 0),
               "IR" = c(1, 0, 0, 1, 0, 0), "Inosine+IR" = c(1, 1, 0, 1, 1, 0), "IFA+IR" = c(1, 0, 1, 1, 0, 1))
beta_bar <- apply(base$estimate, c(1, 2), mean)                  # coefficient x taxon
cellmeans <- data.frame(taxon = colnames(beta_bar), t(CELLS %*% beta_bar), check.names = FALSE)
write.table(cellmeans, file.path(TAB_DIR, paste0("ALDEx3_factorial_", TAG, suffix, "_groupmeans.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)
out_path <- file.path(TAB_DIR, paste0("ALDEx3_factorial_", TAG, suffix, ".tsv"))
write.table(res, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
message("[aldex3] wrote ", out_path)

summ <- do.call(rbind, lapply(split(res, factor(paste(res$row, res$stratum), levels = unique(paste(res$row, res$stratum)))), function(d)
  data.frame(row = d$row[1], stratum = d$stratum[1], calls = sum(d$q < ALPHA),
             raw_p = sum(d$p < ALPHA), min_q = signif(min(d$q), 2),
             median_log2FC = round(median(d$log2FC), 2))))
rownames(summ) <- NULL
print(summ, right = FALSE)
