#!/usr/bin/env Rscript
# LinDA (zhouhj1994/LinDA 0.2.0) on the Diet x Condition factorial, fitted to
# exactly the taxa seen in every cell: run_ancombc2_factorial.R --export-only
# writes them to tables/tested_counts_<db><suffix>.tsv.
#
# The three diets are peers. Every report row is read as a direct coefficient
# of a fit whose reference levels make it one. LinDA subtracts a separate mode
# from each coefficient, so the same row read as a contrast of other
# coefficients (linda.wald.test) is shifted by a constant; that version is
# written beside it for comparison.
#
# Stage 1 omnibus ("all three diets equal", within each condition, and "diet
# matters in either condition") is computed under every reference diet, since
# the mode correction makes it depend on which diet is the reference.
#
# The pseudo-count check refits at pseudo.cnt 0.1, 0.5 and 1 with the adaptive
# switch off; a call passes only if q < ALPHA in the primary fit and all three.
#
# Adapted from ahmad260430/additional_comparison_scripts/run_linda_factorial.R.
# Needs containers/linda-env.sif:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/linda-env.sif \
#     Rscript run_linda_factorial.R [--db Greengenes2] [--min-presence 6]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "linda_common.R"))
fits <- lapply(settings, function(z) {
  f <- lapply(seq_len(nrow(refs)), function(i) fit(refs$d_ref[i], refs$c_ref[i], z))
  names(f) <- paste(refs$d_ref, refs$c_ref); f })
names(fits) <- c("primary", paste0("pc", PSEUDO))
primary_zero <- fits$primary[["Cellulose Control"]]$zero_handling
message("[linda] primary fit zero handling: ", primary_zero)

rows_for <- function(f) {
  out <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
    o <- f[[paste(ROWS$d_ref[k], ROWS$c_ref[k])]]$output[[ROWS$coef[k]]]
    data.frame(ROWS[k, c("row", "stratum", "family")], taxon = rownames(o),
               log2FC = o$log2FoldChange, lfcSE = o$lfcSE, df = o$df, p = o$pvalue,
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  out$q <- ave(out$p, out$family, FUN = function(p) p.adjust(p, "BH"))
  out
}
res <- rows_for(fits$primary)
for (s in names(fits)[-1]) {
  r <- rows_for(fits[[s]])
  stopifnot(identical(r$taxon, res$taxon), identical(r$row, res$row))
  res[[paste0("q_", s)]] <- r$q
}
res$passed_pseudo <- Reduce(`&`, lapply(c("q", paste0("q_pc", PSEUDO)), function(c) res[[c]] < ALPHA))

base <- fits$primary[["Cellulose Control"]]
beta_c <- base$wald$beta - base$wald$bias
stopifnot(identical(rownames(beta_c), c("(Intercept)", "DietInulin", "DietRS", "ConditionIR",
                                        "DietInulin:ConditionIR", "DietRS:ConditionIR")))
wald <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  w <- linda.wald.test(base, L[k, , drop = FALSE], model = "LM", alpha = ALPHA)
  data.frame(row = ROWS$row[k], stratum = ROWS$stratum[k], taxon = rownames(w),
             wald_log2FC = as.vector(L[k, , drop = FALSE] %*% beta_c), wald_p = w$pvalue)
}))
res <- merge(res, wald, by = c("row", "stratum", "taxon"), sort = FALSE)
res$wald_q <- ave(res$wald_p, res$family, FUN = function(p) p.adjust(p, "BH"))
res$wald_shift <- res$wald_log2FC - res$log2FC
res$zero_handling <- primary_zero
res <- res[order(match(paste(res$row, res$stratum), paste(ROWS$row, ROWS$stratum)), res$p), ]
write.table(res, file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, ".tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Stage 1 omnibus under every reference diet. "diet_<cond>": the two diet
# coefficients of a fit with <cond> as reference (all three diets equal within
# <cond>). "diet_any": those plus both interaction coefficients (4 df).
omni <- do.call(rbind, lapply(seq_len(nrow(refs)), function(i) {
  f <- fits$primary[[paste(refs$d_ref[i], refs$c_ref[i])]]
  cf <- rownames(f$wald$beta)
  tests <- setNames(list(grepl("^Diet", cf) & !grepl(":", cf)), paste0("diet_", refs$c_ref[i]))
  # diet_any is the same hypothesis under either condition reference; take it
  # once per diet reference, from the Control-reference fit.
  if (refs$c_ref[i] == "Control") tests$diet_any <- grepl("^Diet", cf)
  do.call(rbind, lapply(names(tests), function(t) {
    w <- linda.wald.test(f, diag(length(cf))[tests[[t]], , drop = FALSE], model = "LM", alpha = ALPHA)
    data.frame(test = t, ref = refs$d_ref[i], taxon = rownames(w), F = w$Fstat, df1 = w$df1,
               p = w$pvalue, q = p.adjust(w$pvalue, "BH"))
  }))
}))
write.table(omni, file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, "_omnibus.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

omni_summ <- do.call(rbind, lapply(split(omni, omni$test), function(d) {
  hits <- lapply(split(d, d$ref), function(x) sort(x$taxon[x$q < ALPHA]))
  maxp <- tapply(d$p, d$taxon, max)
  data.frame(test = d$test[1],
             calls = paste(sprintf("%s %d", names(hits), lengths(hits)), collapse = ", "),
             identical_sets = length(unique(hits)) == 1,
             in_all = length(Reduce(intersect, hits)), in_any = length(Reduce(union, hits)),
             maxp_bh = sum(p.adjust(maxp, "BH") < ALPHA))
}))
rownames(omni_summ) <- NULL

# Mode diagnostic: the raw (uncorrected) coefficients of each direct fit, the
# mode LinDA subtracts, and two other location estimates. A well-separated
# unimodal peak makes all three agree; a flat or bimodal density does not.
n <- ncol(counts)
diag_tab <- do.call(rbind, lapply(names(fits), function(s) do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  f <- fits[[s]][[paste(ROWS$d_ref[k], ROWS$c_ref[k])]]
  raw <- f$wald$beta[ROWS$coef[k], ]
  data.frame(setting = s, row = ROWS$row[k], stratum = ROWS$stratum[k],
             linda_mode = f$wald$bias[match(ROWS$coef[k], rownames(f$wald$beta))],
             median = median(raw),
             hsm = mlv(sqrt(n) * raw, method = "hsm") / sqrt(n),
             sd_raw = sd(raw))
}))))
write.table(diag_tab, file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, "_mode_diagnostic.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

dens <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  f <- fits$primary[[paste(ROWS$d_ref[k], ROWS$c_ref[k])]]
  data.frame(panel = trimws(paste(ROWS$row[k], ifelse(ROWS$stratum[k] == "", "", paste0("(", ROWS$stratum[k], ")")))),
             raw = f$wald$beta[ROWS$coef[k], ])
}))
dens$panel <- factor(dens$panel, levels = unique(dens$panel))
marks <- data.frame(panel = factor(levels(dens$panel), levels = levels(dens$panel)),
                    mode = diag_tab$linda_mode[diag_tab$setting == "primary"],
                    median = diag_tab$median[diag_tab$setting == "primary"])
p <- ggplot(dens, aes(raw)) +
  geom_density(fill = "grey85", colour = "grey30") + geom_rug(alpha = 0.6) +
  geom_vline(data = marks, aes(xintercept = mode), colour = "#C0392B") +
  geom_vline(data = marks, aes(xintercept = median), colour = "#2471A3", linetype = "dashed") +
  facet_wrap(~ panel, scales = "free", ncol = 3) +
  labs(x = "uncorrected coefficient (log2)", y = "density",
       title = sprintf("LinDA coefficients before bias correction, %d taxa", nrow(counts)),
       subtitle = sprintf("red: mode LinDA subtracts; blue dashed: median; zero handling: %s", primary_zero)) +
  theme_bw(base_size = 9)
ggsave(file.path(FIG_DIR, paste0("LinDA_factorial_", TAG, "_mode_diagnostic.png")), p,
       width = 10, height = 9, dpi = 150)

summ <- do.call(rbind, lapply(split(res, factor(paste(res$row, res$stratum), levels = unique(paste(res$row, res$stratum)))), function(d)
  data.frame(row = d$row[1], stratum = d$stratum[1], q_primary = sum(d$q < ALPHA),
             passed_pseudo = sum(d$passed_pseudo), q_pc0.1 = sum(d$q_pc0.1 < ALPHA),
             q_pc0.5 = sum(d$q_pc0.5 < ALPHA), q_pc1 = sum(d$q_pc1 < ALPHA),
             wald_q = sum(d$wald_q < ALPHA), wald_shift = round(median(d$wald_shift), 3),
             wald_shift_range = round(diff(range(d$wald_shift)), 6))))
rownames(summ) <- NULL
print(summ, right = FALSE)
print(omni_summ, right = FALSE)
