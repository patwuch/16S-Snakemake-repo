#!/usr/bin/env Rscript
# LinDA (zhouhj1994/LinDA 0.2.0) on the Treatment x Radiation factorial, fitted
# to exactly the taxa ANCOM-BC2 tests: run_ancombc2_factorial.R --export-only
# writes them to tables/tested_counts_<db><suffix>.tsv.
#
# Every report row is read as a direct coefficient of a fit whose reference
# levels make it one. LinDA subtracts a separate mode from each coefficient, so
# the same row read as a contrast of other coefficients (linda.wald.test) is
# shifted by a constant; that version is written beside it for comparison.
#
# The pseudo-count check refits at pseudo.cnt 0.1, 0.5 and 1 with the adaptive
# switch off; a call passes only if q < ALPHA in the primary fit and all three.
#
# Needs containers/linda-env.sif:
#   apptainer exec --bind /srv/projects/16S-Snakemake-repo:/srv/projects/16S-Snakemake-repo \
#     /srv/projects/16S-Snakemake-repo/containers/linda-env.sif \
#     Rscript run_linda_factorial.R [--db Greengenes2] [--min-presence 6]

suppressPackageStartupMessages({ library(LinDA); library(modeest); library(ggplot2) })
stopifnot(packageVersion("LinDA") == "0.2.0")

THIS_DIR  <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
STUDY_DIR <- dirname(THIS_DIR)
OUT_DIR   <- file.path(STUDY_DIR, "comparison-Treatment-x-Radiation-factorial")
TAB_DIR   <- file.path(OUT_DIR, "tables")
FIG_DIR   <- file.path(OUT_DIR, "figures")
dir.create(FIG_DIR, showWarnings = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_val <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1L]
}
DB     <- arg_val("--db", "Greengenes2")
MINPRES <- as.integer(arg_val("--min-presence", "6"))
ALPHA  <- 0.05
PSEUDO <- c(0.1, 0.5, 1)
FORMULA <- "~ Treatment * Radiation"
# "--prev-cut F" replaces both upstream filters with LinDA's own prevalence rule
# (present in >= F of samples, as linda()'s prev.cut) on the collapsed table
# (run_ancombc2_factorial.R --export-only), so taxa absent from a whole cell
# can be fitted. Applied here rather than through linda() so every fit and the
# figure see the same taxa.
PREV_CUT <- as.numeric(arg_val("--prev-cut", NA))
if (is.na(PREV_CUT)) {
  TAG <- sprintf("%s_Day17_minpres%d", DB, MINPRES)
  counts <- as.matrix(read.delim(file.path(TAB_DIR, paste0("tested_counts_", TAG, ".tsv")),
                                 row.names = 1, check.names = FALSE))
} else {
  TAG <- sprintf("%s_Day17_prev%s", DB, format(PREV_CUT))
  counts <- as.matrix(read.delim(file.path(TAB_DIR, sprintf("collapsed_counts_%s_Day17.tsv", DB)),
                                 row.names = 1, check.names = FALSE))
  counts <- counts[rowSums(counts > 0) / ncol(counts) >= PREV_CUT, , drop = FALSE]
}
meta <- read.delim(file.path(STUDY_DIR, "metadata.tsv"), row.names = 1, check.names = FALSE,
                   stringsAsFactors = FALSE)[colnames(counts), c("Group", "Treatment", "Radiation")]
stopifnot(!anyNA(meta$Treatment), all(table(meta$Treatment, meta$Radiation) == 4))
message(sprintf("[linda] %s: %d taxa x %d samples", TAG, nrow(counts), ncol(counts)))

# One fit per reference pair; `adaptive = TRUE` is the primary fit, a number is
# a fixed pseudo-count with the adaptive switch off.
fit <- function(t_ref, r_ref, zero) {
  m <- meta
  # Explicit orders: alphabetical would put IFA before Inosine and IR before None.
  m$Treatment <- relevel(factor(m$Treatment, levels = c("Control", "Inosine", "IFA")), t_ref)
  m$Radiation <- relevel(factor(m$Radiation, levels = c("None", "IR")), r_ref)
  msg <- NULL
  o <- if (identical(zero, "adaptive")) {
    msg <- capture.output(r <- linda(counts, m, FORMULA, type = "count", adaptive = TRUE,
                                     prev.cut = 0, lib.cut = 1, winsor.quan = NULL))
    r
  } else {
    capture.output(r <- linda(counts, m, FORMULA, type = "count", adaptive = FALSE,
                              imputation = FALSE, pseudo.cnt = zero,
                              prev.cut = 0, lib.cut = 1, winsor.quan = NULL))
    r
  }
  o$zero_handling <- if (is.null(msg)) paste("pseudo-count", zero) else
    if (any(grepl("Imputation", msg))) "imputation" else "pseudo-count 0.5"
  o
}

# Report rows: the fit (treatment ref, radiation ref) and coefficient that
# estimate each one directly, and the same row as a contrast L of the
# Control/None fit's coefficients (Intercept, Inosine, IFA, IR, Ino:IR, IFA:IR).
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
L <- rbind(c(0, 1, 0, 0, 0, 0), c(0, 0, 1, 0, 0, 0), c(0, -1, 1, 0, 0, 0),
           c(0, 1, 0, 0, 1, 0), c(0, 0, 1, 0, 0, 1), c(0, -1, 1, 0, -1, 1),
           c(0, 0, 0, 1, 0, 0), c(0, 0, 0, 1, 1, 0), c(0, 0, 0, 1, 0, 1),
           c(0, 0, 0, 0, 1, 0), c(0, 0, 0, 0, 0, 1), c(0, 0, 0, 0, -1, 1))

refs <- unique(ROWS[, c("t_ref", "r_ref")])
settings <- c(list("adaptive"), as.list(PSEUDO))
fits <- lapply(settings, function(z) {
  f <- lapply(seq_len(nrow(refs)), function(i) fit(refs$t_ref[i], refs$r_ref[i], z))
  names(f) <- paste(refs$t_ref, refs$r_ref); f })
names(fits) <- c("primary", paste0("pc", PSEUDO))
primary_zero <- fits$primary[["Control None"]]$zero_handling
message("[linda] primary fit zero handling: ", primary_zero)

# Direct estimates for every row under one setting, BH within each family.
rows_for <- function(f) {
  out <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
    o <- f[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]$output[[ROWS$coef[k]]]
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

# The same rows as Wald contrasts on the Control/None primary fit.
base <- fits$primary[["Control None"]]
beta_c <- base$wald$beta - base$wald$bias
stopifnot(identical(rownames(beta_c), c("(Intercept)", "TreatmentInosine", "TreatmentIFA", "RadiationIR",
                                        "TreatmentInosine:RadiationIR", "TreatmentIFA:RadiationIR")))
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

# Mode diagnostic: the raw (uncorrected) coefficients of each direct fit, the
# mode LinDA subtracts, and two other location estimates. A well-separated
# unimodal peak makes all three agree; a flat or bimodal density does not.
n <- ncol(counts)
diag <- do.call(rbind, lapply(names(fits), function(s) do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  f <- fits[[s]][[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]
  raw <- f$wald$beta[ROWS$coef[k], ]
  data.frame(setting = s, row = ROWS$row[k], stratum = ROWS$stratum[k],
             linda_mode = f$wald$bias[match(ROWS$coef[k], rownames(f$wald$beta))],
             median = median(raw),
             hsm = mlv(sqrt(n) * raw, method = "hsm") / sqrt(n),
             sd_raw = sd(raw))
}))))
write.table(diag, file.path(TAB_DIR, paste0("LinDA_factorial_", TAG, "_mode_diagnostic.tsv")),
            sep = "\t", quote = FALSE, row.names = FALSE)

dens <- do.call(rbind, lapply(seq_len(nrow(ROWS)), function(k) {
  f <- fits$primary[[paste(ROWS$t_ref[k], ROWS$r_ref[k])]]
  data.frame(panel = trimws(paste(ROWS$row[k], ifelse(ROWS$stratum[k] == "", "", paste0("(", ROWS$stratum[k], ")")))),
             raw = f$wald$beta[ROWS$coef[k], ])
}))
dens$panel <- factor(dens$panel, levels = unique(dens$panel))
marks <- data.frame(panel = factor(levels(dens$panel), levels = levels(dens$panel)),
                    mode = diag$linda_mode[diag$setting == "primary"],
                    median = diag$median[diag$setting == "primary"])
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

# Summary per row: calls in the primary fit, those surviving the pseudo-count
# check, and the Wald version's calls and constant shift.
summ <- do.call(rbind, lapply(split(res, factor(paste(res$row, res$stratum), levels = unique(paste(res$row, res$stratum)))), function(d)
  data.frame(row = d$row[1], stratum = d$stratum[1], q_primary = sum(d$q < ALPHA),
             passed_pseudo = sum(d$passed_pseudo), q_pc0.1 = sum(d$q_pc0.1 < ALPHA),
             q_pc0.5 = sum(d$q_pc0.5 < ALPHA), q_pc1 = sum(d$q_pc1 < ALPHA),
             wald_q = sum(d$wald_q < ALPHA), wald_shift = round(median(d$wald_shift), 3),
             wald_shift_range = round(diff(range(d$wald_shift)), 6))))
rownames(summ) <- NULL
print(summ, right = FALSE)
