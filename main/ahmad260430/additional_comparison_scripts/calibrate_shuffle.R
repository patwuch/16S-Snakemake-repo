#!/usr/bin/env Rscript
# How many abundance calls MaAsLin3 and limma make when there is nothing to find
# (rule fixed in ANALYSIS_LOG.md, 2026-10-23). The six group labels are shuffled
# across the 24 mice, 200 times (seed 1); each mouse keeps its counts and read
# depth. Each shuffle runs the real pipeline's abundance step: same features,
# all-groups rule (on the shuffled groups), rows, families and BH. Robustness
# (q < 0.05, same sign, in all 24 leave-one-mouse-out refits) is computed only
# for shuffles with a call, since a shuffle without one cannot have a robust one.
#
# Before shuffling, the unshuffled labels must reproduce the real runs' call
# counts.
#
# Needs containers/maaslin3-env.sif with limma:
#   Rscript calibrate_shuffle.R [--shuffles 200]

THIS_DIR <- dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))[1]))
args <- commandArgs(trailingOnly = TRUE)
source(file.path(THIS_DIR, "maaslin3_common.R"))
source(file.path(THIS_DIR, "limma_common.R"))
i <- match("--shuffles", args)
B <- if (is.na(i)) 200L else as.integer(args[i + 1L])

# run_maaslin3_rows.R's abundance table: NA where the taxon is missing from a group
# or the fit errored, BH within families.
maaslin3_rows <- function(rows) {
  fits <- lapply(setNames(FITS, FITS), function(fn) fit(fn, rows, abundance_only = TRUE)$abundance)
  out <- do.call(rbind, lapply(seq_along(ROWS), function(k) {
    x <- fits[[ROWS[[k]]$fit]]; x <- x[x$name == ROWS[[k]]$term, ]
    p <- ifelse(is.na(x$error), x$pval_individual, NA)
    p[!identifiable(x$feature, k, rows)] <- NA
    data.frame(ROWTAB[k, c("row", "stratum", "family")], feature = x$feature, coef = x$coef, p = p,
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  out$q <- ave(out$p, out$family, FUN = function(p) p.adjust(p, "BH"))
  out
}
METHODS <- list(maaslin3 = maaslin3_rows, limma = limma_rows)
key <- function(d) paste(d$row, d$stratum, d$feature)

score <- function(method) {
  full <- METHODS[[method]](seq_len(n))
  calls <- full[!is.na(full$q) & full$q < ALPHA, ]
  hold <- rep(TRUE, nrow(calls))
  for (s in seq_len(n)) {
    if (!any(hold)) break
    r <- METHODS[[method]](setdiff(seq_len(n), s))
    r <- r[match(key(calls), key(r)), ]
    hold <- hold & !is.na(r$q) & r$q < ALPHA & sign(r$coef) == sign(calls$coef)
  }
  c(calls = nrow(calls), robust = sum(hold), taxa_robust = length(unique(calls$feature[hold])))
}

labels <- meta[, c("Group", "Treatment", "Radiation")]
real <- sapply(names(METHODS), score)
message("[calibration] real labels: ", paste(names(METHODS), real["calls", ], real["robust", ], sep = " calls/robust ", collapse = "; "))
stopifnot(real["calls", "maaslin3"] == 26, real["robust", "maaslin3"] == 14, real["calls", "limma"] == 38)

set.seed(1)
perms <- replicate(B, sample(n), simplify = FALSE)
t0 <- Sys.time()
out <- do.call(rbind, lapply(seq_len(B), function(b) {
  meta[, c("Group", "Treatment", "Radiation")] <<- labels[perms[[b]], ]
  sc <- sapply(names(METHODS), score)
  if (b %% 10 == 0) message(sprintf("[calibration] %d/%d (%.1f min)", b, B, as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  data.frame(shuffle = b, method = names(METHODS), calls = sc["calls", ], robust = sc["robust", ],
             taxa_robust = sc["taxa_robust", ], row.names = NULL)
}))
meta[, c("Group", "Treatment", "Radiation")] <- labels
write.table(out, file.path(MA_DIR, paste0("calibration_shuffle_", TAG, ".tsv")), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(split(out, out$method), function(d) data.frame(
  method = d$method[1], shuffles = nrow(d), mean_calls = mean(d$calls), any_call = mean(d$calls > 0),
  mean_robust = mean(d$robust), any_robust = mean(d$robust > 0), max_robust = max(d$robust))))
summ$real_calls <- real["calls", summ$method]; summ$real_robust <- real["robust", summ$method]
summ$real_robust_beaten <- sapply(summ$method, function(m) mean(out$robust[out$method == m] >= real["robust", m]))
write.table(summ, file.path(MA_DIR, paste0("calibration_shuffle_", TAG, "_summary.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, digits = 3, right = FALSE)
