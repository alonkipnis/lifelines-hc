## Numerical parity: R port (hchg) vs the reference Python implementation
## (lifelines-hc). Reference values were produced by the Python package on the
## same input files and exported as CSV.
##
## Run after installing hchg:  Rscript parity_check.R

library(hchg)
dir <- "/mnt/user-data/uploads/lifelines-hc/.parity_tmp"

ref <- read.csv(file.path(dir, "python_stats.csv"), stringsAsFactors = FALSE)
funcs <- list(hc = higher_criticism_test, bj = berk_jones_test,
              fisher = fisher_combination_test, minp = min_p_test)

dat <- list()
for (cs in unique(ref$case)) {
  dat[[cs]] <- list(A = read.csv(file.path(dir, paste0(cs, "_A.csv"))),
                    B = read.csv(file.path(dir, paste0(cs, "_B.csv"))))
}

ref$R <- NA_real_
for (i in seq_len(nrow(ref))) {
  d <- dat[[ref$case[i]]]
  K <- if (ref$K[i] < 0) NULL else ref$K[i]
  ref$R[i] <- unname(funcs[[ref$stat[i]]](d$A$t, d$B$t, d$A$e, d$B$e,
                     alternative = ref$alt[i], n_intervals_to_pool = K)$statistic)
}

## Both implementations return Inf when the exact Berk-Jones p-value underflows
## on unpooled continuous data; agreeing infinities count as agreement.
both_inf <- is.infinite(ref$R) & is.infinite(ref$python) &
            sign(ref$R) == sign(ref$python)

ref$rel <- ifelse(abs(ref$python) > 1e-12,
                  abs(ref$R - ref$python) / abs(ref$python),
                  abs(ref$R - ref$python))
ref$rel[both_inf] <- 0

cat("=== test statistics: R vs Python ===\n")
print(ref[ref$K == 50 & ref$alt == "both", c("case", "stat", "R", "python", "rel")],
      row.names = FALSE, digits = 10)

cat("\nComparisons:", nrow(ref), "\n")
cat("Agreeing infinities (Berk-Jones underflow, unpooled):", sum(both_inf), "\n")
cat("Worst relative difference over the", sum(!both_inf), "finite comparisons:",
    format(max(ref$rel[!both_inf]), scientific = TRUE), "\n\n")

pref <- read.csv(file.path(dir, "python_pvals.csv"), stringsAsFactors = FALSE)
worst_pv <- 0
for (cs in unique(pref$case)) {
  d <- dat[[cs]]
  pv <- event_pvalues(d$A$t, d$B$t, d$A$e, d$B$e,
                      alternative = "greater", n_intervals_to_pool = 50)
  py <- pref$python[pref$case == cs]
  stopifnot(length(pv) == length(py))
  worst_pv <- max(worst_pv, max(abs(pv - py)))
  cat(sprintf("per-interval p-values [%-8s] n=%d  max abs diff = %s\n",
              cs, length(pv), format(max(abs(pv - py)), scientific = TRUE)))
}
cat("\nWorst p-value difference:", format(worst_pv, scientific = TRUE), "\n")

stopifnot(max(ref$rel) < 1e-10, worst_pv < 1e-12)
cat("\nPARITY CONFIRMED.\n")
