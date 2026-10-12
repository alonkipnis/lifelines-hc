## Numerical parity: R port (survHC) vs the reference Python implementation
## (lifelines-hc). Reference values were produced by the Python package on the
## same input files and exported as CSV.
##
## Run after installing survHC:  Rscript parity_check.R

library(survHC)
dir <- "/mnt/user-data/uploads/lifelines-hc/.parity_tmp"

ref <- read.csv(file.path(dir, "python_stats.csv"), stringsAsFactors = FALSE)
funcs <- list(hc = higher_criticism_test, bj = berk_jones_test,
              fisher = fisher_combination_test, minp = min_p_test)

dat <- list()
for (cs in unique(ref$case)) {
  dat[[cs]] <- list(A = read.csv(file.path(dir, paste0(cs, "_A.csv"))),
                    B = read.csv(file.path(dir, paste0(cs, "_B.csv"))))
}

## The same two samples arranged as one data frame, for the formula interface.
## Group A must come first so that it is the first factor level.
frames <- lapply(dat, function(d) data.frame(
  t = c(d$A$t, d$B$t),
  e = c(d$A$e, d$B$e),
  arm = factor(rep(c("A", "B"), c(nrow(d$A), nrow(d$B))), levels = c("A", "B"))
))

ref$R <- NA_real_
ref$R_formula <- NA_real_
for (i in seq_len(nrow(ref))) {
  d <- dat[[ref$case[i]]]
  K <- if (ref$K[i] < 0) NULL else ref$K[i]
  ## lifelines-hc requests hc_dj2008, which is also this package's default.
  ref$R[i] <- unname(funcs[[ref$stat[i]]](d$A$t, d$B$t, d$A$e, d$B$e,
                     alternative = ref$alt[i], n_intervals_to_pool = K)$statistic)
  ## The formula interface must reach the same number on the same data.
  ref$R_formula[i] <- unname(funcs[[ref$stat[i]]](
                        cbind(t, e) ~ arm, data = frames[[ref$case[i]]],
                        alternative = ref$alt[i],
                        n_intervals_to_pool = K)$statistic)
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
  pv_f <- event_pvalues(cbind(t, e) ~ arm, data = frames[[cs]],
                        alternative = "greater", n_intervals_to_pool = 50)
  stopifnot(identical(pv, pv_f))
  py <- pref$python[pref$case == cs]
  stopifnot(length(pv) == length(py))
  worst_pv <- max(worst_pv, max(abs(pv - py)))
  cat(sprintf("per-interval p-values [%-8s] n=%d  max abs diff = %s\n",
              cs, length(pv), format(max(abs(pv - py)), scientific = TRUE)))
}
cat("\nWorst p-value difference:", format(worst_pv, scientific = TRUE), "\n")

## The formula interface is a pure rearrangement of its inputs, so it must
## agree with the default method bit for bit, not merely to tolerance.
iface <- identical(ref$R, ref$R_formula)
cat("Formula interface identical to the default method on all",
    nrow(ref), "comparisons:", iface, "\n")

stopifnot(max(ref$rel) < 1e-10, worst_pv < 1e-12, iface)
cat("\nPARITY CONFIRMED (default normalization, donoho-jin2008),",
    "through both interfaces.\n")

## Size of the gap to the other standardisation, for reference
cat("\n=== hc_version variants ===\n")
for (cs in unique(ref$case)) {
  d <- dat[[cs]]
  v <- sapply(c("dj2008", "dj2004", "beta", "star"), function(hv)
    unname(higher_criticism_test(d$A$t, d$B$t, d$A$e, d$B$e,
           n_intervals_to_pool = 50, hc_version = hv)$statistic))
  cat(sprintf("  %-9s dj2008=%.6f  dj2004=%.6f  beta=%.6f  star=%.6f\n",
              cs, v[["dj2008"]], v[["dj2004"]], v[["beta"]], v[["star"]]))
}
