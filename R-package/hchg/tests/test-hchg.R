## Base-R test suite (no testthat dependency).
## Any stopifnot() failure aborts R CMD check with a non-zero status.

library(hchg)

ok <- function(label, expr) {
  stopifnot(isTRUE(expr))
  cat("PASS:", label, "\n")
}

set.seed(20261005)

## ---------------------------------------------------------------- structure

r <- higher_criticism_test(rexp(200, 1 / 10), rexp(200, 1 / 10),
                           n_intervals_to_pool = 50)
ok("returns an htest", inherits(r, "htest"))
ok("has a finite statistic", is.finite(unname(r$statistic)))
ok("p-value is NA without permutations", is.na(r$p.value))

r2 <- higher_criticism_test(rexp(150, 1 / 10), rexp(150, 1 / 10),
                            n_intervals_to_pool = 40,
                            n_permutations = 100, seed = 1)
ok("p-value finite with permutations", is.finite(r2$p.value))
ok("p-value in (0, 1]", r2$p.value > 0 && r2$p.value <= 1)
ok("permutation statistics retained", length(r2$permutation_statistics) == 100)

## ------------------------------------------------------------ reproducible

a <- rexp(150, 1 / 10); b <- rexp(150, 1 / 10)
p1 <- higher_criticism_test(a, b, n_intervals_to_pool = 40,
                            n_permutations = 100, seed = 7)$p.value
p2 <- higher_criticism_test(a, b, n_intervals_to_pool = 40,
                            n_permutations = 100, seed = 7)$p.value
ok("seeded p-values reproduce exactly", identical(p1, p2))

## the caller's RNG stream must be left untouched
set.seed(99); before <- runif(1)
set.seed(99); invisible(higher_criticism_test(a, b, n_intervals_to_pool = 40,
                                              n_permutations = 20, seed = 3))
after <- runif(1)
ok("caller's RNG state is restored", isTRUE(all.equal(before, after)))

## ------------------------------------------------------------- directions

g <- higher_criticism_test(a, b, alternative = "greater", n_intervals_to_pool = 40)
l <- higher_criticism_test(a, b, alternative = "less", n_intervals_to_pool = 40)
bo <- higher_criticism_test(a, b, alternative = "both", n_intervals_to_pool = 40)
ok("two-sided is the max of the one-sided statistics",
   isTRUE(all.equal(unname(bo$statistic),
                    max(unname(g$statistic), unname(l$statistic)))))

## ------------------------------------------------------------- censoring

eA <- rbinom(150, 1, 0.8); eB <- rbinom(150, 1, 0.8)
rc <- higher_criticism_test(a, b, eA, eB, n_intervals_to_pool = 40)
ok("handles censoring", is.finite(unname(rc$statistic)))

## censoring everything after a cutoff must equal the t_0 restriction
cut <- stats::median(c(a, b))
eA2 <- ifelse(a > cut, 0, 1); eB2 <- ifelse(b > cut, 0, 1)
a2 <- pmin(a, cut); b2 <- pmin(b, cut)
s_t0 <- higher_criticism_test(a, b, t_0 = cut, n_intervals_to_pool = 30)$statistic
s_manual <- higher_criticism_test(a2, b2, eA2, eB2,
                                  n_intervals_to_pool = 30)$statistic
ok("t_0 restriction matches manual censoring",
   isTRUE(all.equal(unname(s_t0), unname(s_manual))))

## -------------------------------------------------------------- behaviour

## power: a sparse departure should give a larger HC score than the null
null_stat <- higher_criticism_test(rexp(400, 1 / 10), rexp(400, 1 / 10),
                                   n_intervals_to_pool = 60)$statistic
alt_stat <- higher_criticism_test(
  rexp(400, 1 / 10), c(rexp(340, 1 / 10), rexp(60, 1 / 1.2)),
  n_intervals_to_pool = 60)$statistic
ok("HC score is larger under a sparse alternative",
   unname(alt_stat) > unname(null_stat))

## the permutation p-value should not be extreme under the null
pnull <- higher_criticism_test(rexp(300, 1 / 10), rexp(300, 1 / 10),
                               n_intervals_to_pool = 50,
                               n_permutations = 200, seed = 11)$p.value
ok("null p-value is not tiny", pnull > 0.01)

## ------------------------------------------------------- per-interval API

pv <- event_pvalues(a, b, alternative = "greater", n_intervals_to_pool = 40)
ok("event_pvalues returns a numeric vector", is.numeric(pv) && length(pv) > 0)
ok("p-values lie in [0, 1]", all(pv >= 0 & pv <= 1))

pv_both <- event_pvalues(a, b, alternative = "both", n_intervals_to_pool = 40)
ok("alternative='both' returns both directions",
   is.list(pv_both) && all(c("greater", "less") %in% names(pv_both)))

n1 <- length(event_pvalues(a, b, alternative = "greater", n_intervals_to_pool = 20))
n2 <- length(event_pvalues(a, b, alternative = "greater", n_intervals_to_pool = 60))
ok("pooling controls the number of intervals", n1 <= 20 && n2 <= 60 && n2 > n1)

## ------------------------------------------------------------ diagnostics

d <- suspected_deviations(a, b, n_intervals_to_pool = 40)
ok("suspected_deviations returns a data frame", is.data.frame(d))
ok("diagnostic columns present",
   all(c("at_risk_A", "at_risk_B", "observed_A", "observed_B",
         "hypergeom_pvalue", "suspected", "hc_threshold") %in% names(d)))
ok("flagging agrees with the HC threshold",
   all(d$suspected == (d$hypergeom_pvalue <= d$hc_threshold)))

d2 <- suspected_deviations(a, b, alternative = "both", n_intervals_to_pool = 40)
ok("two-sided diagnostics add reverse columns",
   all(c("hypergeom_pvalue_rev", "suspected_rev") %in% names(d2)))

## ------------------------------------------------------------ other tests

for (f in list(berk_jones_test, fisher_combination_test, min_p_test)) {
  res <- f(a, b, n_intervals_to_pool = 40)
  ok("sibling test returns a finite htest statistic",
     inherits(res, "htest") && is.finite(unname(res$statistic)))
}

## ------------------------------------------------------------- validation

ok("rejects mismatched event vector",
   inherits(try(higher_criticism_test(a, b, event_observed_A = c(1, 0)),
                silent = TRUE), "try-error"))
ok("rejects gamma outside (0, 1]",
   inherits(try(higher_criticism_test(a, b, gamma = 0), silent = TRUE),
            "try-error"))
ok("rejects negative durations",
   inherits(try(higher_criticism_test(c(-1, 2, 3), b), silent = TRUE),
            "try-error"))

cat("\nAll tests passed.\n")
