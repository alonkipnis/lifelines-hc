## Base-R test suite (no testthat dependency).
## Any stopifnot() failure aborts R CMD check with a non-zero status.

library(survHC)

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

## --------------------------------------------------------- HC variants

vals <- sapply(c("dj2008","dj2004","beta","star"), function(v)
  unname(higher_criticism_test(a, b, n_intervals_to_pool = 40, hc_version = v)$statistic))
ok("all four hc_version values give finite statistics", all(is.finite(vals)))
ok("dj2008 is the default",
   isTRUE(all.equal(unname(vals[["dj2008"]]),
     unname(higher_criticism_test(a, b, n_intervals_to_pool = 40)$statistic))))
ok("the variants are not all identical", length(unique(round(vals, 10))) > 1)
ok("hc_version is recorded on the result",
   identical(higher_criticism_test(a, b, n_intervals_to_pool = 40)$hc_version, "dj2008"))
ok("rejects an unknown hc_version",
   inherits(try(higher_criticism_test(a, b, hc_version = "nope"), silent = TRUE),
            "try-error"))

## closed-form checks of the two Donoho-Jin definitions
pv <- c(0.001, 0.02, 0.3, 0.4, 0.55, 0.6, 0.7, 0.8, 0.9, 0.95)
N <- length(pv); spv <- sort(pv); EPS <- 1 / (1e4 + N^2)
uu <- seq_len(N) / N; uu[N] <- uu[N] - EPS
imax <- as.integer(0.2 * N + 0.5)
ok("dj2008 matches its closed form",
   isTRUE(all.equal(unname(survHC:::.mt_hc(pv, 0.2, "dj2008")[1]),
                    max(((uu - spv) / sqrt(uu * (1 - uu) / N))[seq_len(imax)]))))
ok("dj2004 matches its closed form (paper eq. 2)",
   isTRUE(all.equal(unname(survHC:::.mt_hc(pv, 0.2, "dj2004")[1]),
                    max(((uu - spv) / sqrt(spv * (1 - spv) / N))[seq_len(imax)]))))

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

## ======================================================= formula interface
##
## The formula method must be a pure rearrangement of its inputs: every
## statistic it produces has to be bit-identical to the one the default method
## gives on the same data, split by hand. Anything less means the two
## documented interfaces disagree.

set.seed(20261012)
nA <- 220L; nB <- 260L
tA <- rexp(nA, 1 / 10)
tB <- c(rexp(nB - 60L, 1 / 10), rexp(60L, 1 / 2))   # excess hazard in B
eA <- rbinom(nA, 1, 0.85)
eB <- rbinom(nB, 1, 0.85)

trial <- data.frame(
  time = c(tA, tB),
  status = c(eA, eB),
  arm = factor(rep(c("ctrl", "trt"), c(nA, nB)), levels = c("ctrl", "trt"))
)

## A structural stand-in for survival::Surv(type = "right"), so the Surv code
## path is exercised even where 'survival' is not installed. The real package
## is used as well, further down, whenever it is available.
as_surv <- function(time, status) {
  structure(cbind(time = as.numeric(time), status = as.numeric(status)),
            type = "right", class = "Surv")
}

## ------------------------------------------- the four tests, exactly equal

for (nm in c("higher_criticism_test", "berk_jones_test",
             "fisher_combination_test", "min_p_test")) {
  f <- get(nm)
  v <- f(tA, tB, eA, eB, n_intervals_to_pool = 60,
         n_permutations = 200, seed = 7)
  m <- f(cbind(time, status) ~ arm, data = trial, n_intervals_to_pool = 60,
         n_permutations = 200, seed = 7)
  s <- f(as_surv(time, status) ~ arm, data = trial, n_intervals_to_pool = 60,
         n_permutations = 200, seed = 7)
  ok(paste(nm, "statistic identical through the formula interface"),
     identical(unname(v$statistic), unname(m$statistic)))
  ok(paste(nm, "permutation p-value identical through the formula interface"),
     identical(v$p.value, m$p.value))
  ok(paste(nm, "Surv() and cbind() left-hand sides agree"),
     identical(unname(m$statistic), unname(s$statistic)) &&
       identical(m$p.value, s$p.value))
  ok(paste(nm, "records the group levels in order (A, B)"),
     identical(m$group_levels, c("ctrl", "trt")))
  ok(paste(nm, "relabels data.name with the formula"),
     grepl("by arm", m$data.name, fixed = TRUE))
}

## every hc_version, every direction
for (v in c("dj2008", "dj2004", "beta", "star")) {
  for (alt in c("both", "greater", "less")) {
    sv <- higher_criticism_test(tA, tB, eA, eB, alternative = alt,
                                hc_version = v, n_intervals_to_pool = 60)
    sf <- higher_criticism_test(cbind(time, status) ~ arm, data = trial,
                                alternative = alt, hc_version = v,
                                n_intervals_to_pool = 60)
    ok(paste0("identical for hc_version=", v, ", alternative=", alt),
       identical(unname(sv$statistic), unname(sf$statistic)))
  }
}

## unpooled, and with the t_0 restriction
ok("identical without pooling",
   identical(unname(higher_criticism_test(tA, tB, eA, eB)$statistic),
             unname(higher_criticism_test(cbind(time, status) ~ arm,
                                          data = trial)$statistic)))
ok("identical under the t_0 restriction",
   identical(unname(higher_criticism_test(tA, tB, eA, eB, t_0 = 12,
                                          n_intervals_to_pool = 40)$statistic),
             unname(higher_criticism_test(cbind(time, status) ~ arm,
                                          data = trial, t_0 = 12,
                                          n_intervals_to_pool = 40)$statistic)))

## ------------------------------------------------ per-interval and tables

ok("event_pvalues identical through the formula interface",
   identical(event_pvalues(tA, tB, eA, eB, alternative = "greater",
                           n_intervals_to_pool = 60),
             event_pvalues(cbind(time, status) ~ arm, data = trial,
                           alternative = "greater", n_intervals_to_pool = 60)))
ok("event_pvalues identical with alternative='both'",
   identical(event_pvalues(tA, tB, eA, eB, n_intervals_to_pool = 60),
             event_pvalues(cbind(time, status) ~ arm, data = trial,
                           n_intervals_to_pool = 60)))
ok("suspected_deviations identical through the formula interface",
   identical(suspected_deviations(tA, tB, eA, eB, n_intervals_to_pool = 60),
             suspected_deviations(cbind(time, status) ~ arm, data = trial,
                                  n_intervals_to_pool = 60)))
ok("suspected_deviations identical with alternative='both'",
   identical(suspected_deviations(tA, tB, eA, eB, alternative = "both",
                                  n_intervals_to_pool = 60),
             suspected_deviations(cbind(time, status) ~ arm, data = trial,
                                  alternative = "both",
                                  n_intervals_to_pool = 60)))

## the plot method returns the same table it would from vectors
pdf(file = tempfile(fileext = ".pdf"))
pv_v <- plot_pvalue_profile(tA, tB, eA, eB, n_intervals_to_pool = 60)
pv_f <- plot_pvalue_profile(cbind(time, status) ~ arm, data = trial,
                            n_intervals_to_pool = 60)
dev.off()
ok("plot_pvalue_profile returns the same table from either interface",
   identical(pv_v, pv_f))

## -------------------------------------------------- group order semantics

g <- higher_criticism_test(cbind(time, status) ~ arm, data = trial,
                           alternative = "greater", n_intervals_to_pool = 60)
trial_rev <- trial
trial_rev$arm <- factor(trial$arm, levels = c("trt", "ctrl"))
l_rev <- higher_criticism_test(cbind(time, status) ~ arm, data = trial_rev,
                               alternative = "less", n_intervals_to_pool = 60)
ok("the first level is group A: reversing the levels swaps greater and less",
   identical(unname(g$statistic), unname(l_rev$statistic)))
ok("reversed levels are reported in the new order",
   identical(l_rev$group_levels, c("trt", "ctrl")))
ok("the planted excess hazard in the second level is found by 'greater'",
   unname(g$statistic) >
     unname(higher_criticism_test(cbind(time, status) ~ arm, data = trial,
                                  alternative = "less",
                                  n_intervals_to_pool = 60)$statistic))

## unordered group variables level by sort(unique(.))
trial_chr <- trial; trial_chr$arm <- as.character(trial$arm)
trial_num <- trial; trial_num$arm <- rep(c(0, 1), c(nA, nB))
trial_lgl <- trial; trial_lgl$arm <- rep(c(FALSE, TRUE), c(nA, nB))
base_stat <- unname(g$statistic)
for (d in list(character = trial_chr, numeric = trial_num, logical = trial_lgl)) {
  ok("grouping by a non-factor gives the same split",
     identical(base_stat,
               unname(higher_criticism_test(cbind(time, status) ~ arm, data = d,
                                            alternative = "greater",
                                            n_intervals_to_pool = 60)$statistic)))
}

## numeric levels sort numerically, not alphabetically
trial_num2 <- trial; trial_num2$arm <- rep(c(2, 10), c(nA, nB))
ok("numeric group levels sort numerically",
   identical(base_stat,
             unname(higher_criticism_test(cbind(time, status) ~ arm,
                                          data = trial_num2,
                                          alternative = "greater",
                                          n_intervals_to_pool = 60)$statistic)))

## ------------------------------------------------- a bare numeric response

ok("'time ~ group' treats every observation as an event",
   identical(unname(higher_criticism_test(tA, tB,
                                          n_intervals_to_pool = 60)$statistic),
             unname(higher_criticism_test(time ~ arm, data = trial,
                                          n_intervals_to_pool = 60)$statistic)))

## --------------------------------------------------- subset and na.action

keepA <- tA <= 15; keepB <- tB <= 15
ok("subset= matches the hand-made subset",
   identical(unname(higher_criticism_test(tA[keepA], tB[keepB],
                                          eA[keepA], eB[keepB],
                                          n_intervals_to_pool = 40)$statistic),
             unname(higher_criticism_test(cbind(time, status) ~ arm,
                                          data = trial, subset = time <= 15,
                                          n_intervals_to_pool = 40)$statistic)))

## na.action is the one documented difference between the interfaces: the
## formula method drops incomplete rows, the default method refuses them
trial_na <- trial
trial_na$time[5L] <- NA            # an observation in group A
trial_na$status[nA + 180L] <- NA   # an observation in group B
dropA <- seq_len(nA) != 5L
dropB <- seq_len(nB) != 180L
ok("na.omit drops exactly the incomplete rows",
   identical(unname(higher_criticism_test(tA[dropA], tB[dropB],
                                          eA[dropA], eB[dropB],
                                          n_intervals_to_pool = 40)$statistic),
             unname(higher_criticism_test(cbind(time, status) ~ arm,
                                          data = trial_na,
                                          n_intervals_to_pool = 40)$statistic)))
ok("na.fail still refuses incomplete rows",
   inherits(try(higher_criticism_test(cbind(time, status) ~ arm,
                                      data = trial_na,
                                      na.action = na.fail), silent = TRUE),
            "try-error"))
ok("the default method still rejects NA durations",
   inherits(try(higher_criticism_test(c(NA, tA[-1]), tB), silent = TRUE),
            "try-error"))

## -------------------------------------------------------------- misuse

ok("rejects a one-sided formula",
   inherits(try(higher_criticism_test(~arm, data = trial), silent = TRUE),
            "try-error"))
ok("rejects two grouping variables",
   inherits(try(higher_criticism_test(cbind(time, status) ~ arm + status,
                                      data = trial), silent = TRUE),
            "try-error"))
ok("rejects a grouping variable with three levels",
   inherits(try(higher_criticism_test(
     cbind(time, status) ~ g3,
     data = transform(trial, g3 = rep(c("a", "b", "c"), length.out = nA + nB))),
     silent = TRUE), "try-error"))
ok("rejects a grouping variable with one level",
   inherits(try(higher_criticism_test(cbind(time, status) ~ g1,
                                      data = transform(trial, g1 = "only")),
                silent = TRUE), "try-error"))
ok("rejects event_observed_A alongside a formula",
   inherits(try(higher_criticism_test(cbind(time, status) ~ arm, data = trial,
                                      event_observed_A = eA), silent = TRUE),
            "try-error"))
ok("rejects a three-column matrix response",
   inherits(try(higher_criticism_test(cbind(time, status, status) ~ arm,
                                      data = trial), silent = TRUE),
            "try-error"))
ok("rejects a Surv object that is not right-censored",
   inherits(try(higher_criticism_test(
     structure(cbind(1, 2, trial$status), type = "counting",
               class = "Surv") ~ arm, data = trial), silent = TRUE),
     "try-error"))
ok("rejects a Surv object passed as durations_A",
   inherits(try(higher_criticism_test(as_surv(tA, eA), tB), silent = TRUE),
            "try-error"))
ok("rejects a matrix passed as durations_A",
   inherits(try(higher_criticism_test(cbind(tA, tA), tB), silent = TRUE),
            "try-error"))
ok("rejects a misspelled argument instead of ignoring it",
   inherits(try(higher_criticism_test(tA, tB, n_intervals_pooled = 40),
                silent = TRUE), "try-error"))
ok("partial matching still reaches n_intervals_to_pool",
   identical(unname(higher_criticism_test(tA, tB, n_intervals = 40)$statistic),
             unname(higher_criticism_test(tA, tB,
                                          n_intervals_to_pool = 40)$statistic)))

## ------------------------------------------------- against real survival::Surv

if (requireNamespace("survival", quietly = TRUE)) {
  sv <- higher_criticism_test(survival::Surv(time, status) ~ arm, data = trial,
                              n_intervals_to_pool = 60, n_permutations = 200,
                              seed = 7)
  ok("survival::Surv agrees with the default method",
     identical(unname(sv$statistic),
               unname(higher_criticism_test(tA, tB, eA, eB,
                                            n_intervals_to_pool = 60,
                                            n_permutations = 200,
                                            seed = 7)$statistic)))
  ok("survival::Surv agrees with the cbind() form",
     identical(unname(sv$statistic),
               unname(higher_criticism_test(cbind(time, status) ~ arm,
                                            data = trial,
                                            n_intervals_to_pool = 60,
                                            n_permutations = 200,
                                            seed = 7)$statistic)))

  ## survival normalises 1/2 status coding to 0/1
  trial12 <- trial; trial12$status12 <- trial$status + 1
  ok("survival::Surv normalises 1/2 status coding",
     identical(unname(sv$statistic),
               unname(higher_criticism_test(
                 survival::Surv(time, status12) ~ arm, data = trial12,
                 n_intervals_to_pool = 60, n_permutations = 200,
                 seed = 7)$statistic)))

  ## na.omit must survive [.Surv
  ok("survival::Surv drops incomplete rows like the cbind() form",
     identical(unname(higher_criticism_test(survival::Surv(time, status) ~ arm,
                                            data = trial_na,
                                            n_intervals_to_pool = 40)$statistic),
               unname(higher_criticism_test(cbind(time, status) ~ arm,
                                            data = trial_na,
                                            n_intervals_to_pool = 40)$statistic)))

  ok("rejects counting-process Surv objects",
     inherits(try(higher_criticism_test(
       survival::Surv(rep(0, nA + nB), time, status) ~ arm, data = trial),
       silent = TRUE), "try-error"))
} else {
  cat("SKIP: 'survival' not installed; Surv() checked against a",
      "structural stand-in only\n")
}

cat("\nAll tests passed.\n")
