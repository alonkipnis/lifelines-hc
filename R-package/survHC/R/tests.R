## Public API: the four omnibus tests, the per-interval p-values, and the
## interval-level diagnostic table.
##
## Each is an S3 generic with a "default" method taking the two samples as
## vectors and a "formula" method taking Surv(time, status) ~ group. The
## formula methods do no statistics of their own; they split the data and call
## the default method. See R/formula.R.

.permutation_pvalue <- function(durations_A, durations_B,
                                event_observed_A, event_observed_B,
                                observed_stat, n_permutations, method,
                                alternative, gamma, hc_version, t_0,
                                n_intervals_to_pool) {
  nA <- length(durations_A)
  all_dur <- c(as.numeric(durations_A), as.numeric(durations_B))
  all_evt <- c(
    if (is.null(event_observed_A)) rep(1, nA) else as.numeric(event_observed_A),
    if (is.null(event_observed_B)) rep(1, length(durations_B)) else as.numeric(event_observed_B)
  )
  n <- length(all_dur)

  perm_stats <- numeric(n_permutations)
  for (i in seq_len(n_permutations)) {
    perm <- sample.int(n)
    iA <- perm[seq_len(nA)]
    iB <- perm[(nA + 1L):n]
    cnt <- .survival_table_counts(all_dur[iA], all_dur[iB],
                                  all_evt[iA], all_evt[iB],
                                  t_0, n_intervals_to_pool)
    perm_stats[i] <- .stat_from_counts(cnt$Nt1, cnt$Nt2, cnt$Ot1, cnt$Ot2,
                                       method, alternative, gamma,
                                       hc_version)
  }
  list(p_value = (sum(perm_stats >= observed_stat) + 1) / (n_permutations + 1),
       perm_stats = perm_stats)
}

.run_test <- function(durations_A, durations_B, event_observed_A,
                      event_observed_B, method, test_name, alternative,
                      gamma, hc_version, t_0, n_intervals_to_pool,
                      n_permutations, seed, data_name) {
  alternative <- match.arg(alternative, c("both", "greater", "less"))
  hc_version <- match.arg(hc_version, .HC_VERSIONS)
  if (gamma <= 0 || gamma > 1) {
    stop("'gamma' must lie in (0, 1]", call. = FALSE)
  }
  if (n_permutations < 0) {
    stop("'n_permutations' must be non-negative", call. = FALSE)
  }

  cnt <- .survival_table_counts(durations_A, durations_B,
                                event_observed_A, event_observed_B,
                                t_0, n_intervals_to_pool)
  stat <- .stat_from_counts(cnt$Nt1, cnt$Nt2, cnt$Ot1, cnt$Ot2,
                            method, alternative, gamma, hc_version)

  p_value <- NA_real_
  perm_stats <- NULL
  if (n_permutations > 0) {
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = globalenv()), add = TRUE)
      }
      set.seed(seed)
    }
    perm <- .permutation_pvalue(durations_A, durations_B,
                                event_observed_A, event_observed_B,
                                stat, n_permutations, method, alternative,
                                gamma, hc_version, t_0,
                                n_intervals_to_pool)
    p_value <- perm$p_value
    perm_stats <- perm$perm_stats
  }

  statistic <- stat
  names(statistic) <- switch(method,
    hc = "HC", berk_jones = "-log(BJ)", fisher = "Fisher", min_p = "-log(min p)")

  parameter <- c(gamma = gamma, n_intervals = length(cnt$Nt1),
                 n_permutations = n_permutations)

  structure(list(
    statistic = statistic,
    parameter = parameter,
    p.value = p_value,
    alternative = alternative,
    method = test_name,
    data.name = data_name,
    permutation_statistics = perm_stats,
    n_intervals_to_pool = n_intervals_to_pool,
    hc_version = hc_version
  ), class = "htest")
}

#' Higher Criticism test for non-proportional hazard deviations
#'
#' Compares two right-censored survival samples using Higher Criticism applied
#' to exact per-interval hypergeometric p-values. The test is powerful against
#' \emph{sparse} departures, where the hazard difference is concentrated in a
#' small and a priori unknown subset of the follow-up, a regime in which the
#' log-rank test has little power.
#'
#' Two interfaces are available and return the same numbers. The default
#' method takes the two samples as vectors, mirroring the Python package
#' \code{lifelines-hc}. The formula method takes the usual survival formula
#' \code{Surv(time, status) ~ group}, which needs the \pkg{survival} package
#' for \code{Surv} itself; without it, \code{cbind(time, status) ~ group} and
#' \code{time ~ group} also work, the latter treating every observation as an
#' event.
#'
#' The grouping variable must take exactly two values. Its first level is
#' group A and its second is group B, so \code{alternative = "greater"} tests
#' for excess hazard in the second level. Levels are taken in factor order for
#' a factor, and otherwise in the order \code{sort(unique(.))} gives:
#' numerically for numbers, alphabetically for strings.
#'
#' The two interfaces differ only in their treatment of missing values. The
#' default method rejects \code{NA} durations; the formula method applies
#' \code{na.action}, which by default drops incomplete rows.
#'
#' @param durations_A,durations_B Numeric vectors of event or censoring times.
#'   \code{durations_A} may instead be a \code{\link{formula}}, which selects
#'   the formula method.
#' @param event_observed_A,event_observed_B Optional 0/1 vectors, 1 indicating
#'   an observed event and 0 a censored observation. Default: all events.
#' @param formula A two-sided formula \code{Surv(time, status) ~ group}; see
#'   Details.
#' @param data An optional data frame, list or environment holding the
#'   variables named in \code{formula}.
#' @param subset An optional vector specifying a subset of observations.
#' @param na.action A function describing what should happen when the data
#'   contain \code{NA}s. Defaults to \code{getOption("na.action")}.
#' @param alternative One of \code{"both"} (default), \code{"greater"} (excess
#'   hazard in group B) or \code{"less"} (excess hazard in group A).
#' @param gamma Fraction of the smallest ordered p-values scanned, in (0, 1].
#'   Default 0.2.
#' @param hc_version Which Higher Criticism variant to use, matching the named
#'   methods of the Python package \code{multitest} (>= 0.2.2) and the
#'   \code{hc_version} argument of \code{lifelines-hc}:
#'   \code{"dj2008"} (default) uses \eqn{u_i = i/N} with denominator
#'   \eqn{\sqrt{u_i(1-u_i)/N}}; \code{"dj2004"} uses \eqn{u_i = i/N} with
#'   denominator \eqn{\sqrt{p_{(i)}(1-p_{(i)})/N}}, which is equation (2) of
#'   Kipnis, Galili and Yakhini (2026); \code{"beta"} uses
#'   \eqn{u_i = i/(N+1)} with denominator \eqn{\sqrt{u_i(1-u_i)/(N+2)}}; and
#'   \code{"star"} uses the \code{"beta"} denominator but starts the scan at
#'   the first ordered p-value exceeding \eqn{(1-\epsilon)/N}. Ignored by
#'   \code{berk_jones_test}, \code{fisher_combination_test} and
#'   \code{min_p_test}.
#' @param t_0 Restrict the analysis to events at or before \code{t_0}. A
#'   negative value (default) applies no restriction.
#' @param n_intervals_to_pool Number of equal-width intervals to pool event
#'   times into. Strongly recommended for continuous survival data, where each
#'   distinct time carries at most one event and the per-time hypergeometric
#'   p-values are too coarse to be informative. Values between 50 and 200
#'   usually work well. \code{NULL} (default) applies no pooling.
#' @param n_permutations Number of label permutations used to calibrate the
#'   p-value. \code{0} (default) skips calibration and returns \code{NA}.
#' @param seed Optional integer seed for the permutation test. When supplied,
#'   the caller's random seed is restored on exit.
#' @param ... Further arguments passed to the default method.
#'
#' @return An object of class \code{"htest"}, with the Higher Criticism score
#'   in \code{statistic} and the permutation p-value in \code{p.value}. The
#'   permutation statistics are returned in \code{permutation_statistics}. The
#'   formula method additionally records the two group levels, in the order
#'   (A, B), in \code{group_levels}.
#'
#' @references
#' Kipnis, A., Galili, B. and Yakhini, Z. (2026). Higher criticism for rare and
#' weak non-proportional hazard deviations in survival analysis.
#' \emph{Biometrika} \bold{113}(1), asaf075. \doi{10.1093/biomet/asaf075}
#'
#' Donoho, D. and Jin, J. (2004). Higher criticism for detecting sparse
#' heterogeneous mixtures. \emph{The Annals of Statistics} \bold{32}(3),
#' 962--994. \doi{10.1214/009053604000000265}
#'
#' @seealso \code{\link{suspected_deviations}} to recover the intervals driving
#'   a detection; \code{\link{berk_jones_test}},
#'   \code{\link{fisher_combination_test}} and \code{\link{min_p_test}} for
#'   alternative aggregations of the same per-interval p-values.
#'
#' @examples
#' set.seed(1)
#' T_A <- rexp(200, rate = 1 / 10)
#' T_B <- rexp(200, rate = 1 / 10)
#' higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50)
#'
#' # with permutation calibration
#' higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50,
#'                       n_permutations = 200, seed = 42)
#'
#' # the same data through the formula interface
#' trial <- data.frame(time = c(T_A, T_B), status = 1,
#'                     arm = rep(c("A", "B"), c(200, 200)))
#' higher_criticism_test(cbind(time, status) ~ arm, data = trial,
#'                       n_intervals_to_pool = 50)
#'
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   higher_criticism_test(survival::Surv(time, status) ~ arm, data = trial,
#'                         n_intervals_to_pool = 50)
#' }
#' @export
higher_criticism_test <- function(durations_A, ...) {
  UseMethod("higher_criticism_test")
}

#' @rdname higher_criticism_test
#' @export
higher_criticism_test.default <- function(durations_A, durations_B,
                                          event_observed_A = NULL,
                                          event_observed_B = NULL,
                                          alternative = c("both", "greater", "less"),
                                          gamma = 0.2,
                                          hc_version = c("dj2008", "dj2004", "beta", "star"),
                                          t_0 = -1,
                                          n_intervals_to_pool = NULL,
                                          n_permutations = 0, seed = NULL, ...) {
  .reject_unused(...)
  dn <- .data_name(substitute(durations_A), substitute(durations_B))
  .run_test(durations_A, durations_B, event_observed_A, event_observed_B,
            method = "hc",
            test_name = "Higher Criticism test for non-proportional hazards",
            alternative = alternative, gamma = gamma, hc_version = hc_version, t_0 = t_0,
            n_intervals_to_pool = n_intervals_to_pool,
            n_permutations = n_permutations, seed = seed, data_name = dn)
}

#' @rdname higher_criticism_test
#' @export
higher_criticism_test.formula <- function(formula, data, subset, na.action, ...) {
  .formula_htest(higher_criticism_test.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}

#' Berk-Jones test for non-proportional hazard deviations
#'
#' Aggregates the same per-interval hypergeometric p-values as
#' \code{\link{higher_criticism_test}} using the exact Berk-Jones statistic,
#' returned as \eqn{-\log(BJ)} so that large values are significant.
#'
#' The statistic is \code{Inf} when the exact Berk-Jones p-value underflows to
#' zero, which happens routinely on unpooled continuous survival data where the
#' smallest per-interval p-value is extreme. Set \code{n_intervals_to_pool} to
#' obtain a finite value. This matches the behaviour of the reference Python
#' implementation.
#'
#' @inheritParams higher_criticism_test
#' @return An object of class \code{"htest"}.
#' @references
#' Moscovich, A., Nadler, B. and Spiegelman, C. (2016). On the exact Berk-Jones
#' statistics and their p-value calculation. \emph{Electronic Journal of
#' Statistics} \bold{10}(2), 2329--2354.
#' @seealso \code{\link{higher_criticism_test}}
#' @examples
#' set.seed(1)
#' berk_jones_test(rexp(200, 1 / 10), rexp(200, 1 / 10),
#'                 n_intervals_to_pool = 50)
#' @export
berk_jones_test <- function(durations_A, ...) {
  UseMethod("berk_jones_test")
}

#' @rdname berk_jones_test
#' @export
berk_jones_test.default <- function(durations_A, durations_B,
                                    event_observed_A = NULL, event_observed_B = NULL,
                                    alternative = c("both", "greater", "less"),
                                    gamma = 0.2,
                                    hc_version = c("dj2008", "dj2004", "beta", "star"),
                                    t_0 = -1,
                                    n_intervals_to_pool = NULL,
                                    n_permutations = 0, seed = NULL, ...) {
  .reject_unused(...)
  dn <- .data_name(substitute(durations_A), substitute(durations_B))
  .run_test(durations_A, durations_B, event_observed_A, event_observed_B,
            method = "berk_jones",
            test_name = "Berk-Jones test for non-proportional hazards",
            alternative = alternative, gamma = gamma, hc_version = hc_version, t_0 = t_0,
            n_intervals_to_pool = n_intervals_to_pool,
            n_permutations = n_permutations, seed = seed, data_name = dn)
}

#' @rdname berk_jones_test
#' @export
berk_jones_test.formula <- function(formula, data, subset, na.action, ...) {
  .formula_htest(berk_jones_test.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}

#' Fisher combination test for non-proportional hazard deviations
#'
#' Aggregates the same per-interval hypergeometric p-values as
#' \code{\link{higher_criticism_test}} using Fisher's method,
#' \eqn{-2 \sum \log p_i}.
#'
#' @inheritParams higher_criticism_test
#' @return An object of class \code{"htest"}.
#' @seealso \code{\link{higher_criticism_test}}
#' @examples
#' set.seed(1)
#' fisher_combination_test(rexp(200, 1 / 10), rexp(200, 1 / 10),
#'                         n_intervals_to_pool = 50)
#' @export
fisher_combination_test <- function(durations_A, ...) {
  UseMethod("fisher_combination_test")
}

#' @rdname fisher_combination_test
#' @export
fisher_combination_test.default <- function(durations_A, durations_B,
                                            event_observed_A = NULL,
                                            event_observed_B = NULL,
                                            alternative = c("both", "greater", "less"),
                                            gamma = 0.2,
                                            hc_version = c("dj2008", "dj2004", "beta", "star"),
                                            t_0 = -1,
                                            n_intervals_to_pool = NULL,
                                            n_permutations = 0, seed = NULL, ...) {
  .reject_unused(...)
  dn <- .data_name(substitute(durations_A), substitute(durations_B))
  .run_test(durations_A, durations_B, event_observed_A, event_observed_B,
            method = "fisher",
            test_name = "Fisher combination test for non-proportional hazards",
            alternative = alternative, gamma = gamma, hc_version = hc_version, t_0 = t_0,
            n_intervals_to_pool = n_intervals_to_pool,
            n_permutations = n_permutations, seed = seed, data_name = dn)
}

#' @rdname fisher_combination_test
#' @export
fisher_combination_test.formula <- function(formula, data, subset, na.action, ...) {
  .formula_htest(fisher_combination_test.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}

#' Minimum-p test for non-proportional hazard deviations
#'
#' Aggregates the same per-interval hypergeometric p-values as
#' \code{\link{higher_criticism_test}} by taking the minimum, reported as
#' \eqn{-\log(\min_i p_i)}.
#'
#' @inheritParams higher_criticism_test
#' @return An object of class \code{"htest"}.
#' @seealso \code{\link{higher_criticism_test}}
#' @examples
#' set.seed(1)
#' min_p_test(rexp(200, 1 / 10), rexp(200, 1 / 10), n_intervals_to_pool = 50)
#' @export
min_p_test <- function(durations_A, ...) {
  UseMethod("min_p_test")
}

#' @rdname min_p_test
#' @export
min_p_test.default <- function(durations_A, durations_B,
                               event_observed_A = NULL, event_observed_B = NULL,
                               alternative = c("both", "greater", "less"),
                               gamma = 0.2,
                               hc_version = c("dj2008", "dj2004", "beta", "star"),
                               t_0 = -1, n_intervals_to_pool = NULL,
                               n_permutations = 0, seed = NULL, ...) {
  .reject_unused(...)
  dn <- .data_name(substitute(durations_A), substitute(durations_B))
  .run_test(durations_A, durations_B, event_observed_A, event_observed_B,
            method = "min_p",
            test_name = "Minimum-p test for non-proportional hazards",
            alternative = alternative, gamma = gamma, hc_version = hc_version, t_0 = t_0,
            n_intervals_to_pool = n_intervals_to_pool,
            n_permutations = n_permutations, seed = seed, data_name = dn)
}

#' @rdname min_p_test
#' @export
min_p_test.formula <- function(formula, data, subset, na.action, ...) {
  .formula_htest(min_p_test.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}

#' Per-interval hypergeometric p-values
#'
#' Returns the exact one-sided hypergeometric p-value for each event time, or
#' for each pooled interval when \code{n_intervals_to_pool} is supplied. Each
#' p-value quantifies the evidence against equal hazards at that point in the
#' follow-up.
#'
#' @inheritParams higher_criticism_test
#' @param alternative One of \code{"both"} (default), returning a list with
#'   components \code{greater} and \code{less}; \code{"greater"}; or
#'   \code{"less"}.
#' @return A numeric vector, or a list of two numeric vectors when
#'   \code{alternative = "both"}.
#' @seealso \code{\link{suspected_deviations}}
#' @examples
#' set.seed(1)
#' p <- event_pvalues(rexp(200, 1 / 10), rexp(200, 1 / 10),
#'                    alternative = "greater", n_intervals_to_pool = 50)
#' summary(p)
#' @export
event_pvalues <- function(durations_A, ...) {
  UseMethod("event_pvalues")
}

#' @rdname event_pvalues
#' @export
event_pvalues.default <- function(durations_A, durations_B,
                                  event_observed_A = NULL, event_observed_B = NULL,
                                  alternative = c("both", "greater", "less"),
                                  t_0 = -1, n_intervals_to_pool = NULL, ...) {
  .reject_unused(...)
  alternative <- match.arg(alternative)
  cnt <- .survival_table_counts(durations_A, durations_B,
                                event_observed_A, event_observed_B,
                                t_0, n_intervals_to_pool)
  if (alternative == "both") {
    return(list(
      greater = .one_direction_pvals(cnt$Nt1, cnt$Nt2, cnt$Ot1, cnt$Ot2),
      less    = .one_direction_pvals(cnt$Nt2, cnt$Nt1, cnt$Ot2, cnt$Ot1)
    ))
  }
  if (alternative == "greater") {
    .one_direction_pvals(cnt$Nt1, cnt$Nt2, cnt$Ot1, cnt$Ot2)
  } else {
    .one_direction_pvals(cnt$Nt2, cnt$Nt1, cnt$Ot2, cnt$Ot1)
  }
}

#' @rdname event_pvalues
#' @export
event_pvalues.formula <- function(formula, data, subset, na.action, ...) {
  .formula_plain(event_pvalues.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}

#' Intervals with suspected non-proportional hazard deviations
#'
#' Computes per-interval hypergeometric p-values and flags those at or below
#' the Higher Criticism threshold. This localises a detection to specific
#' windows of the follow-up, which is the diagnostic information a
#' single-number test cannot provide, and is suitable for shading on a
#' Kaplan-Meier plot.
#'
#' @inheritParams higher_criticism_test
#' @param alternative One of \code{"greater"} (default), \code{"less"} or
#'   \code{"both"}. With \code{"both"}, columns for the reverse direction are
#'   added and \code{suspected} flags either direction.
#'
#' @return A data frame with one row per event time or pooled interval and
#'   columns \code{time} (or \code{time_lower} and \code{time_upper} when
#'   pooled), \code{at_risk_A}, \code{at_risk_B}, \code{observed_A},
#'   \code{observed_B}, \code{hypergeom_pvalue}, \code{suspected} and
#'   \code{hc_threshold}.
#'
#' @seealso \code{\link{higher_criticism_test}},
#'   \code{\link{plot_pvalue_profile}}
#' @examples
#' set.seed(1)
#' # a sparse departure confined to one window
#' T_A <- rexp(300, 1 / 10)
#' T_B <- c(rexp(250, 1 / 10), rexp(50, 1 / 2))
#' d <- suspected_deviations(T_A, T_B, n_intervals_to_pool = 50)
#' head(d[d$suspected, ])
#' @export
suspected_deviations <- function(durations_A, ...) {
  UseMethod("suspected_deviations")
}

#' @rdname suspected_deviations
#' @export
suspected_deviations.default <- function(durations_A, durations_B,
                                         event_observed_A = NULL,
                                         event_observed_B = NULL,
                                         alternative = c("greater", "less", "both"),
                                         gamma = 0.2,
                                         hc_version = c("dj2008", "dj2004", "beta", "star"),
                                         t_0 = -1,
                                         n_intervals_to_pool = NULL, ...) {
  .reject_unused(...)
  alternative <- match.arg(alternative)
  hc_version <- match.arg(hc_version)
  cnt <- .survival_table_counts(durations_A, durations_B,
                                event_observed_A, event_observed_B,
                                t_0, n_intervals_to_pool)

  pvals <- if (alternative == "less") {
    .one_direction_pvals(cnt$Nt2, cnt$Nt1, cnt$Ot2, cnt$Ot1)
  } else {
    .one_direction_pvals(cnt$Nt1, cnt$Nt2, cnt$Ot1, cnt$Ot2)
  }

  usable <- pvals <= 1
  hc_thresh <- if (any(usable)) {
    unname(.mt_hc(pvals[usable], gamma, hc_version)[2L])
  } else {
    0
  }
  flagged <- pvals <= hc_thresh

  out <- data.frame(
    at_risk_A = cnt$Nt1,
    at_risk_B = cnt$Nt2,
    observed_A = cnt$Ot1,
    observed_B = cnt$Ot2,
    hypergeom_pvalue = pvals,
    suspected = flagged,
    hc_threshold = hc_thresh,
    stringsAsFactors = FALSE
  )

  if (is.null(cnt$bin_width)) {
    out <- cbind(time = round(cnt$times, 4), out)
  } else {
    half <- cnt$bin_width / 2
    out <- cbind(time_lower = cnt$times - half,
                 time_upper = cnt$times + half, out)
  }

  if (alternative == "both") {
    pv_rev <- .one_direction_pvals(cnt$Nt2, cnt$Nt1, cnt$Ot2, cnt$Ot1)
    usable_rev <- pv_rev <= 1
    thr_rev <- if (any(usable_rev)) {
      unname(.mt_hc(pv_rev[usable_rev], gamma, hc_version)[2L])
    } else {
      0
    }
    flagged_rev <- pv_rev <= thr_rev
    out$hypergeom_pvalue_rev <- pv_rev
    out$suspected_rev <- flagged_rev
    out$hc_threshold_rev <- thr_rev
    out$suspected <- flagged | flagged_rev
  }

  rownames(out) <- NULL
  out
}

#' @rdname suspected_deviations
#' @export
suspected_deviations.formula <- function(formula, data, subset, na.action, ...) {
  .formula_plain(suspected_deviations.default,
                 .formula_split(formula, match.call(expand.dots = FALSE),
                                parent.frame(), list(...)),
                 list(...))
}
