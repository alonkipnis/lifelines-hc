## Aggregation statistics for a vector of p-values.
##
## These reproduce the conventions of the Python 'multitest' package used by
## the reference implementation, which differ in detail from the formulae as
## displayed in the paper:
##
##   * the ordered p-values are compared against the beta-normalised means
##     u_i = i / (N + 1) rather than i / N;
##   * the stabilised denominator is sqrt(u_i (1 - u_i) / (N + 2)) rather than
##     sqrt(p_(i) (1 - p_(i)) / N).
##
## Matching these exactly is what makes the statistics agree with the Python
## implementation; see tests/test-parity.R.

## Shared setup: sorted p-values and the standardised deviations z_i.
.mt_prepare <- function(pvals, stbl = TRUE) {
  N <- length(pvals)
  if (N == 0L) stop("'pvals' must be non-empty", call. = FALSE)
  spv <- sort(as.numeric(pvals))
  uu <- seq(1 / (N + 1), 1 - 1 / (N + 1), length.out = N)
  denom <- if (stbl) sqrt(uu * (1 - uu) / (N + 2)) else sqrt(spv * (1 - spv))
  list(N = N, spv = spv, zz = (uu - spv) / denom)
}

## Higher Criticism statistic and the p-value threshold attaining it.
.mt_hc <- function(pvals, gamma = .DEFAULT_GAMMA, stbl = TRUE) {
  mt <- .mt_prepare(pvals, stbl)
  imax <- max(0L, as.integer(gamma * mt$N + 0.5))
  istar <- if (imax <= 0L) 1L else which.max(mt$zz[seq_len(imax)])
  c(score = mt$zz[istar], threshold = mt$spv[istar])
}

## Exact Berk-Jones statistic, returned as -log(BJ) so that large values are
## significant. Follows Moscovich, Nadler and Spiegelman (2016).
.mt_berkjones <- function(pvals, gamma = .DEFAULT_GAMMA) {
  N <- length(pvals)
  if (N == 0L) return(NA_real_)
  spv <- sort(as.numeric(pvals))
  max_i <- max(1L, as.integer(gamma * N))
  s <- spv[seq_len(min(max_i, N))]
  ii <- seq_along(s)
  BJpv <- stats::pbeta(s, ii, N - ii + 1)
  -log(min(min(BJpv), min(1 - BJpv)))
}

## Fisher's combination statistic, -2 sum log p.
.mt_fisher <- function(pvals) sum(-2 * log(as.numeric(pvals)))

## Bonferroni-type statistic, -log(min p).
.mt_minp <- function(pvals) -log(min(as.numeric(pvals)))

.aggregate <- function(pvals, method, gamma, stbl) {
  if (!length(pvals)) return(0)
  switch(method,
    hc         = unname(.mt_hc(pvals, gamma, stbl)[1L]),
    berk_jones = .mt_berkjones(pvals, gamma),
    fisher     = .mt_fisher(pvals),
    min_p      = .mt_minp(pvals),
    stop(sprintf("unknown aggregation method: '%s'", method), call. = FALSE)
  )
}

## Aggregated statistic from survival-table counts, honouring 'alternative'.
.stat_from_counts <- function(Nt1, Nt2, Ot1, Ot2, method, alternative,
                              gamma, stbl) {
  if (!length(Nt1)) return(0)
  if (alternative == "both") {
    g <- .aggregate(.one_direction_pvals(Nt1, Nt2, Ot1, Ot2), method, gamma, stbl)
    l <- .aggregate(.one_direction_pvals(Nt2, Nt1, Ot2, Ot1), method, gamma, stbl)
    return(max(g, l))
  }
  pv <- if (alternative == "greater") {
    .one_direction_pvals(Nt1, Nt2, Ot1, Ot2)
  } else {
    .one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
  }
  .aggregate(pv, method, gamma, stbl)
}
