## Aggregation statistics for a vector of p-values.
##
## Two normalisations of the ordered p-values are supported.
##
## "donoho-jin2008" (default) is the classical Higher Criticism normalisation:
##
##     u_i   = i / N
##     denom = sqrt(u_i (1 - u_i) / N)      (stbl = TRUE)
##     HC_i  = (u_i - p_(i)) / denom
##           = sqrt(N) (i/N - p_(i)) / sqrt((i/N)(1 - i/N))
##
## The largest u_i equals 1, which would make the denominator vanish, so it is
## pulled back by EPS = 1 / (1e4 + N^2) exactly as in the Python 'multitest'
## package. With gamma <= 1 - 1/N the affected index is never the maximiser.
##
## "beta" matches the mean and standard deviation of the ordered p-values:
##
##     u_i   = i / (N + 1)
##     denom = sqrt(u_i (1 - u_i) / (N + 2))
##
## With stbl = FALSE the denominator is sqrt(p_(i) (1 - p_(i))) under either
## normalisation, which drops the sqrt(N) scaling; this mirrors the Python
## implementation rather than correcting it.
##
## Note: the Python package 'multitest' (0.2.1) defaults to "beta", and its
## own "donoho-jin2008" branch builds a constant vector rather than i/N. The
## "beta" option here reproduces that package's default exactly; see
## parity_check.R.

.mt_prepare <- function(pvals, stbl = TRUE,
                        normalization = c("donoho-jin2008", "beta")) {
  normalization <- match.arg(normalization)
  N <- length(pvals)
  if (N == 0L) stop("'pvals' must be non-empty", call. = FALSE)
  spv <- sort(as.numeric(pvals))

  if (normalization == "donoho-jin2008") {
    uu <- seq_len(N) / N
    uu[N] <- uu[N] - 1 / (1e4 + N^2)   # keep the last denominator positive
    std <- sqrt(uu * (1 - uu) / N)
  } else {
    uu <- seq(1 / (N + 1), 1 - 1 / (N + 1), length.out = N)
    std <- sqrt(uu * (1 - uu) / (N + 2))
  }

  denom <- if (stbl) std else sqrt(spv * (1 - spv))
  list(N = N, spv = spv, zz = (uu - spv) / denom)
}

## Higher Criticism statistic and the p-value threshold attaining it.
.mt_hc <- function(pvals, gamma = .DEFAULT_GAMMA, stbl = TRUE,
                   normalization = "donoho-jin2008") {
  mt <- .mt_prepare(pvals, stbl, normalization)
  imax <- max(0L, as.integer(gamma * mt$N + 0.5))
  istar <- if (imax <= 0L) 1L else which.max(mt$zz[seq_len(imax)])
  c(score = mt$zz[istar], threshold = mt$spv[istar])
}

## Exact Berk-Jones statistic, returned as -log(BJ) so that large values are
## significant. Follows Moscovich, Nadler and Spiegelman (2016). Independent
## of the normalisation above.
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

.aggregate <- function(pvals, method, gamma, stbl, normalization) {
  if (!length(pvals)) return(0)
  switch(method,
    hc         = unname(.mt_hc(pvals, gamma, stbl, normalization)[1L]),
    berk_jones = .mt_berkjones(pvals, gamma),
    fisher     = .mt_fisher(pvals),
    min_p      = .mt_minp(pvals),
    stop(sprintf("unknown aggregation method: '%s'", method), call. = FALSE)
  )
}

## Aggregated statistic from survival-table counts, honouring 'alternative'.
.stat_from_counts <- function(Nt1, Nt2, Ot1, Ot2, method, alternative,
                              gamma, stbl, normalization) {
  if (!length(Nt1)) return(0)
  if (alternative == "both") {
    g <- .aggregate(.one_direction_pvals(Nt1, Nt2, Ot1, Ot2),
                    method, gamma, stbl, normalization)
    l <- .aggregate(.one_direction_pvals(Nt2, Nt1, Ot2, Ot1),
                    method, gamma, stbl, normalization)
    return(max(g, l))
  }
  pv <- if (alternative == "greater") {
    .one_direction_pvals(Nt1, Nt2, Ot1, Ot2)
  } else {
    .one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
  }
  .aggregate(pv, method, gamma, stbl, normalization)
}
