## Aggregation statistics for a vector of p-values.
##
## The Higher Criticism variants mirror the named methods of the Python
## package 'multitest' (>= 0.2.2), so that hc_version here and hc_version in
## 'lifelines-hc' select the same statistic.
##
##   "dj2008" (default)  u_i = i/N,       denom = sqrt(u_i (1 - u_i) / N)
##   "dj2004"            u_i = i/N,       denom = sqrt(p_(i) (1 - p_(i)) / N)
##   "beta"              u_i = i/(N + 1), denom = sqrt(u_i (1 - u_i) / (N + 2))
##   "star"              beta denominator, but the scan starts at the first
##                       ordered p-value exceeding (1 - EPS)/N
##
## In every case HC_i = (u_i - p_(i)) / denom and the statistic is the maximum
## over the scanned range.
##
## "dj2004" is equation (2) of Kipnis, Galili and Yakhini (2026). "dj2008" is
## the classical Donoho-Jin normalisation and the default throughout.
##
## For "dj2008" and "dj2004" the largest u_i equals 1, which would make the
## denominator vanish, so it is pulled back by EPS = 1 / (1e4 + N^2) exactly as
## in the Python implementation.

.HC_VERSIONS <- c("dj2008", "dj2004", "beta", "star")

.mt_prepare <- function(pvals, hc_version = "dj2008") {
  N <- length(pvals)
  if (N == 0L) stop("'pvals' must be non-empty", call. = FALSE)
  spv <- sort(as.numeric(pvals))
  EPS <- 1 / (1e4 + N^2)

  if (hc_version %in% c("dj2008", "dj2004")) {
    uu <- seq_len(N) / N
    uu[N] <- uu[N] - EPS
    std <- if (hc_version == "dj2008") {
      sqrt(uu * (1 - uu) / N)
    } else {
      sqrt(spv * (1 - spv) / N)
    }
  } else {                       # "beta" and "star" share the standardisation
    uu <- seq(1 / (N + 1), 1 - 1 / (N + 1), length.out = N)
    std <- sqrt(uu * (1 - uu) / (N + 2))
  }

  list(N = N, spv = spv, zz = (uu - spv) / std, EPS = EPS)
}

## Higher Criticism statistic and the p-value threshold attaining it.
##
## Index bookkeeping follows the Python implementation, which scans the
## half-open range [imin, imax) over zero-based indices.
.mt_hc <- function(pvals, gamma = .DEFAULT_GAMMA, hc_version = "dj2008") {
  mt <- .mt_prepare(pvals, hc_version)
  N <- mt$N

  if (hc_version == "star") {
    gt <- mt$spv > (1 - mt$EPS) / N
    imin <- if (any(gt)) which(gt)[1L] - 1L else 0L
    imax <- max(imin + 1L, as.integer(floor(gamma * N + 0.5)))
  } else {
    imin <- 0L
    imax <- max(imin, as.integer(gamma * N + 0.5))
  }
  imax <- min(imax, N)

  if (imin > imax) return(c(score = NA_real_, threshold = NA_real_))
  istar <- if (imin == imax) {
    imin
  } else {
    imin + which.max(mt$zz[(imin + 1L):imax]) - 1L
  }

  i <- istar + 1L
  c(score = mt$zz[i], threshold = mt$spv[i])
}

## Exact Berk-Jones statistic, returned as -log(BJ) so that large values are
## significant. Follows Moscovich, Nadler and Spiegelman (2016). Independent
## of hc_version.
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

.aggregate <- function(pvals, method, gamma, hc_version) {
  if (!length(pvals)) return(0)
  switch(method,
    hc         = unname(.mt_hc(pvals, gamma, hc_version)[1L]),
    berk_jones = .mt_berkjones(pvals, gamma),
    fisher     = .mt_fisher(pvals),
    min_p      = .mt_minp(pvals),
    stop(sprintf("unknown aggregation method: '%s'", method), call. = FALSE)
  )
}

## Aggregated statistic from survival-table counts, honouring 'alternative'.
.stat_from_counts <- function(Nt1, Nt2, Ot1, Ot2, method, alternative,
                              gamma, hc_version) {
  if (!length(Nt1)) return(0)
  if (alternative == "both") {
    g <- .aggregate(.one_direction_pvals(Nt1, Nt2, Ot1, Ot2),
                    method, gamma, hc_version)
    l <- .aggregate(.one_direction_pvals(Nt2, Nt1, Ot2, Ot1),
                    method, gamma, hc_version)
    return(max(g, l))
  }
  pv <- if (alternative == "greater") {
    .one_direction_pvals(Nt1, Nt2, Ot1, Ot2)
  } else {
    .one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
  }
  .aggregate(pv, method, gamma, hc_version)
}
