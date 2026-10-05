% Generated from R/tests.R
\name{berk_jones_test}
\alias{berk_jones_test}
\title{Berk-Jones test for non-proportional hazard deviations}
\usage{
berk_jones_test(
  durations_A,
  durations_B,
  event_observed_A = NULL,
  event_observed_B = NULL,
  alternative = c("both", "greater", "less"),
  gamma = 0.2,
  stbl = TRUE,
  t_0 = -1,
  n_intervals_to_pool = NULL,
  n_permutations = 0,
  seed = NULL
)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.}

\item{event_observed_A, event_observed_B}{Optional 0/1 vectors, 1 indicating an
observed event and 0 a censored observation. Default: all events.}

\item{alternative}{One of \code{"both"} (default), \code{"greater"} (excess
hazard in group B) or \code{"less"} (excess hazard in group A).}

\item{gamma}{Fraction of the smallest ordered p-values scanned, in (0, 1].
Default 0.2.}

\item{stbl}{Logical; use the variance-stabilised denominator. Default
\code{TRUE}.}

\item{t_0}{Restrict the analysis to events at or before \code{t_0}. A negative
value (default) applies no restriction.}

\item{n_intervals_to_pool}{Number of equal-width intervals to pool event times
into. Strongly recommended for continuous survival data. Values between 50 and
200 usually work well. \code{NULL} (default) applies no pooling.}

\item{n_permutations}{Number of label permutations used to calibrate the
p-value. \code{0} (default) skips calibration and returns \code{NA}.}

\item{seed}{Optional integer seed for the permutation test. When supplied, the
caller's random seed is restored on exit.}
}
\value{
An object of class \code{"htest"}.
}
\description{
Aggregates the same per-interval hypergeometric p-values as \\code{\\link{higher_criticism_test}} using the exact Berk-Jones statistic, returned as \\eqn{-\\log(BJ)} so that large values are significant.

The statistic is \code{Inf} when the exact Berk-Jones p-value underflows to
zero, which happens routinely on unpooled continuous survival data where the
smallest per-interval p-value is extreme. Set \code{n_intervals_to_pool} to
obtain a finite value. This matches the behaviour of the reference Python
implementation.
}
\examples{
set.seed(1)
berk_jones_test(rexp(200, 1 / 10), rexp(200, 1 / 10), n_intervals_to_pool = 50)
}
\seealso{
\code{\link{higher_criticism_test}}
}
