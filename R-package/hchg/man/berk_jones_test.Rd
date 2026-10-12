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
  hc_version = c("dj2008", "dj2004", "beta", "star"),
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

\item{hc_version}{Which Higher Criticism variant to use, matching the named
methods of the Python package \code{multitest} (>= 0.2.2) and the
\code{hc_version} argument of \code{lifelines-hc}. \code{"dj2008"} (default)
uses \eqn{u_i = i/N} with denominator \eqn{\sqrt{u_i(1-u_i)/N}};
\code{"dj2004"} uses \eqn{u_i = i/N} with denominator
\eqn{\sqrt{p_{(i)}(1-p_{(i)})/N}}, which is equation (2) of Kipnis, Galili and
Yakhini (2026); \code{"beta"} uses \eqn{u_i = i/(N+1)} with denominator
\eqn{\sqrt{u_i(1-u_i)/(N+2)}}; \code{"star"} uses the \code{"beta"}
denominator but starts the scan at the first ordered p-value exceeding
\eqn{(1-\epsilon)/N}. Ignored by \code{berk_jones_test},
\code{fisher_combination_test} and \code{min_p_test}.}

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
\references{
Moscovich, A., Nadler, B. and Spiegelman, C. (2016). On the exact Berk-Jones
statistics and their p-value calculation. \emph{Electronic Journal of
Statistics} \bold{10}(2), 2329--2354.

Kipnis, A., Galili, B. and Yakhini, Z. (2026). Higher criticism for rare and
weak non-proportional hazard deviations in survival analysis.
\emph{Biometrika} \bold{113}(1), asaf075. \doi{10.1093/biomet/asaf075}
}
\seealso{
\code{\link{higher_criticism_test}}
}
