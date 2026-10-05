% Generated from R/tests.R
\name{higher_criticism_test}
\alias{higher_criticism_test}
\title{Higher Criticism test for non-proportional hazard deviations}
\usage{
higher_criticism_test(
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
into. Strongly recommended for continuous survival data, where each distinct
time carries at most one event and the per-time hypergeometric p-values are too
coarse to be informative. Values between 50 and 200 usually work well.
\code{NULL} (default) applies no pooling.}

\item{n_permutations}{Number of label permutations used to calibrate the
p-value. \code{0} (default) skips calibration and returns \code{NA}.}

\item{seed}{Optional integer seed for the permutation test. When supplied, the
caller's random seed is restored on exit.}
}
\value{
An object of class \code{"htest"}, with the Higher Criticism score in
\code{statistic} and the permutation p-value in \code{p.value}. The permutation
statistics are returned in \code{permutation_statistics}.
}
\description{
Compares two right-censored survival samples using Higher Criticism applied to
exact per-interval hypergeometric p-values. The test is powerful against
\emph{sparse} departures, where the hazard difference is concentrated in a
small and a priori unknown subset of the follow-up, a regime in which the
log-rank test has little power.
}
\examples{
set.seed(1)
T_A <- rexp(200, rate = 1 / 10)
T_B <- rexp(200, rate = 1 / 10)
higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50)

# with permutation calibration
higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50,
                      n_permutations = 200, seed = 42)
}
\references{
Kipnis, A., Galili, B. and Yakhini, Z. (2026). Higher criticism for rare and
weak non-proportional hazard deviations in survival analysis.
\emph{Biometrika} \bold{113}(1), asaf075. \doi{10.1093/biomet/asaf075}

Donoho, D. and Jin, J. (2004). Higher criticism for detecting sparse
heterogeneous mixtures. \emph{The Annals of Statistics} \bold{32}(3),
962--994. \doi{10.1214/009053604000000265}
}
\seealso{
\code{\link{suspected_deviations}} to recover the intervals driving a
detection; \code{\link{berk_jones_test}}, \code{\link{fisher_combination_test}}
and \code{\link{min_p_test}} for alternative aggregations of the same
per-interval p-values.
}
