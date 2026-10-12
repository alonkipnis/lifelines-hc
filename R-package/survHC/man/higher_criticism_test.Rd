% Generated from R/tests.R
\name{higher_criticism_test}
\alias{higher_criticism_test}
\alias{higher_criticism_test.default}
\alias{higher_criticism_test.formula}
\title{Higher Criticism test for non-proportional hazard deviations}
\usage{
higher_criticism_test(durations_A, ...)

\method{higher_criticism_test}{default}(
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
  seed = NULL,
  ...
)

\method{higher_criticism_test}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.
\code{durations_A} may instead be a \code{\link{formula}}, which selects the
formula method.}

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
into. Strongly recommended for continuous survival data, where each distinct
time carries at most one event and the per-time hypergeometric p-values are too
coarse to be informative. Values between 50 and 200 usually work well.
\code{NULL} (default) applies no pooling.}

\item{n_permutations}{Number of label permutations used to calibrate the
p-value. \code{0} (default) skips calibration and returns \code{NA}.}

\item{seed}{Optional integer seed for the permutation test. When supplied, the
caller's random seed is restored on exit.}

\item{formula}{A two-sided formula \code{Surv(time, status) ~ group}; see
Details.}

\item{data}{An optional data frame, list or environment holding the variables
named in \code{formula}.}

\item{subset}{An optional vector specifying a subset of observations to use.}

\item{na.action}{A function describing what should happen when the data contain
\code{NA}s. Defaults to \code{getOption("na.action")}, normally
\code{\link[stats]{na.omit}}.}

\item{...}{Further arguments. The formula method passes these on to the default
method, which accepts none and signals an error if any are supplied, so that a
misspelled argument name cannot be silently ignored.}
}
\value{
An object of class \code{"htest"}, with the Higher Criticism score in
\code{statistic} and the permutation p-value in \code{p.value}. The permutation
statistics are returned in \code{permutation_statistics}. The formula method
additionally records the two group levels, in the order (A, B), in
\code{group_levels}.
}
\description{
Compares two right-censored survival samples using Higher Criticism applied to
exact per-interval hypergeometric p-values. The test is powerful against
\emph{sparse} departures, where the hazard difference is concentrated in a
small and a priori unknown subset of the follow-up, a regime in which the
log-rank test has little power.
}
\details{
Two interfaces are available and return the same numbers. The default method
takes the two samples as vectors, mirroring the Python package
\code{lifelines-hc}. The formula method takes the usual survival formula
\code{Surv(time, status) ~ group}, which needs the \pkg{survival} package for
\code{Surv} itself; without it, \code{cbind(time, status) ~ group} and
\code{time ~ group} also work, the latter treating every observation as an
event.

The grouping variable must take exactly two values. Its first level is group A
and its second is group B, so \code{alternative = "greater"} tests for excess
hazard in the second level. Levels are taken in factor order for a factor, and
otherwise in the order \code{sort(unique(.))} gives: numerically for numbers,
alphabetically for strings.

The two interfaces differ only in their treatment of missing values. The
default method rejects \code{NA} durations; the formula method applies
\code{na.action}, which by default drops incomplete rows.
}
\examples{
set.seed(1)
T_A <- rexp(200, rate = 1 / 10)
T_B <- rexp(200, rate = 1 / 10)
higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50)

# with permutation calibration
higher_criticism_test(T_A, T_B, n_intervals_to_pool = 50,
                      n_permutations = 200, seed = 42)

# the same data through the formula interface
trial <- data.frame(time = c(T_A, T_B), status = 1,
                    arm = rep(c("A", "B"), c(200, 200)))
higher_criticism_test(cbind(time, status) ~ arm, data = trial,
                      n_intervals_to_pool = 50)

# Surv() needs the survival package, which is only suggested
if (requireNamespace("survival", quietly = TRUE)) {
  higher_criticism_test(survival::Surv(time, status) ~ arm, data = trial,
                        n_intervals_to_pool = 50)
}
}
\references{
Kipnis, A., Galili, B. and Yakhini, Z. (2026). Higher criticism for rare and
weak non-proportional hazard deviations in survival analysis.
\emph{Biometrika} \bold{113}(1), asaf075. \doi{10.1093/biomet/asaf075}

Donoho, D. and Jin, J. (2004). Higher criticism for detecting sparse
heterogeneous mixtures. \emph{The Annals of Statistics} \bold{32}(3),
962--994. \doi{10.1214/009053604000000265}

Kipnis, A., Galili, B. and Yakhini, Z. (2026). \code{lifelines-hc}: Higher
Criticism testing for sparse non-proportional hazard departures in Python.
\doi{10.48550/arXiv.2610.07268}
}
\seealso{
\code{\link{suspected_deviations}} to recover the intervals driving a
detection; \code{\link{berk_jones_test}}, \code{\link{fisher_combination_test}}
and \code{\link{min_p_test}} for alternative aggregations of the same
per-interval p-values.
}
