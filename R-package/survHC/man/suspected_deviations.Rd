% Generated from R/tests.R
\name{suspected_deviations}
\alias{suspected_deviations}
\title{Intervals with suspected non-proportional hazard deviations}
\usage{
suspected_deviations(
  durations_A,
  durations_B,
  event_observed_A = NULL,
  event_observed_B = NULL,
  alternative = c("greater", "less", "both"),
  gamma = 0.2,
  hc_version = c("dj2008", "dj2004", "beta", "star"),
  t_0 = -1,
  n_intervals_to_pool = NULL
)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.}

\item{event_observed_A, event_observed_B}{Optional 0/1 vectors, 1 indicating an
observed event and 0 a censored observation. Default: all events.}

\item{alternative}{One of \code{"greater"} (default), \code{"less"} or
\code{"both"}. With \code{"both"}, columns for the reverse direction are added
and \code{suspected} flags either direction.}

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
into. \code{NULL} (default) applies no pooling.}
}
\value{
A data frame with one row per event time or pooled interval and columns
\code{time} (or \code{time_lower} and \code{time_upper} when pooled),
\code{at_risk_A}, \code{at_risk_B}, \code{observed_A}, \code{observed_B},
\code{hypergeom_pvalue}, \code{suspected} and \code{hc_threshold}.
}
\description{
Computes per-interval hypergeometric p-values and flags those at or below the
Higher Criticism threshold. This localises a detection to specific windows of
the follow-up, which is the diagnostic information a single-number test cannot
provide, and is suitable for shading on a Kaplan-Meier plot.
}
\examples{
set.seed(1)
# a sparse departure confined to one window
T_A <- rexp(300, 1 / 10)
T_B <- c(rexp(250, 1 / 10), rexp(50, 1 / 2))
d <- suspected_deviations(T_A, T_B, n_intervals_to_pool = 50)
head(d[d$suspected, ])
}
\references{
Kipnis, A., Galili, B. and Yakhini, Z. (2026). Higher criticism for rare and
weak non-proportional hazard deviations in survival analysis.
\emph{Biometrika} \bold{113}(1), asaf075. \doi{10.1093/biomet/asaf075}

Kipnis, A., Galili, B. and Yakhini, Z. (2026). \code{lifelines-hc}: Higher
Criticism testing for sparse non-proportional hazard departures in Python.
\doi{10.48550/arXiv.2610.07268}
}
\seealso{
\code{\link{higher_criticism_test}}, \code{\link{plot_pvalue_profile}}
}
