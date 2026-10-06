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
  stbl = TRUE,
  normalization = c("donoho-jin2008", "beta"),
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

\item{stbl}{Logical; use the variance-stabilised denominator. Default
\code{TRUE}.}

\item{normalization}{Normalisation of the ordered p-values used by the Higher
Criticism statistic. \code{"donoho-jin2008"} (default) uses \eqn{u_i = i/N}
with denominator \eqn{\sqrt{u_i(1-u_i)/N}}; \code{"beta"} uses
\eqn{u_i = i/(N+1)} with denominator \eqn{\sqrt{u_i(1-u_i)/(N+2)}}, which is
what the Python package \code{multitest} applies by default. The choice changes
the value of the statistic and therefore the permutation p-value; it is ignored
by \code{berk_jones_test}, \code{fisher_combination_test} and
\code{min_p_test}.}

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
\seealso{
\code{\link{higher_criticism_test}}, \code{\link{plot_pvalue_profile}}
}
