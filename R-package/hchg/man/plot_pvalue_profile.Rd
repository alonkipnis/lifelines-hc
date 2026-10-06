% Generated from R/plot.R
\name{plot_pvalue_profile}
\alias{plot_pvalue_profile}
\title{Plot the per-interval p-value profile}
\usage{
plot_pvalue_profile(
  durations_A,
  durations_B,
  event_observed_A = NULL,
  event_observed_B = NULL,
  gamma = 0.2,
  stbl = TRUE,
  normalization = c("donoho-jin2008", "beta"),
  t_0 = -1,
  n_intervals_to_pool = NULL,
  col_flagged = "firebrick",
  col_plain = "steelblue",
  xlab = "Time",
  ylab = "signed -log10(p)",
  main = NULL,
  ...
)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.}

\item{event_observed_A, event_observed_B}{Optional 0/1 vectors, 1 indicating an
observed event and 0 a censored observation. Default: all events.}

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

\item{col_flagged}{Colour for intervals at or below the HC threshold.}

\item{col_plain}{Colour for the remaining intervals.}

\item{xlab, ylab, main}{Standard graphical parameters.}

\item{...}{Further arguments passed to \code{\link[graphics]{plot.default}}.}
}
\value{
Invisibly, the data frame returned by \code{\link{suspected_deviations}} with
\code{alternative = "both"}.
}
\description{
Draws the signed \eqn{-\log_{10}} per-interval hypergeometric p-value profile,
with bars above the axis marking excess events in group B and bars below
marking excess in group A. Intervals flagged by the Higher Criticism threshold
are highlighted, giving a visual reading of which windows of the follow-up
drive a detection.
}
\examples{
set.seed(1)
T_A <- rexp(300, 1 / 10)
T_B <- c(rexp(250, 1 / 10), rexp(50, 1 / 2))
plot_pvalue_profile(T_A, T_B, n_intervals_to_pool = 50)
}
\seealso{
\code{\link{suspected_deviations}}
}
