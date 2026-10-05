% Generated from R/tests.R
\name{event_pvalues}
\alias{event_pvalues}
\title{Per-interval hypergeometric p-values}
\usage{
event_pvalues(
  durations_A,
  durations_B,
  event_observed_A = NULL,
  event_observed_B = NULL,
  alternative = c("both", "greater", "less"),
  t_0 = -1,
  n_intervals_to_pool = NULL
)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.}

\item{event_observed_A, event_observed_B}{Optional 0/1 vectors, 1 indicating an
observed event and 0 a censored observation. Default: all events.}

\item{alternative}{One of \code{"both"} (default), returning a list with
components \code{greater} and \code{less}; \code{"greater"}; or \code{"less"}.}

\item{t_0}{Restrict the analysis to events at or before \code{t_0}. A negative
value (default) applies no restriction.}

\item{n_intervals_to_pool}{Number of equal-width intervals to pool event times
into. \code{NULL} (default) applies no pooling.}
}
\value{
A numeric vector, or a list of two numeric vectors when
\code{alternative = "both"}.
}
\description{
Returns the exact one-sided hypergeometric p-value for each event time, or for
each pooled interval when \code{n_intervals_to_pool} is supplied. Each p-value
quantifies the evidence against equal hazards at that point in the follow-up.
}
\examples{
set.seed(1)
p <- event_pvalues(rexp(200, 1 / 10), rexp(200, 1 / 10),
                   alternative = "greater", n_intervals_to_pool = 50)
summary(p)
}
\seealso{
\code{\link{suspected_deviations}}
}
