% Generated from R/tests.R
\name{event_pvalues}
\alias{event_pvalues}
\alias{event_pvalues.default}
\alias{event_pvalues.formula}
\title{Per-interval hypergeometric p-values}
\usage{
event_pvalues(durations_A, ...)

\method{event_pvalues}{default}(
  durations_A,
  durations_B,
  event_observed_A = NULL,
  event_observed_B = NULL,
  alternative = c("both", "greater", "less"),
  t_0 = -1,
  n_intervals_to_pool = NULL,
  ...
)

\method{event_pvalues}{formula}(formula, data, subset, na.action, ...)
}
\arguments{
\item{durations_A, durations_B}{Numeric vectors of event or censoring times.
\code{durations_A} may instead be a \code{\link{formula}}, which selects the
formula method.}

\item{event_observed_A, event_observed_B}{Optional 0/1 vectors, 1 indicating an
observed event and 0 a censored observation. Default: all events.}

\item{alternative}{One of \code{"both"} (default), returning a list with
components \code{greater} and \code{less}; \code{"greater"}; or
\code{"less"}.}

\item{t_0}{Restrict the analysis to events at or before \code{t_0}. A negative
value (default) applies no restriction.}

\item{n_intervals_to_pool}{Number of equal-width intervals to pool event times
into. Strongly recommended for continuous survival data, where each distinct
time carries at most one event and the per-time hypergeometric p-values are too
coarse to be informative. Values between 50 and 200 usually work well.
\code{NULL} (default) applies no pooling.}

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
A numeric vector, or a list of two numeric vectors when
\code{alternative = "both"}.
}
\description{
Returns the exact one-sided hypergeometric p-value for each event time, or for
each pooled interval when \code{n_intervals_to_pool} is supplied. Each p-value
quantifies the evidence against equal hazards at that point in the follow-up.
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
p <- event_pvalues(rexp(200, 1 / 10), rexp(200, 1 / 10),
                   alternative = "greater", n_intervals_to_pool = 50)
summary(p)

trial <- data.frame(time = rexp(400, 1 / 10), status = 1,
                    arm = rep(c("A", "B"), each = 200))
summary(event_pvalues(cbind(time, status) ~ arm, data = trial,
                      alternative = "greater", n_intervals_to_pool = 50))
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
\code{\link{suspected_deviations}}
}
