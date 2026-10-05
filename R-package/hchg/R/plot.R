#' Plot the per-interval p-value profile
#'
#' Draws the signed \eqn{-\log_{10}} per-interval hypergeometric p-value
#' profile, with bars above the axis marking excess events in group B and bars
#' below marking excess in group A. Intervals flagged by the Higher Criticism
#' threshold are highlighted, giving a visual reading of which windows of the
#' follow-up drive a detection.
#'
#' @inheritParams higher_criticism_test
#' @param col_flagged Colour for intervals at or below the HC threshold.
#' @param col_plain Colour for the remaining intervals.
#' @param xlab,ylab,main Standard graphical parameters.
#' @param ... Further arguments passed to \code{\link[graphics]{plot.default}}.
#'
#' @return Invisibly, the data frame returned by
#'   \code{\link{suspected_deviations}} with \code{alternative = "both"}.
#'
#' @seealso \code{\link{suspected_deviations}}
#' @examples
#' set.seed(1)
#' T_A <- rexp(300, 1 / 10)
#' T_B <- c(rexp(250, 1 / 10), rexp(50, 1 / 2))
#' plot_pvalue_profile(T_A, T_B, n_intervals_to_pool = 50)
#' @export
plot_pvalue_profile <- function(durations_A, durations_B,
                                event_observed_A = NULL,
                                event_observed_B = NULL,
                                gamma = 0.2, stbl = TRUE, t_0 = -1,
                                n_intervals_to_pool = NULL,
                                col_flagged = "firebrick",
                                col_plain = "steelblue",
                                xlab = "Time", ylab = "signed -log10(p)",
                                main = NULL, ...) {
  d <- suspected_deviations(durations_A, durations_B,
                            event_observed_A, event_observed_B,
                            alternative = "both", gamma = gamma, stbl = stbl,
                            t_0 = t_0, n_intervals_to_pool = n_intervals_to_pool)

  tt <- if ("time" %in% names(d)) d$time else (d$time_lower + d$time_upper) / 2

  ## Signed profile: positive where group B shows excess events, negative
  ## where group A does. The direction with the smaller p-value wins.
  up <- -log10(pmax(d$hypergeom_pvalue, .Machine$double.xmin))
  dn <- -log10(pmax(d$hypergeom_pvalue_rev, .Machine$double.xmin))
  signed <- ifelse(d$hypergeom_pvalue <= d$hypergeom_pvalue_rev, up, -dn)

  flagged <- d$suspected
  cols <- ifelse(flagged, col_flagged, col_plain)

  if (is.null(main)) main <- "Per-interval hypergeometric p-values"

  graphics::plot.default(tt, signed, type = "n", xlab = xlab, ylab = ylab,
                         main = main, ...)
  graphics::abline(h = 0, col = "grey40")
  graphics::segments(tt, 0, tt, signed, col = cols, lwd = 2)

  thr <- d$hc_threshold[1L]
  thr_rev <- d$hc_threshold_rev[1L]
  if (is.finite(thr) && thr > 0) {
    graphics::abline(h = -log10(thr), lty = 2, col = "grey30")
  }
  if (is.finite(thr_rev) && thr_rev > 0) {
    graphics::abline(h = log10(thr_rev), lty = 2, col = "grey30")
  }

  graphics::legend("topright", bty = "n",
                   legend = c("HC-flagged", "not flagged"),
                   col = c(col_flagged, col_plain), lwd = 2)

  invisible(d)
}
