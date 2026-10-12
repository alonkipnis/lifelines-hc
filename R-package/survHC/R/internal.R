## Internal helpers: survival-table construction, interval pooling, and the
## per-interval exact hypergeometric p-values.
##
## These mirror the reference Python implementation (package 'lifelines-hc')
## so that the test statistics agree to numerical precision.

.DEFAULT_GAMMA <- 0.2

## Validate and coerce inputs; apply the time restriction t_0.
.prepare_inputs <- function(durations_A, durations_B,
                            event_observed_A, event_observed_B, t_0) {
  durations_A <- as.numeric(durations_A)
  durations_B <- as.numeric(durations_B)

  if (!length(durations_A) || !length(durations_B)) {
    stop("both 'durations_A' and 'durations_B' must be non-empty", call. = FALSE)
  }
  if (anyNA(durations_A) || anyNA(durations_B)) {
    stop("'durations_A' and 'durations_B' must not contain NA", call. = FALSE)
  }
  if (any(durations_A < 0) || any(durations_B < 0)) {
    stop("durations must be non-negative", call. = FALSE)
  }

  event_observed_A <- if (is.null(event_observed_A)) {
    rep(1, length(durations_A))
  } else {
    as.numeric(event_observed_A)
  }
  event_observed_B <- if (is.null(event_observed_B)) {
    rep(1, length(durations_B))
  } else {
    as.numeric(event_observed_B)
  }

  if (length(event_observed_A) != length(durations_A)) {
    stop("'event_observed_A' must have the same length as 'durations_A'",
         call. = FALSE)
  }
  if (length(event_observed_B) != length(durations_B)) {
    stop("'event_observed_B' must have the same length as 'durations_B'",
         call. = FALSE)
  }

  durations <- c(durations_A, durations_B)
  groups <- c(rep(0L, length(durations_A)), rep(1L, length(durations_B)))
  events <- c(event_observed_A, event_observed_B)

  if (t_0 >= 0) {
    events[durations > t_0] <- 0
    durations <- pmin(durations, t_0)
  }

  list(durations = durations, groups = groups, events = events)
}

## Two-group survival table reduced to per-event-time counts.
##
## Returns at-risk counts (Nt1, Nt2), observed events (Ot1, Ot2), the event
## times, and the pooling bin width (NULL when no pooling was applied).
.survival_table_counts <- function(durations_A, durations_B,
                                   event_observed_A = NULL,
                                   event_observed_B = NULL,
                                   t_0 = -1, n_intervals_to_pool = NULL) {
  p <- .prepare_inputs(durations_A, durations_B,
                       event_observed_A, event_observed_B, t_0)
  durations <- p$durations
  groups <- p$groups
  events <- p$events

  times <- sort(unique(durations))
  nt <- length(times)
  idx <- match(durations, times)

  inA <- groups == 0L
  inB <- groups == 1L
  isev <- events > 0

  ## "Removed" counts everyone leaving the risk set (events and censorings);
  ## "observed" counts events only.
  rem1 <- tabulate(idx[inA], nbins = nt)
  rem2 <- tabulate(idx[inB], nbins = nt)
  obs1 <- tabulate(idx[inA & isev], nbins = nt)
  obs2 <- tabulate(idx[inB & isev], nbins = nt)

  ## At risk at t = group size minus everyone removed strictly before t.
  cum1 <- cumsum(rem1)
  cum2 <- cumsum(rem2)
  Nt1 <- sum(rem1) - c(0, cum1[-nt])
  Nt2 <- sum(rem2) - c(0, cum2[-nt])

  keep <- Nt1 > 0 & Nt2 > 0 & (obs1 + obs2) > 0
  Nt1 <- Nt1[keep]; Nt2 <- Nt2[keep]
  Ot1 <- obs1[keep]; Ot2 <- obs2[keep]
  event_times <- times[keep]

  bin_width <- NULL
  if (!is.null(n_intervals_to_pool) && length(Nt1) > n_intervals_to_pool) {
    binned <- .bin_counts(Nt1, Nt2, Ot1, Ot2, event_times,
                          as.integer(n_intervals_to_pool))
    Nt1 <- binned$Nt1; Nt2 <- binned$Nt2
    Ot1 <- binned$Ot1; Ot2 <- binned$Ot2
    event_times <- binned$times
    bin_width <- binned$bin_width
  }

  list(Nt1 = Nt1, Nt2 = Nt2, Ot1 = Ot1, Ot2 = Ot2,
       times = event_times, bin_width = bin_width)
}

## Aggregate per-event-time counts into equal-width intervals.
##
## Within a bin, events are summed and the at-risk count is taken from the
## earliest event time in that bin. Bins containing no events inherit the
## at-risk counts of the nearest earlier non-empty bin.
.bin_counts <- function(Nt1, Nt2, Ot1, Ot2, event_times, n_bins) {
  t_max <- event_times[length(event_times)]
  edges <- seq(0, t_max * (1 + 1e-10), length.out = n_bins + 1L)
  midpoints <- 0.5 * (edges[-(n_bins + 1L)] + edges[-1L])

  b <- findInterval(event_times, edges)
  b <- pmin(pmax(b, 1L), n_bins)

  b_Nt1 <- integer(n_bins); b_Nt2 <- integer(n_bins)
  b_Ot1 <- integer(n_bins); b_Ot2 <- integer(n_bins)

  for (k in seq_len(n_bins)) {
    m <- which(b == k)
    if (!length(m)) next
    first <- m[1L]
    b_Nt1[k] <- Nt1[first]
    b_Nt2[k] <- Nt2[first]
    b_Ot1[k] <- sum(Ot1[m])
    b_Ot2[k] <- sum(Ot2[m])
  }

  last1 <- Nt1[1L]; last2 <- Nt2[1L]
  for (k in seq_len(n_bins)) {
    if (b_Nt1[k] == 0L && b_Nt2[k] == 0L) {
      b_Nt1[k] <- last1
      b_Nt2[k] <- last2
    } else {
      last1 <- b_Nt1[k]
      last2 <- b_Nt2[k]
    }
  }

  keep <- b_Nt1 > 0L & b_Nt2 > 0L
  list(Nt1 = b_Nt1[keep], Nt2 = b_Nt2[keep],
       Ot1 = b_Ot1[keep], Ot2 = b_Ot2[keep],
       times = midpoints[keep],
       bin_width = edges[2L] - edges[1L])
}

## Exact one-sided hypergeometric p-values for excess events in group B.
##
## Under the null of equal group-specific hazards, the number of group-B
## events among the events observed in an interval is hypergeometric with
## population Nt1 + Nt2, Nt2 "successes", and Ot1 + Ot2 draws.
.one_direction_pvals <- function(Nt1, Nt2, Ot1, Ot2) {
  stats::phyper(Ot2 - 1, m = Nt2, n = Nt1, k = Ot1 + Ot2, lower.tail = FALSE)
}
