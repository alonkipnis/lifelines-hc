## Formula interface.
##
## Every exported entry point is an S3 generic with a "default" method taking
## the two samples as vectors (mirroring the Python package 'lifelines-hc') and
## a "formula" method taking the familiar
##
##     Surv(time, status) ~ group
##
## The formula method only rearranges its inputs: it splits the data into two
## samples and hands them to the default method, so the two interfaces return
## the same numbers by construction. The only behavioural difference is the
## treatment of missing values, which 'na.action' controls here and which the
## default method rejects outright.
##
## 'survival' is in Suggests, not Imports: a Surv object is recognised by its
## structure, which is all this package needs. A user who can build one already
## has 'survival' loaded, and 'cbind(time, status) ~ group' works without it.

## Build the model frame from the matched call, deferring 'data', 'subset' and
## 'na.action' to stats::model.frame exactly as stats::t.test.formula does.
.formula_model_frame <- function(m, env) {
  keep <- match(c("formula", "data", "subset", "na.action"), names(m), 0L)
  m <- m[c(1L, keep)]
  m$drop.unused.levels <- TRUE
  m[[1L]] <- quote(stats::model.frame)
  eval(m, env)
}

## Time and status from the left-hand side of the formula.
##
## Accepts a Surv object of type "right", a two-column matrix of time and
## status (so the interface works without 'survival' attached), or a bare
## numeric vector of times, which is read as "every observation is an event".
.as_time_status <- function(y) {
  if (is.null(y)) {
    stop("'formula' must have a left-hand side", call. = FALSE)
  }

  if (inherits(y, "Surv")) {
    type <- attr(y, "type")
    if (!identical(type, "right")) {
      stop("survHC handles right-censored data only, but this Surv() object ",
           "has type \"", type, "\"", call. = FALSE)
    }
    cn <- colnames(y)
    it <- if (!is.null(cn) && "time" %in% cn) "time" else 1L
    is_ <- if (!is.null(cn) && "status" %in% cn) "status" else 2L
    return(list(time = as.numeric(y[, it]), status = as.numeric(y[, is_])))
  }

  if (is.matrix(y)) {
    if (ncol(y) != 2L) {
      stop("a matrix left-hand side must have two columns, time and status, ",
           "but this one has ", ncol(y), call. = FALSE)
    }
    return(list(time = as.numeric(y[, 1L]), status = as.numeric(y[, 2L])))
  }

  if (is.numeric(y)) {
    return(list(time = as.numeric(y), status = rep(1, length(y))))
  }

  stop("the left-hand side of 'formula' must be a Surv() object, a two-column ",
       "matrix of time and status, or a numeric vector of times", call. = FALSE)
}

## Index the two groups named by the right-hand side.
##
## The first level is group A and the second is group B, so that
## alternative = "greater" tests for excess hazard in the second level.
## Unordered variables are levelled by sort(unique(.)): numerically for
## numbers, alphabetically for strings, FALSE before TRUE for logicals.
.group_index <- function(g, label) {
  if (is.matrix(g) || is.data.frame(g)) {
    stop("the grouping variable '", label, "' must be a vector", call. = FALSE)
  }
  g <- if (is.factor(g)) droplevels(g) else factor(g)
  lev <- levels(g)

  if (length(lev) != 2L) {
    stop("the grouping variable '", label, "' must take exactly two values, ",
         "but it takes ", length(lev),
         if (length(lev)) paste0(" (", paste(lev, collapse = ", "), ")") else "",
         call. = FALSE)
  }

  list(A = which(g == lev[1L]), B = which(g == lev[2L]), levels = lev)
}

## Split the data named by a formula into the two samples the default methods
## expect. Returns the four vectors plus a data.name for the htest and the
## group levels, in the order (A, B).
.formula_split <- function(formula, m, env, dots) {
  if (missing(formula) || !inherits(formula, "formula") ||
        length(formula) != 3L) {
    stop("'formula' must be a two-sided formula such as ",
         "Surv(time, status) ~ group", call. = FALSE)
  }

  clash <- intersect(names(dots),
                     c("durations_A", "durations_B",
                       "event_observed_A", "event_observed_B"))
  if (length(clash)) {
    stop("'", clash[1L], "' cannot be combined with a formula; name the ",
         "times and the status indicator on the left-hand side, as in ",
         "Surv(time, status) ~ group", call. = FALSE)
  }

  mf <- .formula_model_frame(m, env)

  ## One response column and one grouping column. A formula reusing the
  ## response as the grouping variable collapses to a single column, hence the
  ## separate message.
  if (ncol(mf) < 2L) {
    stop("the right-hand side of 'formula' must name one grouping variable, ",
         "distinct from the response", call. = FALSE)
  }
  if (ncol(mf) > 2L) {
    stop("the right-hand side of 'formula' must name exactly one grouping ",
         "variable, but ", ncol(mf) - 1L, " were supplied", call. = FALSE)
  }

  ts <- .as_time_status(stats::model.response(mf))
  gi <- .group_index(mf[[2L]], names(mf)[2L])

  list(durations_A = ts$time[gi$A],
       durations_B = ts$time[gi$B],
       event_observed_A = ts$status[gi$A],
       event_observed_B = ts$status[gi$B],
       group_levels = gi$levels,
       data_name = paste0(deparse(formula[[2L]], nlines = 1L), " by ",
                          names(mf)[2L], " (", gi$levels[1L], ", ",
                          gi$levels[2L], ")"))
}

## Hand a formula split to a default method.
##
## Named arguments are used throughout so that the dots collected by the
## formula method land on the same formals they would have matched in a direct
## call, which is what makes the two interfaces agree exactly.
.formula_plain <- function(fun, s, dots) {
  do.call(fun, c(list(durations_A = s$durations_A,
                      durations_B = s$durations_B,
                      event_observed_A = s$event_observed_A,
                      event_observed_B = s$event_observed_B), dots))
}

## As .formula_plain, but relabels the htest with the formula and the group
## levels, which the default method cannot know.
.formula_htest <- function(fun, s, dots) {
  res <- .formula_plain(fun, s, dots)
  res$data.name <- s$data_name
  res$group_levels <- s$group_levels
  res
}

## The default methods carry '...' because their generic does, which would
## otherwise let a misspelled argument through unnoticed. Silently dropping
## 'n_intervals_to_pool' would change the statistic, so reject the leftovers.
.reject_unused <- function(...) {
  n <- ...length()
  if (n) {
    nm <- names(list(...))
    nm <- if (is.null(nm)) character() else nm[nzchar(nm)]
    stop(if (length(nm)) {
      paste0("unused argument(s): ", paste(nm, collapse = ", "))
    } else {
      paste0(n, " unused argument(s)")
    }, call. = FALSE)
  }
  invisible(NULL)
}

## data.name for the default methods.
##
## Under do.call() substitute() yields the vector itself rather than a symbol;
## the fallbacks keep a whole data vector out of the deparse, and nlines caps
## a long inline expression at one line so data.name stays a single string.
.data_name <- function(a, b) {
  one <- function(e, fallback) {
    if (is.name(e) || is.call(e)) {
      paste(deparse(e, nlines = 1L), collapse = " ")
    } else {
      fallback
    }
  }
  paste(one(a, "durations_A"), "and", one(b, "durations_B"))
}
