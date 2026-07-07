"""
Higher Criticism and related tests for non-proportional hazard deviations.

These tests extend the ``lifelines`` survival analysis library with omnibus
tests that are powerful against *rare and weak* departures from the null
hypothesis of equal hazard functions.

Reference
---------
Kipnis, A., Galili, B., and Yakhini, Z. (2025). Higher criticism for
rare and weak non-proportional hazard deviations in survival analysis.
*Biometrika*, asaf075.
"""

import numpy as np
import pandas as pd
from lifelines.statistics import StatisticalResult
from lifelines.utils import group_survival_table_from_events
from multiHGtest import hypergeom_test
from multitest import MultiTest


__all__ = [
    "higher_criticism_test",
    "berk_jones_test",
    "fisher_combination_test",
    "min_p_test",
    "event_pvalues",
    "suspected_deviations",
]

_DEFAULT_GAMMA = 0.2


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _prepare_inputs(durations_A, durations_B, event_observed_A, event_observed_B, t_0):
    """Validate and coerce inputs, apply time restriction ``t_0``."""
    durations_A = np.asarray(durations_A, dtype=float)
    durations_B = np.asarray(durations_B, dtype=float)

    if event_observed_A is None:
        event_observed_A = np.ones(len(durations_A))
    else:
        event_observed_A = np.asarray(event_observed_A, dtype=float)

    if event_observed_B is None:
        event_observed_B = np.ones(len(durations_B))
    else:
        event_observed_B = np.asarray(event_observed_B, dtype=float)

    durations = np.concatenate([durations_A, durations_B])
    groups = np.concatenate([
        np.zeros(len(durations_A), dtype=int),
        np.ones(len(durations_B), dtype=int),
    ])
    events = np.concatenate([event_observed_A, event_observed_B])

    if t_0 >= 0:
        events = events.copy()
        events[durations > t_0] = 0
        durations = np.minimum(durations, t_0)

    return durations, groups, events


def _survival_table_counts(durations_A, durations_B,
                           event_observed_A, event_observed_B,
                           t_0=-1, n_intervals_to_pool=None):
    """Build the two-group survival table and return per-event-time counts.

    Parameters
    ----------
    n_intervals_to_pool : int, optional
        If given, aggregate event times into *n_intervals_to_pool* equal-width intervals.
        Recommended for continuous survival data where individual event times
        rarely have more than one event.

    Returns
    -------
    Nt1, Nt2 : ndarray
        At-risk counts in groups A and B at each (binned) event time.
    Ot1, Ot2 : ndarray
        Observed event counts in groups A and B at each (binned) event time.
    """
    durations, groups, events = _prepare_inputs(
        durations_A, durations_B, event_observed_A, event_observed_B, t_0
    )

    _, rm, obs, _ = group_survival_table_from_events(
        groups, durations, events
    )

    # At-risk: N_j(t) = total_ever_removed_j - cumulative_removed_before_t_j
    at_risk = rm.sum(0).values - rm.cumsum(0).shift(1).fillna(0)

    Nt1 = at_risk.iloc[:, 0].values.astype(int)
    Nt2 = at_risk.iloc[:, 1].values.astype(int)
    Ot1 = obs.iloc[:, 0].values.astype(int)
    Ot2 = obs.iloc[:, 1].values.astype(int)
    event_times = obs.index.values.astype(float)

    # Keep only times with events and positive at-risk in both groups
    valid = (Nt1 > 0) & (Nt2 > 0) & ((Ot1 + Ot2) > 0)
    Nt1, Nt2 = Nt1[valid], Nt2[valid]
    Ot1, Ot2 = Ot1[valid], Ot2[valid]
    event_times = event_times[valid]

    bin_width = None
    if n_intervals_to_pool is not None and len(Nt1) > n_intervals_to_pool:
        Nt1, Nt2, Ot1, Ot2, event_times, bin_width = _bin_counts(
            Nt1, Nt2, Ot1, Ot2, event_times, n_intervals_to_pool,
        )

    return Nt1, Nt2, Ot1, Ot2, event_times, bin_width


def _bin_counts(Nt1, Nt2, Ot1, Ot2, event_times, n_intervals_to_pool):
    """Aggregate per-event-time counts into *n_intervals_to_pool* equal-width intervals.

    Intervals span from 0 to the last event time.  Within each bin the
    observed events are summed and the at-risk count is taken from the
    first event time in the bin.  Empty bins (no events) are kept and
    their at-risk counts are forward-filled from the nearest earlier
    non-empty bin (or from the first event time for bins that precede
    all events).

    Returns the binned counts **and** the bin midpoints.
    """
    t_max = event_times[-1]
    edges = np.linspace(0, t_max * (1 + 1e-10), n_intervals_to_pool + 1)
    midpoints = 0.5 * (edges[:-1] + edges[1:])
    bin_idx = np.digitize(event_times, edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_intervals_to_pool - 1)

    b_Nt1 = np.zeros(n_intervals_to_pool, dtype=int)
    b_Nt2 = np.zeros(n_intervals_to_pool, dtype=int)
    b_Ot1 = np.zeros(n_intervals_to_pool, dtype=int)
    b_Ot2 = np.zeros(n_intervals_to_pool, dtype=int)

    for b in range(n_intervals_to_pool):
        mask = bin_idx == b
        if not mask.any():
            continue
        first = np.argmax(mask)
        b_Nt1[b] = Nt1[first]
        b_Nt2[b] = Nt2[first]
        b_Ot1[b] = Ot1[mask].sum()
        b_Ot2[b] = Ot2[mask].sum()

    # Forward-fill at-risk for empty bins.  Bins before the first event
    # get the at-risk count from the earliest event time (≈ total
    # subjects).  Gaps between events inherit the last known at-risk.
    last_n1, last_n2 = Nt1[0], Nt2[0]
    for b in range(n_intervals_to_pool):
        if b_Nt1[b] == 0 and b_Nt2[b] == 0:
            b_Nt1[b] = last_n1
            b_Nt2[b] = last_n2
        else:
            last_n1 = b_Nt1[b]
            last_n2 = b_Nt2[b]

    bin_width = edges[1] - edges[0]
    keep = (b_Nt1 > 0) & (b_Nt2 > 0)
    return (b_Nt1[keep], b_Nt2[keep], b_Ot1[keep], b_Ot2[keep],
            midpoints[keep], bin_width)


def _one_direction_pvals(Nt1, Nt2, Ot1, Ot2):
    """Per-event-time hypergeometric p-values testing excess events in group 2."""
    return hypergeom_test(
        Ot2, Nt2 + Nt1, Nt2, Ot1 + Ot2,
        randomize=False, alternative="greater",
    )


def _aggregate(pvals, method, gamma, stbl):
    """Aggregate a vector of p-values into a single test statistic."""
    if len(pvals) == 0:
        return 0.0
    mt = MultiTest(pvals, stbl=stbl)
    if method == "hc":
        return mt.hc(gamma=gamma)[0]
    if method == "berk_jones":
        return mt.berkjones(gamma=gamma)
    if method == "fisher":
        return mt.fisher()[0]
    if method == "min_p":
        return mt.minp()
    raise ValueError(f"Unknown aggregation method: {method!r}")


def _stat_from_counts(Nt1, Nt2, Ot1, Ot2, method, alternative, gamma, stbl):
    """Compute the aggregated test statistic from survival-table counts."""
    if len(Nt1) == 0:
        return 0.0

    if alternative == "both":
        stat_g = _aggregate(
            _one_direction_pvals(Nt1, Nt2, Ot1, Ot2), method, gamma, stbl
        )
        stat_l = _aggregate(
            _one_direction_pvals(Nt2, Nt1, Ot2, Ot1), method, gamma, stbl
        )
        return max(stat_g, stat_l)

    if alternative == "greater":
        pvals = _one_direction_pvals(Nt1, Nt2, Ot1, Ot2)
    else:
        pvals = _one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
    return _aggregate(pvals, method, gamma, stbl)


def _permutation_pvalue(durations_A, durations_B,
                        event_observed_A, event_observed_B,
                        observed_stat, n_permutations, method,
                        alternative, gamma, stbl, t_0, n_intervals_to_pool, rng):
    """Estimate a p-value by permuting group labels.

    Returns ``(p_value, perm_stats)`` where *perm_stats* is the array of
    test statistics computed under the permuted data.
    """
    durations_A = np.asarray(durations_A, dtype=float)
    durations_B = np.asarray(durations_B, dtype=float)
    event_observed_A = (
        np.asarray(event_observed_A, dtype=float)
        if event_observed_A is not None
        else np.ones(len(durations_A))
    )
    event_observed_B = (
        np.asarray(event_observed_B, dtype=float)
        if event_observed_B is not None
        else np.ones(len(durations_B))
    )

    all_dur = np.concatenate([durations_A, durations_B])
    all_evt = np.concatenate([event_observed_A, event_observed_B])
    n_A = len(durations_A)

    perm_stats = np.empty(n_permutations)
    for i in range(n_permutations):
        perm = rng.permutation(len(all_dur))
        d_A, d_B = all_dur[perm[:n_A]], all_dur[perm[n_A:]]
        e_A, e_B = all_evt[perm[:n_A]], all_evt[perm[n_A:]]
        Nt1, Nt2, Ot1, Ot2, *_ = _survival_table_counts(
            d_A, d_B, e_A, e_B, t_0, n_intervals_to_pool,
        )
        perm_stats[i] = _stat_from_counts(
            Nt1, Nt2, Ot1, Ot2, method, alternative, gamma, stbl,
        )

    p_value = (np.sum(perm_stats >= observed_stat) + 1) / (n_permutations + 1)
    return p_value, perm_stats


# ---------------------------------------------------------------------------
# Generic dispatcher (shared logic for all four public tests)
# ---------------------------------------------------------------------------

def _run_test(durations_A, durations_B, event_observed_A, event_observed_B,
              method, test_name, alternative, gamma, stbl, t_0, n_intervals_to_pool,
              n_permutations, seed, **extra_kw):
    """Shared implementation behind every public test function."""
    Nt1, Nt2, Ot1, Ot2, *_ = _survival_table_counts(
        durations_A, durations_B, event_observed_A, event_observed_B,
        t_0, n_intervals_to_pool,
    )

    stat = _stat_from_counts(Nt1, Nt2, Ot1, Ot2, method, alternative, gamma, stbl)

    perm_stats = None
    if n_permutations > 0:
        rng = np.random.default_rng(seed)
        p_value, perm_stats = _permutation_pvalue(
            durations_A, durations_B, event_observed_A, event_observed_B,
            stat, n_permutations, method, alternative, gamma, stbl,
            t_0, n_intervals_to_pool, rng,
        )
    else:
        p_value = np.nan

    result = StatisticalResult(
        p_value=p_value,
        test_statistic=stat,
        test_name=test_name,
        alternative=alternative,
        gamma=gamma,
        stbl=stbl,
        n_intervals_to_pool=n_intervals_to_pool,
        n_permutations=n_permutations,
        **extra_kw,
    )
    result.permutation_statistics = perm_stats
    return result


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def higher_criticism_test(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="both",
    gamma=_DEFAULT_GAMMA,
    stbl=True,
    t_0=-1,
    n_intervals_to_pool=None,
    n_permutations=0,
    seed=None,
    **kwargs,
) -> StatisticalResult:
    r"""Test for non-proportional hazard deviations using Higher Criticism.

    At every event time a hypergeometric p-value measures the evidence for
    unequal hazards; these p-values are then aggregated by the Higher
    Criticism (HC) statistic.  HC is especially powerful when hazard
    differences are *rare* (occur at few event times) and *weak* (small
    effect at each time).

    Parameters
    ----------
    durations_A : array-like
        Event/censoring times for group A.
    durations_B : array-like
        Event/censoring times for group B.
    event_observed_A : array-like, optional
        1 = event, 0 = censored for group A.  Default: all events.
    event_observed_B : array-like, optional
        1 = event, 0 = censored for group B.  Default: all events.
    alternative : ``{'both', 'greater', 'less'}``
        * ``'greater'`` — test for excess hazard in group B.
        * ``'less'`` — test for excess hazard in group A.
        * ``'both'`` (default) — maximum of both directions.
    gamma : float
        HC fraction parameter (default 0.2).  Only ordered p-values with
        rank :math:`\le \gamma n` are considered.
    stbl : bool
        Use the variance-stabilised HC denominator (default ``True``).
    t_0 : float
        Restrict to events before ``t_0`` (``-1`` = no restriction).
    n_intervals_to_pool : int, optional
        Number of equal-width time bins to aggregate event times into.
        **Recommended for continuous survival data** where individual event
        times rarely have more than one event.  When ``None`` (default) no
        binning is applied.
    n_permutations : int
        Number of label permutations for estimating the p-value.
        ``0`` (default) skips the permutation test and returns ``p_value=NaN``.
    seed : int, optional
        Random seed for the permutation test.

    Returns
    -------
    StatisticalResult
        ``.test_statistic`` is the HC score; ``.p_value`` is the
        permutation p-value (or ``NaN`` when ``n_permutations=0``).

    Notes
    -----
    The per-event hypergeometric test has very coarse resolution when each
    event time contains at most one event (typical of continuous survival
    data).  In that regime, set ``n_intervals_to_pool`` to a moderate value (e.g. 50–200)
    to pool events into time intervals and restore statistical power.

    References
    ----------
    Kipnis, A., Galili, B., and Yakhini, Z. (2025). Higher criticism
    for rare and weak non-proportional hazard deviations in survival
    analysis. *Biometrika*, asaf075.

    Examples
    --------
    >>> from lifelines_hc import higher_criticism_test
    >>> import numpy as np
    >>> rng = np.random.default_rng(0)
    >>> T_A = rng.exponential(10, size=200)
    >>> T_B = rng.exponential(10, size=200)
    >>> result = higher_criticism_test(T_A, T_B, n_intervals_to_pool=50)
    >>> result.print_summary()
    """
    return _run_test(
        durations_A, durations_B, event_observed_A, event_observed_B,
        method="hc",
        test_name="Higher Criticism test for non-proportional hazards",
        alternative=alternative, gamma=gamma, stbl=stbl, t_0=t_0,
        n_intervals_to_pool=n_intervals_to_pool, n_permutations=n_permutations, seed=seed, **kwargs,
    )


def berk_jones_test(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="both",
    gamma=_DEFAULT_GAMMA,
    stbl=True,
    t_0=-1,
    n_intervals_to_pool=None,
    n_permutations=0,
    seed=None,
    **kwargs,
) -> StatisticalResult:
    r"""Test for non-proportional hazard deviations using Berk-Jones.

    Same per-event hypergeometric p-values as
    :func:`higher_criticism_test`, aggregated via the Berk-Jones
    goodness-of-fit statistic.

    Parameters
    ----------
    durations_A, durations_B, event_observed_A, event_observed_B,
    alternative, gamma, stbl, t_0, n_intervals_to_pool, n_permutations, seed
        See :func:`higher_criticism_test`.

    Returns
    -------
    StatisticalResult
    """
    return _run_test(
        durations_A, durations_B, event_observed_A, event_observed_B,
        method="berk_jones",
        test_name="Berk-Jones test for non-proportional hazards",
        alternative=alternative, gamma=gamma, stbl=stbl, t_0=t_0,
        n_intervals_to_pool=n_intervals_to_pool, n_permutations=n_permutations, seed=seed, **kwargs,
    )


def fisher_combination_test(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="both",
    gamma=_DEFAULT_GAMMA,
    stbl=True,
    t_0=-1,
    n_intervals_to_pool=None,
    n_permutations=0,
    seed=None,
    **kwargs,
) -> StatisticalResult:
    r"""Test for non-proportional hazard deviations using Fisher's combination.

    Same per-event hypergeometric p-values as
    :func:`higher_criticism_test`, aggregated via Fisher's combination
    method (:math:`-2 \sum \log p_i`).

    Parameters
    ----------
    durations_A, durations_B, event_observed_A, event_observed_B,
    alternative, gamma, stbl, t_0, n_intervals_to_pool, n_permutations, seed
        See :func:`higher_criticism_test`.

    Returns
    -------
    StatisticalResult
    """
    return _run_test(
        durations_A, durations_B, event_observed_A, event_observed_B,
        method="fisher",
        test_name="Fisher combination test for non-proportional hazards",
        alternative=alternative, gamma=gamma, stbl=stbl, t_0=t_0,
        n_intervals_to_pool=n_intervals_to_pool, n_permutations=n_permutations, seed=seed, **kwargs,
    )


def min_p_test(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="both",
    gamma=_DEFAULT_GAMMA,
    stbl=True,
    t_0=-1,
    n_intervals_to_pool=None,
    n_permutations=0,
    seed=None,
    **kwargs,
) -> StatisticalResult:
    r"""Test for non-proportional hazard deviations using the minimum p-value.

    Same per-event hypergeometric p-values as
    :func:`higher_criticism_test`, aggregated by taking the minimum.

    Parameters
    ----------
    durations_A, durations_B, event_observed_A, event_observed_B,
    alternative, gamma, stbl, t_0, n_intervals_to_pool, n_permutations, seed
        See :func:`higher_criticism_test`.

    Returns
    -------
    StatisticalResult
    """
    return _run_test(
        durations_A, durations_B, event_observed_A, event_observed_B,
        method="min_p",
        test_name="Minimum-p test for non-proportional hazards",
        alternative=alternative, gamma=gamma, stbl=stbl, t_0=t_0,
        n_intervals_to_pool=n_intervals_to_pool, n_permutations=n_permutations, seed=seed, **kwargs,
    )


def event_pvalues(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="both",
    t_0=-1,
    n_intervals_to_pool=None,
):
    """Return per-event-time hypergeometric p-values.

    Useful for diagnostics and visualisation: the p-value at each event
    time quantifies the evidence against equal hazards at that time.

    Parameters
    ----------
    durations_A, durations_B, event_observed_A, event_observed_B, t_0, n_intervals_to_pool
        See :func:`higher_criticism_test`.
    alternative : ``{'both', 'greater', 'less'}``
        * ``'greater'`` — one-sided p-values for excess hazard in B.
        * ``'less'`` — one-sided p-values for excess hazard in A.
        * ``'both'`` — returns a pair ``(pvals_greater, pvals_less)``.

    Returns
    -------
    ndarray or tuple of ndarray
    """
    Nt1, Nt2, Ot1, Ot2, *_ = _survival_table_counts(
        durations_A, durations_B, event_observed_A, event_observed_B,
        t_0, n_intervals_to_pool,
    )
    if alternative == "both":
        return (
            _one_direction_pvals(Nt1, Nt2, Ot1, Ot2),
            _one_direction_pvals(Nt2, Nt1, Ot2, Ot1),
        )
    if alternative == "greater":
        return _one_direction_pvals(Nt1, Nt2, Ot1, Ot2)
    return _one_direction_pvals(Nt2, Nt1, Ot2, Ot1)


def suspected_deviations(
    durations_A,
    durations_B,
    event_observed_A=None,
    event_observed_B=None,
    alternative="greater",
    gamma=_DEFAULT_GAMMA,
    stbl=True,
    t_0=-1,
    n_intervals_to_pool=None,
):
    r"""Identify time intervals with suspected non-proportional hazard deviations.

    Computes per-interval hypergeometric p-values and flags those at or
    below the Higher Criticism threshold.  The returned DataFrame is
    suitable for visualization — gray-shading the flagged intervals on a
    Kaplan-Meier plot, for example.

    Parameters
    ----------
    durations_A, durations_B, event_observed_A, event_observed_B,
    alternative, gamma, stbl, t_0, n_intervals_to_pool
        See :func:`higher_criticism_test`.

    Returns
    -------
    pandas.DataFrame
        Indexed by event time (or bin midpoint when ``n_intervals_to_pool`` is set),
        with columns:

        * ``at_risk_A``, ``at_risk_B`` — subjects at risk in each group.
        * ``observed_A``, ``observed_B`` — observed events.
        * ``pvalue`` — one-sided hypergeometric p-value.
        * ``suspected`` — ``True`` where the p-value ≤ HC threshold.
        * ``hc_threshold`` — the HC p-value threshold (same for all rows).
    """
    Nt1, Nt2, Ot1, Ot2, times, bw = _survival_table_counts(
        durations_A, durations_B, event_observed_A, event_observed_B,
        t_0, n_intervals_to_pool,
    )

    if alternative == "less":
        pvals = _one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
    else:
        pvals = _one_direction_pvals(Nt1, Nt2, Ot1, Ot2)

    usable = pvals <= 1
    if usable.any():
        mt = MultiTest(pvals[usable], stbl=stbl)
        hc_score, hc_thresh = mt.hc(gamma=gamma)
    else:
        hc_score, hc_thresh = 0.0, 0.0

    flagged = pvals <= hc_thresh

    # Build index: interval ranges when binned, plain times otherwise
    if bw is not None:
        half = bw / 2.0
        left = times - half
        right = times + half
        prec = max(0, int(-np.floor(np.log10(bw))) + 2) if bw > 0 else 2
        idx = pd.Index(
            [f"{l:.{prec}f}-{r:.{prec}f}" for l, r in zip(left, right)],
            name="time_interval",
        )
    else:
        idx = pd.Index(np.round(times, 4), name="time")

    df = pd.DataFrame({
        "at_risk_A": Nt1,
        "at_risk_B": Nt2,
        "observed_A": Ot1,
        "observed_B": Ot2,
        "hypergeom_pvalue": np.round(pvals, 6),
        "suspected": flagged,
        "hc_threshold": hc_thresh,
    }, index=idx)

    if alternative == "both":
        pvals_rev = _one_direction_pvals(Nt2, Nt1, Ot2, Ot1)
        usable_rev = pvals_rev <= 1
        if usable_rev.any():
            mt_rev = MultiTest(pvals_rev[usable_rev], stbl=stbl)
            _, hc_thresh_rev = mt_rev.hc(gamma=gamma)
        else:
            hc_thresh_rev = 0.0
        flagged_rev = pvals_rev <= hc_thresh_rev
        df["hypergeom_pvalue_rev"] = np.round(pvals_rev, 6)
        df["suspected_rev"] = flagged_rev
        df["hc_threshold_rev"] = hc_thresh_rev
        df["suspected"] = flagged | flagged_rev

    return df
