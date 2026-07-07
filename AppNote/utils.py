"""
Shared utilities for the lifelines-hc Application Note analyses.

Provides:
- run_all_tests()       : run HC, BJ, Fisher, MinP, and Log-rank
- pvalue_profile()      : per-interval hypergeometric p-values + HC threshold
- plot_km_with_hc()     : KM curves with HC-flagged intervals shaded
- plot_pvalue_profile() : bar chart of -log10(p) with HC threshold line
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

from lifelines_hc import (
    higher_criticism_test,
    fisher_combination_test,
    min_p_test,
    suspected_deviations,
    event_pvalues,
)

DEFAULT_N_INTERVALS = 100
DEFAULT_N_PERMS = 500
DEFAULT_SEED = 42
DEFAULT_GAMMA = 0.2


# ---------------------------------------------------------------------------
# Omnibus test runner
# ---------------------------------------------------------------------------

def run_all_tests(
    T_A, T_B,
    E_A=None, E_B=None,
    n_intervals: int = DEFAULT_N_INTERVALS,
    n_permutations: int = DEFAULT_N_PERMS,
    seed: int = DEFAULT_SEED,
    gamma: float = DEFAULT_GAMMA,
    label_A: str = "A",
    label_B: str = "B",
) -> pd.DataFrame:
    """Run all competing survival tests and return a comparison DataFrame.

    Methods compared
    ----------------
    Log-rank (standard)
        Global averaging statistic; powerful for proportional hazards.
    Gehan-Wilcoxon
        Weights by number at risk; emphasises early events.
    Tarone-Ware
        Weights by sqrt(number at risk); intermediate early emphasis.
    Peto-Prentice
        Product-limit weight; emphasises early events like Wilcoxon but
        more robust to tied times.
    Fleming-Harrington (1,1)
        Weights by S(t)(1-S(t)); emphasises middle of the follow-up,
        useful for detecting late effects after an initial null period.
    Higher Criticism (HC)
        Detects rare and weak deviations; optimal for sparse hazard hot-spots.
    Fisher combination
        Combines interval p-values via the chi-squared log-sum.
    Minimum p-value (MinP)
        Bonferroni-corrected single-interval test.

    Parameters
    ----------
    T_A, T_B : array-like
        Event / censoring times for groups A and B.
    E_A, E_B : array-like, optional
        Event indicators (1=event, 0=censored). Default: all events.
    n_intervals : int
        Number of equal-width time bins for the HC-based tests.
    n_permutations : int
        Label-permutation repetitions for calibrating HC-based p-values.
    seed : int
        Random seed for the permutation procedure.
    gamma : float
        HC fraction parameter (top-gamma fraction of ordered p-values).
    label_A, label_B : str
        Group names used in the output DataFrame.

    Returns
    -------
    pd.DataFrame
        One row per method with columns: method, statistic, p_value.
    """
    T_A, T_B = np.asarray(T_A, float), np.asarray(T_B, float)

    # Log-rank (standard)
    lr = logrank_test(T_A, T_B, E_A, E_B)
    rows = [
        dict(method="Log-rank",
             statistic=lr.test_statistic,
             p_value=lr.p_value),
    ]

    # Weighted log-rank variants for non-proportional hazards
    for method_name, weighting, kw in [
        ("Gehan-Wilcoxon",          "wilcoxon",          {}),
        ("Tarone-Ware",             "tarone-ware",        {}),
        ("Peto-Prentice",           "peto",               {}),
        ("Fleming-Harrington (1,1)","fleming-harrington", {"p": 1, "q": 1}),
    ]:
        r = logrank_test(T_A, T_B, E_A, E_B, weightings=weighting, **kw)
        rows.append(dict(method=method_name,
                         statistic=r.test_statistic,
                         p_value=r.p_value))

    shared_kw = dict(
        event_observed_A=E_A, event_observed_B=E_B,
        alternative="both", gamma=gamma,
        n_intervals_to_pool=n_intervals,
        n_permutations=n_permutations, seed=seed,
    )

    # HC
    hc = higher_criticism_test(T_A, T_B, **shared_kw)
    rows.append(dict(method="Higher Criticism (HC)",
                     statistic=hc.test_statistic,
                     p_value=hc.p_value))

    # Fisher
    fi = fisher_combination_test(T_A, T_B, **shared_kw)
    rows.append(dict(method="Fisher combination",
                     statistic=fi.test_statistic,
                     p_value=fi.p_value))

    # MinP
    mp = min_p_test(T_A, T_B, **shared_kw)
    rows.append(dict(method="MinP (Bonferroni)",
                     statistic=mp.test_statistic,
                     p_value=mp.p_value))

    df = pd.DataFrame(rows).set_index("method")
    df.attrs["label_A"] = label_A
    df.attrs["label_B"] = label_B
    return df


# ---------------------------------------------------------------------------
# Per-interval p-value profile
# ---------------------------------------------------------------------------

def pvalue_profile(
    T_A, T_B,
    E_A=None, E_B=None,
    n_intervals: int = DEFAULT_N_INTERVALS,
    gamma: float = DEFAULT_GAMMA,
    alternative: str = "both",
) -> pd.DataFrame:
    """Return per-interval hypergeometric p-values and the HC threshold.

    Useful for the inset diagnostic plots that show *where* in time the
    hazard deviation is occurring.

    Returns
    -------
    pd.DataFrame
        Columns: ``time``, ``pvalue``, ``suspected``, ``hc_threshold``.
    """
    df = suspected_deviations(
        T_A, T_B, E_A, E_B,
        alternative=alternative,
        gamma=gamma,
        n_intervals_to_pool=n_intervals,
    )
    return df


# ---------------------------------------------------------------------------
# Plotting helpers
# ---------------------------------------------------------------------------

def plot_km_with_hc(
    ax,
    T_A, T_B,
    E_A=None, E_B=None,
    n_intervals: int = DEFAULT_N_INTERVALS,
    gamma: float = DEFAULT_GAMMA,
    label_A: str = "Control",
    label_B: str = "Treatment",
    shade_color: str = "steelblue",
    shade_alpha: float = 0.20,
    ci_show: bool = False,
    show_censors: bool = True,
    title: str = None,
    xlabel: str = "Time",
    ylabel: str = "Survival probability",
) -> pd.DataFrame:
    """Plot KM curves and shade HC-suspected intervals on *ax*.

    Returns the suspected_deviations DataFrame for further use.
    """
    kmf_A = KaplanMeierFitter().fit(T_A, E_A, label=label_A)
    kmf_B = KaplanMeierFitter().fit(T_B, E_B, label=label_B)

    kmf_A.plot_survival_function(ax=ax, ci_show=ci_show,
                                  show_censors=show_censors)
    kmf_B.plot_survival_function(ax=ax, ci_show=ci_show,
                                  show_censors=show_censors)

    df_dev = pvalue_profile(T_A, T_B, E_A, E_B,
                             n_intervals=n_intervals, gamma=gamma)
    flagged = df_dev[df_dev["suspected"]]

    first = True
    for label in flagged.index:
        left, right = _parse_interval_label(label)
        kw = dict(color=shade_color, alpha=shade_alpha, zorder=0)
        if first:
            kw["label"] = r"$\Delta^*$ (HC-flagged)"
            first = False
        ax.axvspan(left, right, **kw)

    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    if title:
        ax.set_title(title, fontsize=12)
    ax.legend(fontsize=10, loc="lower left")

    return df_dev


def plot_pvalue_profile(
    ax,
    df_dev: pd.DataFrame,
    direction: str = "greater",          # kept for backward compatibility; unused
    color_suspect: str = "crimson",
    color_null: str = "steelblue",
    xlabel: str = "Time interval",
    ylabel: str = r"signed $-\log_{10}(p)$",
    title: str = "Per-interval hypergeometric p-values",
) -> None:
    """Signed per-interval p-value bar chart, highlighting suspected windows.

    For the two-sided test, each interval is summarised by whichever
    direction is more extreme:

    * **upward** bar  = excess events in group B (treatment), height
      :math:`-\\log_{10}(p_{\\text{fwd}})`;
    * **downward** bar = excess events in group A (control), depth
      :math:`-\\log_{10}(p_{\\text{rev}})` drawn below zero.

    A dashed HC-threshold line is drawn on each side, so a bar crosses a
    line **iff** the interval is flagged as suspected — making this panel
    consistent with the shaded intervals on the companion KM plot.

    Parameters
    ----------
    df_dev : pd.DataFrame
        Output of :func:`pvalue_profile` (with ``alternative='both'`` so the
        reverse-direction columns are present).
    """
    p_fwd = df_dev["hypergeom_pvalue"].values.astype(float).clip(1e-10, 1)
    has_rev = "hypergeom_pvalue_rev" in df_dev.columns
    if has_rev:
        p_rev = df_dev["hypergeom_pvalue_rev"].values.astype(float).clip(1e-10, 1)
        # per interval pick the more extreme direction; sign encodes which arm
        use_fwd = p_fwd <= p_rev
        score = np.where(use_fwd, -np.log10(p_fwd), np.log10(p_rev))
    else:
        score = -np.log10(p_fwd)

    suspected = df_dev.get("suspected", pd.Series(False, index=df_dev.index)).values
    colors = np.where(suspected, color_suspect, color_null)
    x = np.arange(len(df_dev))

    ax.bar(x, score, color=colors, width=0.8, alpha=0.85)
    ax.axhline(0, color="0.5", linewidth=0.8)

    thr = df_dev["hc_threshold"].iloc[0] if "hc_threshold" in df_dev.columns else None
    if thr is not None and 0 < thr < 1:
        ax.axhline(-np.log10(thr), color="k", linestyle="--", linewidth=1.2,
                   label=f"HC threshold ({thr:.3f})")
    if has_rev and "hc_threshold_rev" in df_dev.columns:
        thr_rev = df_dev["hc_threshold_rev"].iloc[0]
        if 0 < thr_rev < 1:
            ax.axhline(np.log10(thr_rev), color="k", linestyle="--", linewidth=1.2)

    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=9, loc="upper right")

    # Replace numeric x-ticks with time labels (every ~10 intervals)
    step = max(1, len(df_dev) // 10)
    tick_pos = x[::step]
    ax.set_xticks(tick_pos)
    ax.set_xticklabels(
        [str(df_dev.index[i]) for i in tick_pos],
        rotation=45, ha="right", fontsize=7,
    )


def print_results_table(df: pd.DataFrame) -> None:
    """Pretty-print the comparison DataFrame from :func:`run_all_tests`."""
    print(f"\n{'Method':<30} {'Statistic':>12} {'p-value':>10}")
    print("-" * 54)
    for method, row in df.iterrows():
        sig = "*" if row["p_value"] < 0.05 else " "
        print(f"{method:<30} {row['statistic']:>12.4f} {row['p_value']:>10.4f} {sig}")
    print()


# ---------------------------------------------------------------------------
# Internal
# ---------------------------------------------------------------------------

def _parse_interval_label(label: str) -> tuple[float, float]:
    """Parse an interval index label into (left, right) floats."""
    s = str(label)
    parts = s.split("-")
    if len(parts) >= 2:
        try:
            right = float(parts[-1])
            left = float("-".join(parts[:-1]))
            return left, right
        except ValueError:
            pass
    t = float(s)
    eps = max(abs(t) * 1e-3, 1e-3)
    return t - eps, t + eps
