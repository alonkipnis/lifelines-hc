"""Kaplan-Meier visualization with Higher Criticism flagged intervals.

Provides :class:`KaplanMeierHCIllustrator`, a lightweight wrapper around
``lifelines.KaplanMeierFitter`` that overlays the suspected non-proportional
hazard intervals identified by :func:`~lifelines_hc.suspected_deviations`.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter

from lifelines_hc.statistics import higher_criticism_test, suspected_deviations


class KaplanMeierHCIllustrator:
    r"""Plot Kaplan-Meier survival curves and shade HC-suspected intervals.

    Parameters
    ----------
    durations_A, durations_B : array-like
        Event / censoring times for the two groups.
    event_observed_A, event_observed_B : array-like, optional
        1 = event, 0 = censored.  Default: all events.
    label_A, label_B : str
        Legend labels for the two curves.

    Examples
    --------
    >>> from lifelines_hc import suspected_deviations
    >>> from lifelines_hc.plotting import KaplanMeierHCIllustrator
    >>> df_dev = suspected_deviations(T_A, T_B, event_A, event_B,
    ...                               n_intervals_to_pool=100)
    >>> ill = KaplanMeierHCIllustrator(T_A, T_B, event_A, event_B)
    >>> ill.plot(df_dev)
    """

    def __init__(
        self,
        durations_A,
        durations_B,
        event_observed_A=None,
        event_observed_B=None,
        label_A=r"$\hat{S}_A$",
        label_B=r"$\hat{S}_B$",
    ):
        self.durations_A = np.asarray(durations_A, dtype=float)
        self.durations_B = np.asarray(durations_B, dtype=float)
        self.event_observed_A = (
            np.asarray(event_observed_A, dtype=float)
            if event_observed_A is not None else None
        )
        self.event_observed_B = (
            np.asarray(event_observed_B, dtype=float)
            if event_observed_B is not None else None
        )

        self.kmf_A = KaplanMeierFitter()
        self.kmf_A.fit(durations_A, event_observed_A, label=label_A)

        self.kmf_B = KaplanMeierFitter()
        self.kmf_B.fit(durations_B, event_observed_B, label=label_B)

    # ------------------------------------------------------------------
    # Plotting
    # ------------------------------------------------------------------

    def plot(
        self,
        deviations,
        ax=None,
        ci_show=False,
        show_censors=True,
        title=None,
        shade_color="k",
        shade_alpha=0.15,
        xlabel="Time Interval",
        ylabel="Survival proportion",
    ):
        """Plot Kaplan-Meier curves and shade suspected intervals.

        Parameters
        ----------
        deviations : pandas.DataFrame
            Output of :func:`~lifelines_hc.suspected_deviations`.
        ax : matplotlib.axes.Axes, optional
            Target axes.  Created automatically if *None*.
        ci_show : bool
            Show confidence intervals on the KM curves.
        show_censors : bool
            Show censor tick marks.
        title : str, optional
            Axes title.  Use :meth:`format_title` for a convenient
            HC / log-rank summary string.
        shade_color : str
            Colour of the shaded rectangles.
        shade_alpha : float
            Alpha (transparency) of the shaded rectangles.
        xlabel, ylabel : str
            Axis labels.

        Returns
        -------
        matplotlib.axes.Axes
        """
        if ax is None:
            _, ax = plt.subplots(figsize=(10, 5))

        self.kmf_A.plot_survival_function(
            ax=ax, ci_show=ci_show, show_censors=show_censors,
        )
        self.kmf_B.plot_survival_function(
            ax=ax, ci_show=ci_show, show_censors=show_censors,
        )

        flagged = deviations[deviations["suspected"]]
        first_span = True
        for label in flagged.index:
            left, right = self._parse_interval(label)
            kw = dict(color=shade_color, alpha=shade_alpha, zorder=0)
            if first_span:
                kw["label"] = "increased hazard suspicion"
                first_span = False
            ax.axvspan(left, right, **kw)

        if title is not None:
            ax.set_title(title, fontsize=14)
        ax.set_xlabel(xlabel, fontsize=14)
        ax.set_ylabel(ylabel, fontsize=14)
        ax.legend(fontsize=14, loc="lower left")

        return ax

    # ------------------------------------------------------------------
    # HC test
    # ------------------------------------------------------------------

    def test(self, n_intervals_to_pool=None, gamma=0.2, alternative="greater",
             stbl=True, t_0=-1, n_permutations=0, seed=None, **kwargs):
        """Run the Higher Criticism test, optionally with a permutation p-value.

        Convenience wrapper around :func:`~lifelines_hc.higher_criticism_test`
        using the survival data stored in this illustrator.

        Parameters
        ----------
        n_intervals_to_pool, gamma, alternative, stbl, t_0, n_permutations, seed
            See :func:`~lifelines_hc.higher_criticism_test`.

        Returns
        -------
        StatisticalResult
        """
        return higher_criticism_test(
            self.durations_A, self.durations_B,
            self.event_observed_A, self.event_observed_B,
            alternative=alternative, gamma=gamma, stbl=stbl, t_0=t_0,
            n_intervals_to_pool=n_intervals_to_pool,
            n_permutations=n_permutations, seed=seed, **kwargs,
        )

    # ------------------------------------------------------------------
    # Permutation histogram
    # ------------------------------------------------------------------

    @staticmethod
    def plot_permutation_histogram(hc_result, ax=None, bins=30):
        """Plot a histogram of permutation HC statistics with the observed value.

        Parameters
        ----------
        hc_result : StatisticalResult
            Must have been produced with ``n_permutations > 0`` so that
            ``hc_result.permutation_statistics`` is available.
        ax : matplotlib.axes.Axes, optional
            Target axes.  Created automatically if *None*.
        bins : int
            Number of histogram bins.

        Returns
        -------
        matplotlib.axes.Axes or None
            *None* when no permutation statistics are available.
        """
        perm_stats = getattr(hc_result, "permutation_statistics", None)
        if perm_stats is None:
            return None

        if ax is None:
            _, ax = plt.subplots(figsize=(5, 3))

        ax.hist(perm_stats, bins=bins, color="#7fafdf", edgecolor="white",
                density=True, label="Permutations")
        ax.axvline(hc_result.test_statistic, color="crimson", linewidth=2,
                   linestyle="--", label=f"Observed = {hc_result.test_statistic:.2f}")
        if np.isfinite(hc_result.p_value):
            ax.set_title(f"Permutation p-value = {hc_result.p_value:.4f}",
                         fontsize=9)
        ax.set_xlabel("HC statistic", fontsize=8)
        ax.set_ylabel("Density", fontsize=8)
        ax.tick_params(labelsize=7)
        ax.legend(fontsize=7)

        return ax

    # ------------------------------------------------------------------
    # Title helper
    # ------------------------------------------------------------------

    @staticmethod
    def format_title(hc_result):
        r"""Build a LaTeX-formatted title from an HC test result.

        Displays the HC statistic and, when a permutation p-value is
        available (i.e. ``n_permutations > 0`` was used), the p-value.

        Parameters
        ----------
        hc_result : StatisticalResult
            From :func:`~lifelines_hc.higher_criticism_test` or
            :meth:`test`.

        Returns
        -------
        str
        """
        title = rf"$\mathrm{{HC}} = {hc_result.test_statistic:.2f}$"
        if np.isfinite(hc_result.p_value):
            title += rf",  $p = {hc_result.p_value:.4f}$"
        return title

    # ------------------------------------------------------------------
    # Table printing
    # ------------------------------------------------------------------

    @staticmethod
    def print_table(deviations):
        """Print a formatted table of all intervals, marking suspected with ``*``.

        Parameters
        ----------
        deviations : pandas.DataFrame
            Output of :func:`~lifelines_hc.suspected_deviations`.

        Returns
        -------
        str
            The formatted table (also printed to stdout).
        """
        display_cols = [
            c for c in ("at_risk_A", "at_risk_B",
                        "observed_A", "observed_B",
                        "hypergeom_pvalue")
            if c in deviations.columns
        ]
        table = deviations[display_cols]

        idx_width = max(
            max((len(str(i)) for i in table.index), default=0),
            len("time_interval"),
        )
        col_widths = {c: max(len(c), 8) for c in display_cols}

        lines = []

        hdr = f"  {'time_interval':<{idx_width}}"
        for c in display_cols:
            hdr += f"  {c:>{col_widths[c]}}"
        lines.append(hdr)
        lines.append("  " + "-" * (len(hdr) - 2))

        for idx, row in table.iterrows():
            prefix = "* " if deviations.at[idx, "suspected"] else "  "
            line = f"{str(idx):<{idx_width}}"
            for c in display_cols:
                v = row[c]
                if "pvalue" in c:
                    line += f"  {v:>{col_widths[c]}.4f}"
                else:
                    line += f"  {int(v):>{col_widths[c]}}"
            lines.append(f"{prefix}{line}")

        text = "\n".join(lines)
        print(text)
        return text

    # ------------------------------------------------------------------
    # Internal
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_interval(label):
        """Extract (left, right) floats from an interval index label.

        Handles both string ranges (``"7.12-7.42"``) produced when
        ``n_intervals_to_pool`` is set and plain numeric index values.
        """
        s = str(label)
        # String range: "left-right" — split on the hyphen that separates
        # the two numbers.  Because the numbers themselves may be negative
        # we split from the right on a hyphen preceded by a digit.
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
