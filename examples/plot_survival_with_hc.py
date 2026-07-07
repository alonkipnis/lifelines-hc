"""
Example: Kaplan-Meier curves with Higher Criticism flagged intervals.

Generates synthetic two-group survival data where group B has elevated
hazard at a sparse subset of times, fits Kaplan-Meier curves, and
highlights the time intervals that HC identifies as having a suspected
non-proportional hazard deviation.

Usage:

    python examples/plot_survival_with_hc.py

Reference:
    Kipnis, A., Galili, B., and Yakhini, Z. (2025). Higher criticism
    for rare and weak non-proportional hazard deviations in survival
    analysis. *Biometrika*, asaf075.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

from lifelines_hc import suspected_deviations, KaplanMeierHCIllustrator

mpl.style.use("ggplot")


# ---------------------------------------------------------------------------
# 1. Generate synthetic survival data
# ---------------------------------------------------------------------------

def generate_data(n_A=500, n_B=500, frac_affected=0.12, censor_time=30.0,
                  seed=42):
    """Create two-group survival data with a rare hazard deviation in B.

    Most subjects in both groups follow Exp(scale=10).  A fraction
    ``frac_affected`` of group-B subjects instead follow Exp(scale=2),
    modelling a sparse subpopulation with elevated hazard.  Administrative
    censoring is applied at ``censor_time``.
    """
    rng = np.random.default_rng(seed)

    T_A = rng.exponential(10, size=n_A)
    T_B = np.where(
        rng.random(n_B) < frac_affected,
        rng.exponential(2, size=n_B),
        rng.exponential(10, size=n_B),
    )

    event_A = (T_A <= censor_time).astype(float)
    event_B = (T_B <= censor_time).astype(float)
    T_A = np.minimum(T_A, censor_time)
    T_B = np.minimum(T_B, censor_time)

    return T_A, T_B, event_A, event_B


# ---------------------------------------------------------------------------
# 2. Main
# ---------------------------------------------------------------------------

def main():
    T_A, T_B, event_A, event_B = generate_data()

    n_intervals_to_pool = 100
    gamma = 0.2

    ill = KaplanMeierHCIllustrator(T_A, T_B, event_A, event_B)

    # --- Compute suspected deviations ---
    df_dev = suspected_deviations(
        T_A, T_B, event_A, event_B,
        alternative="greater", gamma=gamma,
        n_intervals_to_pool=n_intervals_to_pool,
    )

    # --- HC test with permutation p-value ---
    hc_result = ill.test(
        n_intervals_to_pool=n_intervals_to_pool, gamma=gamma,
        alternative="greater", n_permutations=500, seed=0,
    )
    title = ill.format_title(hc_result)

    # --- Plot ---
    fig, ax = plt.subplots(figsize=(10, 5))
    ill.plot(df_dev, ax=ax, title=title)
    fig.tight_layout()
    fig.savefig("km_with_hc.png", dpi=180, bbox_inches="tight",
                pad_inches=0.05)
    print("Saved figure to km_with_hc.png")

    # --- Print table ---
    print()
    ill.print_table(df_dev)

    plt.show()


if __name__ == "__main__":
    main()
