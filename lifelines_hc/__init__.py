"""
lifelines-hc: Higher Criticism tests for non-proportional hazard deviations.

A ``lifelines`` extension implementing the methods from:

    Kipnis, A., Galili, B., and Yakhini, Z. (2025). Higher criticism
    for rare and weak non-proportional hazard deviations in survival
    analysis. *Biometrika*, asaf075.
"""

from lifelines_hc.statistics import (
    higher_criticism_test,
    berk_jones_test,
    fisher_combination_test,
    min_p_test,
    event_pvalues,
    suspected_deviations,
)
from lifelines_hc.plotting import KaplanMeierHCIllustrator

__all__ = [
    "higher_criticism_test",
    "berk_jones_test",
    "fisher_combination_test",
    "min_p_test",
    "event_pvalues",
    "suspected_deviations",
    "KaplanMeierHCIllustrator",
]
