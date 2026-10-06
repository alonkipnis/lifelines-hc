"""Regression test: the four Application Note case studies must reproduce the
HCHG values reported in the manuscript.

These pin the Donoho-Jin 2008 standardization (``hc_version='dj2008'``, the
package default). They fail if a different HC standardization is used by
mistake.
"""
import pathlib

import pandas as pd
import pytest

from lifelines_hc import higher_criticism_test

DATA = pathlib.Path(__file__).resolve().parent.parent / "AppNote" / "data"

# name: (csv, n_intervals, control arm, treatment arm, published HC, published p)
CASES = {
    "checkmate057": ("Checkmate057_1C.csv", 60, "d1", "nivolumab",
                     1.69490342, 0.0009995002498750),
    "comet1": ("COMET1_2A.csv", 80, "prednisone", "cabozantinib",
               1.17937570, 0.0139930034982508),
    "azure": ("AZURE_2A.csv", 80, "control", "zoledronic_acid",
              1.51301477, 0.0054972513743128),
    "csl": ("CSL.csv", 50, "placebo", "prednisone",
            1.25720345, 0.0019990004997501),
}


@pytest.mark.parametrize("name", list(CASES))
def test_published_hchg_values(name):
    csv, n_intervals, arm_a, arm_b, exp_stat, exp_p = CASES[name]
    path = DATA / csv
    if not path.exists():
        pytest.skip(f"{csv} not present (reconstructed IPD not checked out)")

    df = pd.read_csv(path)
    a = df[df.arm == arm_a]
    b = df[df.arm == arm_b]
    assert len(a) and len(b), f"arms not found in {csv}"

    res = higher_criticism_test(
        a.time.values, b.time.values,
        event_observed_A=a.event.values,
        event_observed_B=b.event.values,
        alternative="both",
        gamma=0.2,
        n_intervals_to_pool=n_intervals,
        n_permutations=2000,
        seed=42,
    )
    assert res.test_statistic == pytest.approx(exp_stat, abs=1e-7)
    assert res.p_value == pytest.approx(exp_p, abs=1e-12)
