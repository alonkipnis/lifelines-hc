"""Tests for lifelines_hc.statistics."""

import numpy as np
import pytest
from lifelines.statistics import StatisticalResult

from lifelines_hc import (
    higher_criticism_test,
    berk_jones_test,
    fisher_combination_test,
    min_p_test,
    event_pvalues,
)

N_BINS = 50


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture
def null_data():
    """Two groups drawn from the same exponential — null is true."""
    rng = np.random.default_rng(42)
    return {
        "durations_A": rng.exponential(10, size=500),
        "durations_B": rng.exponential(10, size=500),
    }


@pytest.fixture
def alt_data():
    """Group B has higher hazard at a sparse subset of times (strong effect)."""
    rng = np.random.default_rng(42)
    T_A = rng.exponential(10, size=500)
    T_B = np.where(
        rng.random(500) < 0.25,
        rng.exponential(2, size=500),
        rng.exponential(10, size=500),
    )
    return {"durations_A": T_A, "durations_B": T_B}


@pytest.fixture
def censored_data():
    """Censored survival data under the null."""
    rng = np.random.default_rng(7)
    n = 200
    durations_A = rng.exponential(10, size=n)
    durations_B = rng.exponential(10, size=n)
    censor_time = 12.0
    event_A = (durations_A <= censor_time).astype(float)
    event_B = (durations_B <= censor_time).astype(float)
    durations_A = np.minimum(durations_A, censor_time)
    durations_B = np.minimum(durations_B, censor_time)
    return {
        "durations_A": durations_A,
        "durations_B": durations_B,
        "event_observed_A": event_A,
        "event_observed_B": event_B,
    }


# ---------------------------------------------------------------------------
# Basic return type and shape
# ---------------------------------------------------------------------------

class TestReturnType:
    def test_returns_statistical_result(self, null_data):
        result = higher_criticism_test(**null_data, n_intervals_to_pool=N_BINS)
        assert isinstance(result, StatisticalResult)

    def test_has_test_statistic_and_pvalue(self, null_data):
        result = higher_criticism_test(**null_data, n_intervals_to_pool=N_BINS)
        assert hasattr(result, "test_statistic")
        assert hasattr(result, "p_value")
        assert np.isfinite(result.test_statistic)

    def test_pvalue_nan_without_permutations(self, null_data):
        result = higher_criticism_test(**null_data, n_intervals_to_pool=N_BINS)
        assert np.isnan(result.p_value)

    def test_pvalue_finite_with_permutations(self, null_data):
        result = higher_criticism_test(
            **null_data, n_intervals_to_pool=N_BINS, n_permutations=50, seed=0,
        )
        assert np.isfinite(result.p_value)
        assert 0 < result.p_value <= 1

    def test_works_without_binning(self, null_data):
        result = higher_criticism_test(**null_data)
        assert isinstance(result, StatisticalResult)


# ---------------------------------------------------------------------------
# All four tests run without error
# ---------------------------------------------------------------------------

class TestAllTests:
    @pytest.mark.parametrize("test_func", [
        higher_criticism_test,
        berk_jones_test,
        fisher_combination_test,
        min_p_test,
    ])
    def test_runs_on_null(self, null_data, test_func):
        result = test_func(**null_data, n_intervals_to_pool=N_BINS)
        assert isinstance(result, StatisticalResult)
        assert not np.isnan(result.test_statistic)

    @pytest.mark.parametrize("test_func", [
        higher_criticism_test,
        berk_jones_test,
        fisher_combination_test,
        min_p_test,
    ])
    def test_runs_on_alt(self, alt_data, test_func):
        result = test_func(**alt_data, n_intervals_to_pool=N_BINS)
        assert isinstance(result, StatisticalResult)
        assert not np.isnan(result.test_statistic)


# ---------------------------------------------------------------------------
# Alternative directions
# ---------------------------------------------------------------------------

class TestAlternative:
    def test_greater_vs_less(self, alt_data):
        res_g = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS, alternative="greater")
        res_l = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS, alternative="less")
        assert res_g.test_statistic >= res_l.test_statistic

    def test_both_geq_either(self, alt_data):
        res_both = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS, alternative="both")
        res_g = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS, alternative="greater")
        res_l = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS, alternative="less")
        assert res_both.test_statistic >= res_g.test_statistic - 1e-12
        assert res_both.test_statistic >= res_l.test_statistic - 1e-12


# ---------------------------------------------------------------------------
# Censored data
# ---------------------------------------------------------------------------

class TestCensoring:
    def test_handles_censoring(self, censored_data):
        result = higher_criticism_test(**censored_data, n_intervals_to_pool=N_BINS)
        assert isinstance(result, StatisticalResult)
        assert np.isfinite(result.test_statistic)


# ---------------------------------------------------------------------------
# Power: HC should be larger under the alternative than under the null
# ---------------------------------------------------------------------------

class TestPower:
    def test_hc_higher_under_alt(self, null_data, alt_data):
        hc_null = higher_criticism_test(**null_data, n_intervals_to_pool=N_BINS).test_statistic
        hc_alt = higher_criticism_test(**alt_data, n_intervals_to_pool=N_BINS).test_statistic
        assert hc_alt > hc_null


# ---------------------------------------------------------------------------
# Permutation p-value
# ---------------------------------------------------------------------------

class TestPermutation:
    def test_pvalue_under_null_not_tiny(self, null_data):
        result = higher_criticism_test(
            **null_data, n_intervals_to_pool=N_BINS, n_permutations=199, seed=0,
        )
        assert result.p_value > 0.01

    def test_pvalue_under_alt_is_small(self, alt_data):
        result = higher_criticism_test(
            **alt_data, n_intervals_to_pool=N_BINS, n_permutations=199, seed=0,
        )
        assert result.p_value < 0.15

    def test_seed_reproducibility(self, null_data):
        kw = dict(**null_data, n_intervals_to_pool=N_BINS, n_permutations=99, seed=123)
        r1 = higher_criticism_test(**kw)
        r2 = higher_criticism_test(**kw)
        assert r1.p_value == r2.p_value
        assert r1.test_statistic == r2.test_statistic


# ---------------------------------------------------------------------------
# event_pvalues
# ---------------------------------------------------------------------------

class TestEventPvalues:
    def test_returns_array(self, null_data):
        pv = event_pvalues(**null_data, alternative="greater", n_intervals_to_pool=N_BINS)
        assert isinstance(pv, np.ndarray)
        assert len(pv) > 0
        assert np.all((pv >= 0) & (pv <= 1))

    def test_both_returns_tuple(self, null_data):
        result = event_pvalues(**null_data, alternative="both", n_intervals_to_pool=N_BINS)
        assert isinstance(result, tuple)
        assert len(result) == 2
        assert len(result[0]) == len(result[1])

    def test_n_intervals_to_pool_controls_length(self, null_data):
        pv = event_pvalues(**null_data, alternative="greater", n_intervals_to_pool=N_BINS)
        assert len(pv) <= N_BINS


# ---------------------------------------------------------------------------
# t_0 restriction
# ---------------------------------------------------------------------------

class TestTimeRestriction:
    def test_t0_reduces_events(self, null_data):
        pv_all = event_pvalues(**null_data, alternative="greater")
        pv_restricted = event_pvalues(
            **null_data, alternative="greater", t_0=5.0,
        )
        assert len(pv_restricted) < len(pv_all)
