# lifelines-hc

Higher Criticism and related tests for detecting **non-proportional hazard
deviations** in two-sample survival data — implemented as a
[lifelines](https://lifelines.readthedocs.io/) extension.

These tests are especially powerful when hazard differences are *rare* (occur
at few event times) and *weak* (small effect at each time), a regime where the
classical log-rank test has little power.

## Reference

> Kipnis, A., Galili, B., and Yakhini, Z. (2025). Higher criticism for rare
> and weak non-proportional hazard deviations in survival analysis.
> *Biometrika*, asaf075.

## Installation

```bash
pip install lifelines-hc
```

Or install from a local checkout in development mode:

```bash
cd lifelines-hc
pip install -e ".[dev]"
```

## Quick start

```python
import numpy as np
from lifelines_hc import higher_criticism_test

rng = np.random.default_rng(42)
durations_A = rng.exponential(10, size=300)
durations_B = rng.exponential(10, size=300)

result = higher_criticism_test(durations_A, durations_B, n_intervals_to_pool=50)
result.print_summary()
```

The result is a `lifelines.statistics.StatisticalResult`, so it plugs straight
into any lifelines workflow.

> **Note on `n_intervals_to_pool`:** The per-event hypergeometric test has very coarse
> resolution when each event time has at most one event (typical of continuous
> survival data). Setting `n_intervals_to_pool` pools events into equal-width time intervals,
> restoring statistical power. A value between 50 and 200 usually works well.
> For data with naturally discrete or tied event times this parameter can be
> omitted.

### With a permutation p-value

```python
result = higher_criticism_test(
    durations_A, durations_B,
    n_intervals_to_pool=50,
    n_permutations=1000,
    seed=0,
)
print(f"HC = {result.test_statistic:.3f}, p = {result.p_value:.4f}")
```

### Handling censored data

```python
event_A = rng.binomial(1, 0.8, size=300)
event_B = rng.binomial(1, 0.8, size=300)

result = higher_criticism_test(
    durations_A, durations_B,
    event_observed_A=event_A,
    event_observed_B=event_B,
    n_intervals_to_pool=50,
)
result.print_summary()
```

### Other tests

The same per-event hypergeometric p-values can be aggregated with different
statistics:

```python
from lifelines_hc import berk_jones_test, fisher_combination_test, min_p_test

bj = berk_jones_test(durations_A, durations_B, n_intervals_to_pool=50)
fc = fisher_combination_test(durations_A, durations_B, n_intervals_to_pool=50)
mp = min_p_test(durations_A, durations_B, n_intervals_to_pool=50)
```

### Inspecting per-event p-values

```python
from lifelines_hc import event_pvalues

pvals = event_pvalues(durations_A, durations_B, alternative="greater", n_intervals_to_pool=50)
```

### Identifying and visualising suspected intervals

`suspected_deviations` returns a DataFrame flagging the time bins where
the p-value falls at or below the HC threshold — ready for gray-shading on a
Kaplan-Meier plot:

```python
from lifelines_hc import suspected_deviations

df = suspected_deviations(durations_A, durations_B, n_intervals_to_pool=50)
print(df[df["suspected"]])
```

See `examples/plot_survival_with_hc.py` for a full working example that
produces a figure like:

![KM with HC](examples/km_with_hc.png)

## API

| Function                    | Returns                     |
|-----------------------------|-----------------------------|
| `higher_criticism_test`     | `StatisticalResult` (HC)    |
| `berk_jones_test`           | `StatisticalResult` (BJ)    |
| `fisher_combination_test`   | `StatisticalResult` (Fisher)|
| `min_p_test`                | `StatisticalResult` (min-p) |
| `event_pvalues`             | ndarray of p-values         |
| `suspected_deviations`      | DataFrame with flags        |

All test functions accept the same parameters and return a
`lifelines.statistics.StatisticalResult`.

## License

MIT
