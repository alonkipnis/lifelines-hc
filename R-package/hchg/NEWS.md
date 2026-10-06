# hchg 0.1.0

* First release: an R port of the Python package `lifelines-hc`.
* `higher_criticism_test()`, `berk_jones_test()`, `fisher_combination_test()`
  and `min_p_test()` for two-sample right-censored survival data.
* `event_pvalues()` and `suspected_deviations()` expose the per-interval
  hypergeometric p-values and the intervals flagged by the Higher Criticism
  threshold.
* `plot_pvalue_profile()` draws the signed per-interval p-value profile.
* `normalization` selects the Higher Criticism normalisation:
  `"donoho-jin2008"` (default, u_i = i/N) or `"beta"` (u_i = i/(N+1)), the
  latter reproducing the default of the Python `multitest` package.
* No dependencies beyond base R, `stats` and `graphics`.
