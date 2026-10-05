"""
Domain 3: Adjuvant Bisphosphonate Therapy — delayed menopause-dependent effect.

Analyses the AZURE trial (Coleman et al. 2011, N Engl J Med; updated 2014) comparing
standard adjuvant therapy alone (control) versus the same therapy with added
zoledronic acid (a bisphosphonate) in early-stage breast cancer patients.

Biological context:
  Zoledronic acid targets the bone microenvironment — the primary niche for dormant
  micrometastases — by inhibiting osteoclast-mediated bone resorption.  This effect
  is strongly modulated by estrogen status: in postmenopausal women (low estrogen,
  high baseline bone turnover) the drug substantially reduces recurrence; in
  premenopausal women (high estrogen, low baseline bone turnover) the effect is
  minimal.  Because the trial enrolled both groups, and premenopausal patients
  progressively transition to postmenopause during follow-up, the net hazard
  difference is *temporally concentrated* rather than uniformly elevated — producing
  the characteristic non-proportional-hazards pattern.

Statistical consequence:
  Log-rank (and all four weighted NPH alternatives) average over the entire
  follow-up and are unable to separate the null early period from the beneficial
  late period, yielding p ≈ 0.30.  HC detects the localized intervals of hazard
  divergence: p ≈ 0.005.

Data source: AppNote/data/AZURE_2A.csv  (run 01_get_kmdata.py first)
  Columns: time (months), event (0/1), arm (control / zoledronic_acid)
  n = 3359 (1678 control, 1681 zoledronic acid); endpoint: disease-free survival.

Reference:
  Coleman RE, et al. (2011). Breast-cancer adjuvant therapy with zoledronic acid.
  N Engl J Med, 365(15), 1396-1405.  DOI: 10.1056/NEJMoa1105195
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use("ggplot")

SCRIPT_DIR  = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from utils import (
    run_all_tests, pvalue_profile,
    plot_km_with_hc, plot_pvalue_profile, print_results_table,
)

DATA_DIR    = SCRIPT_DIR / "data"
FIGS_DIR    = SCRIPT_DIR / "figs"
RES_DIR     = SCRIPT_DIR / "results"
N_INTERVALS = 80   # ~1.5-month bins over 120-month follow-up
N_PERMS     = 2000

for d in (FIGS_DIR, RES_DIR):
    d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_azure(data_dir: Path) -> tuple:
    """Load AZURE DFS IPD from CSV; return (T_ctrl, T_trt, E_ctrl, E_trt)."""
    csv_path = data_dir / "AZURE_2A.csv"
    if not csv_path.exists():
        print("[INFO] AZURE_2A.csv not found — using synthetic fallback.")
        print("       Run '01_get_kmdata.py --trials AZURE_2A' to download real data.")
        return _make_azure_synthetic()

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()
    print(f"Loaded AZURE_2A.csv: {len(df)} rows, columns: {list(df.columns)}")

    ctrl = df[df["arm"] == "control"]
    trt  = df[df["arm"] == "zoledronic_acid"]

    T_ctrl = ctrl["time"].values.astype(float)
    T_trt  = trt["time"].values.astype(float)
    E_ctrl = ctrl["event"].values.astype(float)
    E_trt  = trt["event"].values.astype(float)

    print(f"  Control:         n={len(T_ctrl)}, events={int(E_ctrl.sum())}, "
          f"median={np.median(T_ctrl):.1f} mo")
    print(f"  Zoledronic acid: n={len(T_trt)},  events={int(E_trt.sum())}, "
          f"median={np.median(T_trt):.1f} mo")
    return T_ctrl, T_trt, E_ctrl, E_trt


def _make_azure_synthetic(n: int = 1200, seed: int = 7) -> tuple:
    """Synthetic fallback: delayed bisphosphonate benefit model.

    Simulates the AZURE menopause-dependent pattern:
    - Control arm:  proportional-hazards baseline throughout.
    - Treatment arm: identical hazard for first 24 months (null premenopausal
      period), then 25 % reduced hazard (HR = 0.75) from month 24 onward as
      postmenopausal patients accumulate.
    """
    rng  = np.random.default_rng(seed)
    lam0 = np.log(2) / 80.0   # baseline hazard, median ~80 months
    t0   = 24.0               # onset of bisphosphonate benefit
    hr   = 0.75               # hazard ratio in the benefit window
    censor = 120.0

    T_ctrl = rng.exponential(1 / lam0, n)
    E_ctrl = (T_ctrl <= censor).astype(float)
    T_ctrl = np.minimum(T_ctrl, censor)

    S_t0 = np.exp(-lam0 * t0)
    u = rng.uniform(0, 1, n)
    die_early = u > S_t0
    T_late = t0 - np.log(np.clip(u / S_t0, 1e-300, 1)) / (lam0 * hr)
    T_trt = np.where(die_early,
                     -np.log(np.clip(u, 1e-300, 1)) / lam0,
                     T_late)
    E_trt = (T_trt <= censor).astype(float)
    T_trt = np.minimum(T_trt, censor)

    print("[SYNTHETIC] Delayed-benefit model (illustrative; use real data for paper).")
    return T_ctrl, T_trt, E_ctrl, E_trt


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def run_analysis(T_ctrl, T_trt, E_ctrl, E_trt):
    print(f"\n{'='*60}")
    print("ADJUVANT BISPHOSPHONATE — AZURE trial (DFS)")
    print(f"{'='*60}")

    print(f"\nRunning all tests (n_intervals={N_INTERVALS}, "
          f"n_permutations={N_PERMS}) ...")
    results = run_all_tests(
        T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        n_permutations=N_PERMS,
        label_A="Control",
        label_B="Zoledronic acid",
    )
    print_results_table(results)

    out_csv = RES_DIR / "azure_test_results.csv"
    results.to_csv(out_csv)
    print(f"Results saved to {out_csv}")
    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figure(T_ctrl, T_trt, E_ctrl, E_trt):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Adjuvant bisphosphonate therapy: AZURE trial DFS "
        "(zoledronic acid vs control)",
        fontsize=13, fontweight="bold",
    )

    ax_km = axes[0]
    df_dev = plot_km_with_hc(
        ax_km, T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        label_A="Control",
        label_B="Zoledronic acid",
        shade_color="darkorange",
        title="A  Kaplan-Meier curves (DFS)",
        xlabel="Time (months)",
    )
    n_flagged = df_dev["suspected"].sum()
    ax_km.text(0.98, 0.98,
               f"HC-flagged intervals: {n_flagged}/{len(df_dev)}",
               transform=ax_km.transAxes, ha="right", va="top", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    ax_pv = axes[1]
    plot_pvalue_profile(
        ax_pv, df_dev,
        title=r"B  Per-interval hypergeometric $p$-values",
        xlabel="Time interval (months)",
    )

    fig.tight_layout()
    out_path = FIGS_DIR / "azure_km.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    T_ctrl, T_trt, E_ctrl, E_trt = load_azure(DATA_DIR)
    run_analysis(T_ctrl, T_trt, E_ctrl, E_trt)
    make_figure(T_ctrl, T_trt, E_ctrl, E_trt)


if __name__ == "__main__":
    main()
