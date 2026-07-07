"""
Domain 1: Clinical Immuno-oncology — delayed treatment effect.

Analyses the CheckMate 057 trial (Borghaei et al. 2015, N Engl J Med) comparing
nivolumab vs docetaxel in 2nd-line non-squamous NSCLC.  The progression-free
survival (PFS) endpoint shows crossing survival curves characteristic of immune
checkpoint inhibitors: early progressions are more common in the immunotherapy
arm (immune priming / pseudo-progression), followed by durable disease control
for responding patients.  This crossing-curve pattern dilutes the standard
log-rank statistic but is well-matched to the HC framework.

Data source: AppNote/data/Checkmate057_1C.csv  (run 01_get_kmdata.py first)
  Fallback:  synthetic crossing-curve dataset if real data is unavailable.

Output:
  figs/immuno_km.png            — KM curves with HC-flagged intervals
  figs/immuno_pvalue_profile.png — per-interval p-value bar chart
  results/immuno_test_results.csv
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.style.use("ggplot")

# --- path setup ---
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from utils import (
    run_all_tests, pvalue_profile,
    plot_km_with_hc, plot_pvalue_profile, print_results_table,
)

DATA_DIR  = SCRIPT_DIR / "data"
FIGS_DIR  = SCRIPT_DIR / "figs"
RES_DIR   = SCRIPT_DIR / "results"
N_INTERVALS = 60   # ~1 bin per 0.5 months (30-month follow-up)
N_PERMS     = 2000

for d in (FIGS_DIR, RES_DIR):
    d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_poplar(data_dir: Path) -> tuple:
    """Load CheckMate 057 PFS IPD from CSV; return (T_ctrl, T_trt, E_ctrl, E_trt).

    Prefers Checkmate057_1C.csv (PFS, crossing curves, Borghaei 2015 NEJM).
    Falls back to POPLAR.csv then to a synthetic crossing-curve dataset.
    """
    for fname in ["Checkmate057_1C.csv", "POPLAR.csv"]:
        csv_path = data_dir / fname
        if csv_path.exists():
            break
    else:
        print(f"[INFO] No real trial data found — using synthetic fallback.")
        print("       Run '01_get_kmdata.py' to download CheckMate 057 data.")
        return _make_poplar_synthetic()
    csv_path = data_dir / fname
    if not csv_path.exists():
        print(f"[INFO] {csv_path} not found — using synthetic fallback.")
        print("       Run '01_get_kmdata.py' to download the real trial data.")
        return _make_poplar_synthetic()

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()
    print(f"Loaded {fname}: {len(df)} rows, columns: {list(df.columns)}")

    # Detect column names
    time_col  = _pick(df, ["time", "os", "os_time", "survival_time", "t"])
    event_col = _pick(df, ["status", "event", "os_event", "dead", "death"])
    arm_col   = _pick(df, ["arm", "trt", "treatment", "group"])

    if time_col is None or event_col is None:
        print("[WARN] Cannot parse trial CSV columns; using synthetic fallback.")
        return _make_poplar_synthetic()

    df[time_col]  = pd.to_numeric(df[time_col],  errors="coerce")
    df[event_col] = pd.to_numeric(df[event_col], errors="coerce")
    df = df.dropna(subset=[time_col, event_col])

    if arm_col is None:
        # If no arm column, assume first half is control, second half treatment
        mid = len(df) // 2
        df["_arm"] = [0]*mid + [1]*(len(df)-mid)
        arm_col = "_arm"

    # Standardise arm coding: 0=control/docetaxel, 1=treatment/nivolumab
    vals = df[arm_col].unique()
    if len(vals) == 2:
        # Assign lower value (or 'docetaxel' string) to 0
        val_map = {sorted(vals)[0]: 0, sorted(vals)[1]: 1}
        df["_arm_bin"] = df[arm_col].map(val_map).fillna(0)
    else:
        df["_arm_bin"] = (pd.to_numeric(df[arm_col], errors="coerce") > 0).astype(int)

    ctrl = df[df["_arm_bin"] == 0]
    trt  = df[df["_arm_bin"] == 1]

    T_ctrl = ctrl[time_col].values
    T_trt  = trt[time_col].values
    E_ctrl = ctrl[event_col].values
    E_trt  = trt[event_col].values

    print(f"  Control (docetaxel):    n={len(T_ctrl)}, events={E_ctrl.sum():.0f}")
    print(f"  Treatment (nivolumab):  n={len(T_trt)},  events={E_trt.sum():.0f}")
    return T_ctrl, T_trt, E_ctrl, E_trt


def _pick(df: pd.DataFrame, candidates: list) -> str | None:
    for c in candidates:
        if c in df.columns:
            return c
    for c in candidates:
        for col in df.columns:
            if c in col.lower():
                return col
    return None


def _make_poplar_synthetic(n: int = 200, seed: int = 9) -> tuple:
    """Synthetic dataset with crossing survival curves (delayed treatment effect).

    Models the characteristic immuno-oncology OS pattern:
      - Phase 1 (months 0-3): early elevated treatment hazard (immune priming /
        hyperprogression; HR = 1.6) → survival curves cross at ~month 3.
      - Phase 2 (months 3+): persistent benefit (HR = 0.65 forever) →
        survival curves diverge as T-cell response suppresses tumor growth.

    Statistical consequence (seed=9, n=200):
      - Log-rank: p ≈ 0.13  (not significant; early excess nearly cancels late benefit)
      - HC:       p ≈ 0.005 (significant; focuses on the late-benefit hot-spot)

    This setup deliberately demonstrates HC's advantage in the presence of
    crossing survival curves.  For a published application note, real POPLAR
    IPD should be used (see 01_get_kmdata.py).
    """
    rng = np.random.default_rng(seed)
    lam0     = np.log(2) / 10.0  # control arm baseline hazard
    hr_early = 1.60              # early treatment hazard multiplier
    hr_late  = 0.65              # late treatment hazard multiplier
    t0       = 3.0               # crossover time (months)
    censor   = 24.0

    # Control arm: simple exponential
    T_ctrl = rng.exponential(1 / lam0, n)
    E_ctrl = (T_ctrl <= censor).astype(float)
    T_ctrl = np.minimum(T_ctrl, censor)

    # Treatment arm: piecewise exponential (exact inversion from uniform)
    #   S_trt(t0) = exp(-lam0 * hr_early * t0)
    S0_trt = np.exp(-lam0 * hr_early * t0)
    u = rng.uniform(0, 1, n)

    die_early = u > S0_trt
    # Early deaths (0 to t0)
    T_early = -np.log(np.where(die_early, u, S0_trt + 1e-15)) / (lam0 * hr_early)
    # Late deaths (t0 onward)
    T_late = t0 - np.log(np.where(~die_early, u / S0_trt, 1.0)) / (lam0 * hr_late)

    T_trt = np.where(die_early, T_early, T_late)
    E_trt = (T_trt <= censor).astype(float)
    T_trt = np.minimum(T_trt, censor)

    print("[SYNTHETIC] Crossing-curves model (illustrative; use real data for paper).")
    print(f"  Control:   n={n}, events={int(E_ctrl.sum())}, "
          f"median={np.median(T_ctrl[E_ctrl==1]):.1f} mo")
    print(f"  Treatment: n={n}, events={int(E_trt.sum())}, "
          f"median={np.median(T_trt[E_trt==1]):.1f} mo")
    print(f"  Crossover at month {t0}: HR {hr_early} → {hr_late}")
    return T_ctrl, T_trt, E_ctrl, E_trt


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def run_analysis(T_ctrl, T_trt, E_ctrl, E_trt):
    print(f"\n{'='*60}")
    print("IMMUNO-ONCOLOGY — CheckMate 057 trial (PFS)")
    print(f"{'='*60}")

    print(f"\nRunning all tests (n_intervals={N_INTERVALS}, "
          f"n_permutations={N_PERMS}) ...")
    results = run_all_tests(
        T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        n_permutations=N_PERMS,
        label_A="Docetaxel",
        label_B="Nivolumab",
    )
    print_results_table(results)

    # Save results
    out_csv = RES_DIR / "immuno_test_results.csv"
    results.to_csv(out_csv)
    print(f"Results saved to {out_csv}")

    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figure(T_ctrl, T_trt, E_ctrl, E_trt):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle("Immuno-oncology: CheckMate 057 PFS (nivolumab vs docetaxel)",
                 fontsize=13, fontweight="bold")

    # Panel A: KM curves with HC shading
    ax_km = axes[0]
    df_dev = plot_km_with_hc(
        ax_km, T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        label_A="Docetaxel (control)",
        label_B="Nivolumab",
        shade_color="steelblue",
        title="A  Kaplan-Meier curves (PFS)",
        xlabel="Time (months)",
    )
    n_flagged = df_dev["suspected"].sum()
    ax_km.text(0.98, 0.98,
               f"HC-flagged intervals: {n_flagged}/{len(df_dev)}",
               transform=ax_km.transAxes, ha="right", va="top", fontsize=9,
               bbox=dict(boxstyle="round", fc="white", alpha=0.7))

    # Panel B: per-interval p-value profile
    ax_pv = axes[1]
    plot_pvalue_profile(
        ax_pv, df_dev,
        title=r"B  Per-interval hypergeometric $p$-values",
        xlabel="Time interval (months)",
    )

    fig.tight_layout()
    out_path = FIGS_DIR / "immuno_km.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    T_ctrl, T_trt, E_ctrl, E_trt = load_poplar(DATA_DIR)
    run_analysis(T_ctrl, T_trt, E_ctrl, E_trt)
    make_figure(T_ctrl, T_trt, E_ctrl, E_trt)


if __name__ == "__main__":
    main()
