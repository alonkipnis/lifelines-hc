"""
Domain 3: Targeted therapy in metastatic castration-resistant prostate cancer
           (COMET-1 trial, cabozantinib vs prednisone).

Biological context:
  Metastatic castration-resistant prostate cancer (mCRPC) predominantly
  metastasises to bone, making the bone microenvironment central to disease
  progression. Cabozantinib inhibits MET and VEGFR2, kinases that drive
  tumour survival and angiogenesis in the bone niche. In the COMET-1 trial
  (Smith et al. 2016, J Clin Oncol), cabozantinib significantly improved
  bone scan response at week 12 (secondary endpoint), confirming genuine
  biological activity against bone-metastatic disease. However, CRPC is
  sustained by multiple parallel survival pathways; once the targeted
  pathways are suppressed, resistance through alternative routes rapidly
  emerges. The net effect is a *temporally concentrated* hazard signal —
  windows of differential mortality during the bone-response phase —
  rather than a sustained proportional reduction in hazard. This structure
  dilutes the log-rank statistic (which averages over the entire follow-up)
  but is detectable by HC, which scans for any localised excess.

Statistical consequence (COMET-1, OS endpoint):
  Log-rank p = 0.262 (NS)   HCHG p ≈ 0.014 (significant)
  All four NPH-weighted alternatives: p ≥ 0.19 (all NS)

Dataset:
  COMET-1 (Smith et al. 2016).  n = 1028 (cabozantinib 682, prednisone 346).
  Endpoint: overall survival (OS). Follow-up ≈ 22 months.
  2:1 randomisation. Heavily pre-treated patients (post-docetaxel and
  abiraterone/enzalutamide).  IPD reconstructed via the Guyot algorithm
  as implemented in the kmdata R package.

Reference:
  Smith MR, et al. (2016). Phase III study of cabozantinib in previously
  treated metastatic castration-resistant prostate cancer: COMET-1.
  J Clin Oncol, 34(25), 3005-3013. DOI: 10.1200/JCO.2015.65.5597
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
N_INTERVALS = 80
N_PERMS     = 2000

for d in (FIGS_DIR, RES_DIR):
    d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def load_comet(data_dir: Path) -> tuple:
    """Load COMET-1 OS IPD from CSV; return (T_ctrl, T_trt, E_ctrl, E_trt)."""
    csv_path = data_dir / "COMET1_2A.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            "COMET1_2A.csv not found. "
            "Run '01_get_kmdata.py --trials COMET1_2A' to download."
        )

    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()

    ctrl = df[df["arm"] == "prednisone"]
    trt  = df[df["arm"] == "cabozantinib"]

    T_ctrl = ctrl["time"].values.astype(float)
    T_trt  = trt["time"].values.astype(float)
    E_ctrl = ctrl["event"].values.astype(float)
    E_trt  = trt["event"].values.astype(float)

    print(f"Loaded COMET1_2A.csv: {len(df)} rows")
    print(f"  Prednisone (control): n={len(T_ctrl)}, events={int(E_ctrl.sum())}, "
          f"median={np.median(T_ctrl):.1f} mo")
    print(f"  Cabozantinib:         n={len(T_trt)},  events={int(E_trt.sum())}, "
          f"median={np.median(T_trt):.1f} mo")
    return T_ctrl, T_trt, E_ctrl, E_trt


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def run_analysis(T_ctrl, T_trt, E_ctrl, E_trt):
    print(f"\n{'='*60}")
    print("TARGETED THERAPY (mCRPC) — COMET-1 trial (OS)")
    print(f"{'='*60}")
    print(f"\nRunning all tests (n_intervals={N_INTERVALS}, "
          f"n_permutations={N_PERMS}) ...")

    results = run_all_tests(
        T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        n_permutations=N_PERMS,
        label_A="Prednisone",
        label_B="Cabozantinib",
    )
    print_results_table(results)

    out_csv = RES_DIR / "comet_test_results.csv"
    results.to_csv(out_csv)
    print(f"Results saved to {out_csv}")
    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figure(T_ctrl, T_trt, E_ctrl, E_trt):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Targeted therapy in mCRPC: COMET-1 trial OS "
        "(cabozantinib vs prednisone)",
        fontsize=13, fontweight="bold",
    )

    ax_km = axes[0]
    df_dev = plot_km_with_hc(
        ax_km, T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        label_A="Prednisone (control)",
        label_B="Cabozantinib",
        shade_color="seagreen",
        title="A  Kaplan-Meier curves (OS)",
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
    out_path = FIGS_DIR / "comet_km.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    T_ctrl, T_trt, E_ctrl, E_trt = load_comet(DATA_DIR)
    run_analysis(T_ctrl, T_trt, E_ctrl, E_trt)
    make_figure(T_ctrl, T_trt, E_ctrl, E_trt)


if __name__ == "__main__":
    main()
