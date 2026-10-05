"""
Domain 4: Corticosteroid therapy in liver cirrhosis (CSL trial).

This is the one domain in the Application Note whose raw individual patient
data (IPD) is *genuinely publicly available* — not reconstructed from a
published Kaplan-Meier curve. The data ship with the CRAN ``timereg``
package and are redistributed in the open SurvSet repository
(Drysdale 2022).

Biological / clinical context:
  The Copenhagen Study group for Liver diseases (CSL) ran a randomized
  controlled trial (1962-1969) of prednisone versus placebo in patients
  with histologically verified liver cirrhosis. Corticosteroids have a
  *time-varying* effect in cirrhosis: an early period of possible harm
  (immunosuppression, infection risk, fluid retention) is followed by a
  later period of anti-inflammatory benefit in the subset of patients
  whose cirrhosis has an active inflammatory component. The net result is
  a survival curve that crosses more than once, with the treatment effect
  concentrated in a few, non-contiguous time windows rather than a
  constant hazard ratio. This is a canonical non-proportional-hazards
  dataset, widely used to illustrate time-varying-coefficient models.

Statistical consequence:
  Log-rank p = 0.384 (NS); all four weighted alternatives p >= 0.13 (NS);
  HCHG p ~ 0.002 (significant). No single fixed temporal weight captures
  the multi-window, sign-changing departure, but HCHG's interval scan does.

Data:
  CSL liver cirrhosis trial. n = 446 patients (226 prednisone, 220
  placebo), 270 deaths, follow-up up to ~13 years. Endpoint: overall
  survival. The SurvSet CSV is in counting-process (long) format with a
  time-varying prothrombin covariate; we collapse to one row per patient
  (final follow-up time and event) for the two-sample survival comparison.

References:
  Copenhagen Study group for Liver diseases. Prednisone versus placebo in
  cirrhosis of the liver (CSL-1 trial), 1962-1969.
  Martinussen T, Scheike TH (2006). Dynamic Regression Models for Survival
  Data. Springer (the ``timereg`` R package, dataset ``csl``).
  Drysdale E (2022). SurvSet: an open-source time-to-event dataset
  repository.
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
N_INTERVALS = 50
N_PERMS     = 2000

for d in (FIGS_DIR, RES_DIR):
    d.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Data loading
# ---------------------------------------------------------------------------

def _collapse_survset_csl(df: pd.DataFrame) -> pd.DataFrame:
    """Collapse SurvSet counting-process CSL data to one row per patient."""
    tcol = "time2" if "time2" in df.columns else "time"
    df = df.sort_values(["pid", tcol])
    last = df.groupby("pid").tail(1).copy()
    out = pd.DataFrame({
        "time":  last[tcol].astype(float).values,
        "event": last["event"].astype(int).values,
        "arm":   last["fac_treat"].astype(str).values,
    })
    return out[out["time"] > 0].reset_index(drop=True)


def load_csl(data_dir: Path) -> tuple:
    """Load CSL cirrhosis IPD; return (T_ctrl, T_trt, E_ctrl, E_trt).

    Prefers the pre-collapsed AppNote/data/CSL.csv. Falls back to the raw
    SurvSet distribution (counting-process format) installed in the
    environment, collapsing it on the fly.
    """
    csv_path = data_dir / "CSL.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip().str.lower()
        print(f"Loaded CSL.csv: {len(df)} patients")
    else:
        try:
            import SurvSet  # noqa: F401
            survset_csv = (Path(SurvSet.__file__).parent
                           / "_datagen" / "output" / "csl.csv")
            df = _collapse_survset_csl(pd.read_csv(survset_csv))
            print(f"Loaded CSL from SurvSet (collapsed): {len(df)} patients")
        except Exception as e:  # pragma: no cover
            raise FileNotFoundError(
                "CSL.csv not found and SurvSet fallback failed. "
                "Install SurvSet (`pip install SurvSet`) or place CSL.csv "
                f"in {data_dir}. Original error: {e}"
            )

    ctrl = df[df["arm"] == "placebo"]
    trt  = df[df["arm"] == "prednisone"]

    T_ctrl = ctrl["time"].values.astype(float)
    T_trt  = trt["time"].values.astype(float)
    E_ctrl = ctrl["event"].values.astype(float)
    E_trt  = trt["event"].values.astype(float)

    print(f"  Placebo:    n={len(T_ctrl)}, events={int(E_ctrl.sum())}, "
          f"median={np.median(T_ctrl):.1f} yr")
    print(f"  Prednisone: n={len(T_trt)},  events={int(E_trt.sum())}, "
          f"median={np.median(T_trt):.1f} yr")
    return T_ctrl, T_trt, E_ctrl, E_trt


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def run_analysis(T_ctrl, T_trt, E_ctrl, E_trt):
    print(f"\n{'='*60}")
    print("CORTICOSTEROID THERAPY (cirrhosis) — CSL trial (OS)")
    print(f"{'='*60}")
    print(f"\nRunning all tests (n_intervals={N_INTERVALS}, "
          f"n_permutations={N_PERMS}) ...")

    results = run_all_tests(
        T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        n_permutations=N_PERMS,
        label_A="Placebo",
        label_B="Prednisone",
    )
    print_results_table(results)

    out_csv = RES_DIR / "csl_test_results.csv"
    results.to_csv(out_csv)
    print(f"Results saved to {out_csv}")
    return results


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def make_figure(T_ctrl, T_trt, E_ctrl, E_trt):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    fig.suptitle(
        "Corticosteroid therapy in liver cirrhosis: CSL trial "
        "(prednisone vs placebo)",
        fontsize=13, fontweight="bold",
    )

    ax_km = axes[0]
    df_dev = plot_km_with_hc(
        ax_km, T_ctrl, T_trt, E_ctrl, E_trt,
        n_intervals=N_INTERVALS,
        label_A="Placebo (control)",
        label_B="Prednisone",
        shade_color="mediumorchid",
        title="A  Kaplan-Meier curves (OS)",
        xlabel="Time (years)",
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
        xlabel="Time interval (years)",
    )

    fig.tight_layout()
    out_path = FIGS_DIR / "csl_km.png"
    fig.savefig(out_path, dpi=180, bbox_inches="tight")
    print(f"\nFigure saved to {out_path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    T_ctrl, T_trt, E_ctrl, E_trt = load_csl(DATA_DIR)
    run_analysis(T_ctrl, T_trt, E_ctrl, E_trt)
    make_figure(T_ctrl, T_trt, E_ctrl, E_trt)


if __name__ == "__main__":
    main()
