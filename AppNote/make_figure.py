"""
Generate the composite Application Note figure (Figure 1).

Layout (2 rows × 3 columns):
  Row 1 — KM curves with HC-flagged intervals
    Panel A: Immuno-oncology    (CheckMate 057 PFS)
    Panel B: Bisphosphonate     (AZURE trial DFS)
    Panel C: Targeted therapy   (COMET-1 trial OS)
  Row 2 — Per-interval -log10(p) bar charts
    Panel D: CheckMate 057
    Panel E: AZURE
    Panel F: COMET-1

Usage:
    python make_figure.py [--out figs/figure1.png]

Pre-requisites:
  run_immuno_oncology.py  (or Checkmate057_1C.csv in data/)
  run_azure.py            (or AZURE_2A.csv in data/)
  run_comet.py            (or COMET1_2A.csv in data/)
"""

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

mpl.style.use("ggplot")
plt.rcParams.update({
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 9,
})

SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))
from utils import (
    run_all_tests, pvalue_profile,
    plot_km_with_hc, plot_pvalue_profile,
)

DATA_DIR = SCRIPT_DIR / "data"
FIGS_DIR = SCRIPT_DIR / "figs"
FIGS_DIR.mkdir(parents=True, exist_ok=True)

N_INTERVALS_IO    = 60   # immuno-oncology  (~1 bin/month for 30-mo PFS)
N_INTERVALS_AZURE = 80   # bisphosphonate   (~1.5-month bins, 120-mo follow-up)
N_INTERVALS_COMET = 80   # targeted therapy (~0.3-month bins, 22-mo OS)
N_INTERVALS_CSL   = 50   # corticosteroid   (~0.25-year bins, 13-yr OS)
N_PERMS           = 2000


# ---------------------------------------------------------------------------
# Data loaders
# ---------------------------------------------------------------------------

def load_poplar():
    """Load CheckMate 057 PFS data or synthetic fallback."""
    for fname in ["Checkmate057_1C.csv", "POPLAR.csv"]:
        csv_path = DATA_DIR / fname
        if csv_path.exists():
            df = pd.read_csv(csv_path)
            df.columns = df.columns.str.strip().str.lower()
            time_col  = _pick(df, ["time", "os", "os_time", "t"])
            event_col = _pick(df, ["status", "event", "os_event", "dead"])
            arm_col   = _pick(df, ["arm", "trt", "treatment", "group"])
            if time_col and event_col:
                df[time_col]  = pd.to_numeric(df[time_col],  errors="coerce")
                df[event_col] = pd.to_numeric(df[event_col], errors="coerce")
                df = df.dropna(subset=[time_col, event_col])
                if arm_col:
                    vals = sorted(df[arm_col].unique())
                    df["_arm"] = df[arm_col].map({vals[0]: 0, vals[1]: 1}).fillna(0)
                else:
                    df["_arm"] = [0]*(len(df)//2) + [1]*(len(df)-len(df)//2)
                ctrl = df[df["_arm"]==0]; trt = df[df["_arm"]==1]
                return (ctrl[time_col].values, trt[time_col].values,
                        ctrl[event_col].values, trt[event_col].values)

    print("[INFO] Trial CSV not found; using synthetic crossover-curve fallback.")
    rng = np.random.default_rng(9)
    n, lam0, hr_early, hr_late, t0, censor = 200, np.log(2)/10, 1.60, 0.65, 3.0, 24.0
    T_ctrl = rng.exponential(1/lam0, n)
    E_ctrl = (T_ctrl <= censor).astype(float); T_ctrl = np.minimum(T_ctrl, censor)
    S0 = np.exp(-lam0 * hr_early * t0)
    u = rng.uniform(0, 1, n)
    die_early = u > S0
    T_ph1 = -np.log(u) / (lam0 * hr_early)
    T_ph2 = t0 - np.log(np.clip(u / S0, 1e-300, 1.0)) / (lam0 * hr_late)
    T_trt = np.where(die_early, T_ph1, T_ph2)
    E_trt = (T_trt <= censor).astype(float); T_trt = np.minimum(T_trt, censor)
    return T_ctrl, T_trt, E_ctrl, E_trt


def load_azure():
    """Load AZURE trial DFS data or synthetic delayed-benefit fallback."""
    csv_path = DATA_DIR / "AZURE_2A.csv"
    if csv_path.exists():
        df = pd.read_csv(csv_path)
        df.columns = df.columns.str.strip().str.lower()
        ctrl = df[df["arm"] == "control"]
        trt  = df[df["arm"] == "zoledronic_acid"]
        return (ctrl["time"].values.astype(float), trt["time"].values.astype(float),
                ctrl["event"].values.astype(float), trt["event"].values.astype(float))

    print("[INFO] AZURE_2A.csv not found; using synthetic delayed-benefit fallback.")
    rng  = np.random.default_rng(7)
    n, lam0, t0, hr, censor = 1200, np.log(2)/80, 24.0, 0.75, 120.0
    T_ctrl = rng.exponential(1/lam0, n)
    E_ctrl = (T_ctrl <= censor).astype(float); T_ctrl = np.minimum(T_ctrl, censor)
    S_t0 = np.exp(-lam0 * t0); u = rng.uniform(0, 1, n)
    die_early = u > S_t0
    T_trt = np.where(die_early, -np.log(np.clip(u, 1e-300, 1))/lam0,
                     t0 - np.log(np.clip(u/S_t0, 1e-300, 1))/(lam0*hr))
    E_trt = (T_trt <= censor).astype(float); T_trt = np.minimum(T_trt, censor)
    return T_ctrl, T_trt, E_ctrl, E_trt


def load_comet():
    """Load COMET-1 trial OS data."""
    csv_path = DATA_DIR / "COMET1_2A.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            "COMET1_2A.csv not found. "
            "Run '01_get_kmdata.py --trials COMET1_2A' to download."
        )
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()
    ctrl = df[df["arm"] == "prednisone"]
    trt  = df[df["arm"] == "cabozantinib"]
    return (ctrl["time"].values.astype(float), trt["time"].values.astype(float),
            ctrl["event"].values.astype(float), trt["event"].values.astype(float))


def load_csl():
    """Load CSL liver-cirrhosis trial OS data (real public IPD)."""
    csv_path = DATA_DIR / "CSL.csv"
    if not csv_path.exists():
        raise FileNotFoundError(
            "CSL.csv not found. Run 'run_csl.py' once (it can regenerate "
            "from the installed SurvSet package)."
        )
    df = pd.read_csv(csv_path)
    df.columns = df.columns.str.strip().str.lower()
    ctrl = df[df["arm"] == "placebo"]
    trt  = df[df["arm"] == "prednisone"]
    return (ctrl["time"].values.astype(float), trt["time"].values.astype(float),
            ctrl["event"].values.astype(float), trt["event"].values.astype(float))


def _pick(df, candidates):
    for c in candidates:
        if c in df.columns: return c
    for c in candidates:
        for col in df.columns:
            if c in col: return col
    return None


# ---------------------------------------------------------------------------
# Figure
# ---------------------------------------------------------------------------

def build_figure(out_path: Path):
    T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt = load_poplar()
    T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt = load_azure()
    T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt = load_comet()

    dev_io = pvalue_profile(T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
                             n_intervals=N_INTERVALS_IO)
    dev_az = pvalue_profile(T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
                             n_intervals=N_INTERVALS_AZURE)
    dev_co = pvalue_profile(T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
                             n_intervals=N_INTERVALS_COMET)

    res_io = run_all_tests(T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
                            n_intervals=N_INTERVALS_IO, n_permutations=N_PERMS)
    res_az = run_all_tests(T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
                            n_intervals=N_INTERVALS_AZURE, n_permutations=N_PERMS)
    res_co = run_all_tests(T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
                            n_intervals=N_INTERVALS_COMET, n_permutations=N_PERMS)

    # 2 rows × 3 columns: row 1 = KM plots, row 2 = p-value profiles
    fig = plt.figure(figsize=(18, 9))
    gs = GridSpec(2, 3, figure=fig, hspace=0.42, wspace=0.32)

    ax_A = fig.add_subplot(gs[0, 0])
    ax_B = fig.add_subplot(gs[0, 1])
    ax_C = fig.add_subplot(gs[0, 2])
    ax_D = fig.add_subplot(gs[1, 0])
    ax_E = fig.add_subplot(gs[1, 1])
    ax_F = fig.add_subplot(gs[1, 2])

    # --- Row 1: KM plots ---
    plot_km_with_hc(
        ax_A, T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
        n_intervals=N_INTERVALS_IO,
        label_A="Docetaxel (control)",
        label_B="Nivolumab",
        shade_color="steelblue",
        title="A  Immuno-oncology (CheckMate 057 PFS)",
        xlabel="Time (months)",
    )
    _annotate_pvals(ax_A, res_io)

    plot_km_with_hc(
        ax_B, T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
        n_intervals=N_INTERVALS_AZURE,
        label_A="Control",
        label_B="Zoledronic acid",
        shade_color="darkorange",
        title="B  Adjuvant bisphosphonate (AZURE DFS)",
        xlabel="Time (months)",
    )
    _annotate_pvals(ax_B, res_az)

    plot_km_with_hc(
        ax_C, T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
        n_intervals=N_INTERVALS_COMET,
        label_A="Prednisone (control)",
        label_B="Cabozantinib",
        shade_color="seagreen",
        title="C  Targeted therapy (COMET-1 OS)",
        xlabel="Time (months)",
    )
    _annotate_pvals(ax_C, res_co)

    # --- Row 2: p-value profiles ---
    plot_pvalue_profile(
        ax_D, dev_io,
        title=r"D  Interval $p$-values — CheckMate 057",
        xlabel="Time interval (months)",
    )
    plot_pvalue_profile(
        ax_E, dev_az,
        title=r"E  Interval $p$-values — AZURE",
        xlabel="Time interval (months)",
    )
    plot_pvalue_profile(
        ax_F, dev_co,
        title=r"F  Interval $p$-values — COMET-1",
        xlabel="Time interval (months)",
    )

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"\nFigure 1 saved to {out_path}")
    plt.close(fig)


def build_figure_4domain(out_path: Path):
    """4-domain composite (adds CSL). Used by the bullet-outline manuscript;
    the 3-domain build_figure() output is left untouched for main.tex."""
    T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt = load_poplar()
    T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt = load_comet()
    T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt = load_azure()
    T_cs_ctrl, T_cs_trt, E_cs_ctrl, E_cs_trt = load_csl()

    dev_io = pvalue_profile(T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
                             n_intervals=N_INTERVALS_IO)
    dev_co = pvalue_profile(T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
                             n_intervals=N_INTERVALS_COMET)
    dev_az = pvalue_profile(T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
                             n_intervals=N_INTERVALS_AZURE)
    dev_cs = pvalue_profile(T_cs_ctrl, T_cs_trt, E_cs_ctrl, E_cs_trt,
                             n_intervals=N_INTERVALS_CSL)

    res_io = run_all_tests(T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
                            n_intervals=N_INTERVALS_IO, n_permutations=N_PERMS)
    res_co = run_all_tests(T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
                            n_intervals=N_INTERVALS_COMET, n_permutations=N_PERMS)
    res_az = run_all_tests(T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
                            n_intervals=N_INTERVALS_AZURE, n_permutations=N_PERMS)
    res_cs = run_all_tests(T_cs_ctrl, T_cs_trt, E_cs_ctrl, E_cs_trt,
                            n_intervals=N_INTERVALS_CSL, n_permutations=N_PERMS)

    # 2 rows × 4 columns: row 1 = KM plots (A-D), row 2 = p-value profiles (E-H)
    fig = plt.figure(figsize=(22, 9))
    gs = GridSpec(2, 4, figure=fig, hspace=0.42, wspace=0.34)

    ax = [fig.add_subplot(gs[r, c]) for r in range(2) for c in range(4)]
    ax_A, ax_B, ax_C, ax_D, ax_E, ax_F, ax_G, ax_H = ax

    # --- Row 1: KM plots ---
    plot_km_with_hc(ax_A, T_io_ctrl, T_io_trt, E_io_ctrl, E_io_trt,
                    n_intervals=N_INTERVALS_IO,
                    label_A="Docetaxel (control)", label_B="Nivolumab",
                    shade_color="steelblue",
                    title="A  Immuno-oncology (CheckMate 057 PFS)",
                    xlabel="Time (months)")
    _annotate_pvals(ax_A, res_io)

    plot_km_with_hc(ax_B, T_co_ctrl, T_co_trt, E_co_ctrl, E_co_trt,
                    n_intervals=N_INTERVALS_COMET,
                    label_A="Prednisone (control)", label_B="Cabozantinib",
                    shade_color="seagreen",
                    title="B  Targeted therapy (COMET-1 OS)",
                    xlabel="Time (months)")
    _annotate_pvals(ax_B, res_co)

    plot_km_with_hc(ax_C, T_az_ctrl, T_az_trt, E_az_ctrl, E_az_trt,
                    n_intervals=N_INTERVALS_AZURE,
                    label_A="Control", label_B="Zoledronic acid",
                    shade_color="darkorange",
                    title="C  Adjuvant bisphosphonate (AZURE DFS)",
                    xlabel="Time (months)")
    _annotate_pvals(ax_C, res_az)

    plot_km_with_hc(ax_D, T_cs_ctrl, T_cs_trt, E_cs_ctrl, E_cs_trt,
                    n_intervals=N_INTERVALS_CSL,
                    label_A="Placebo (control)", label_B="Prednisone",
                    shade_color="mediumorchid",
                    title="D  Corticosteroid (CSL cirrhosis OS)",
                    xlabel="Time (years)")
    _annotate_pvals(ax_D, res_cs)

    # --- Row 2: p-value profiles ---
    plot_pvalue_profile(ax_E, dev_io,
                        title=r"E  Interval $p$-values — CheckMate 057",
                        xlabel="Time interval (months)")
    plot_pvalue_profile(ax_F, dev_co,
                        title=r"F  Interval $p$-values — COMET-1",
                        xlabel="Time interval (months)")
    plot_pvalue_profile(ax_G, dev_az,
                        title=r"G  Interval $p$-values — AZURE",
                        xlabel="Time interval (months)")
    plot_pvalue_profile(ax_H, dev_cs,
                        title=r"H  Interval $p$-values — CSL",
                        xlabel="Time interval (years)")

    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    print(f"\n4-domain figure saved to {out_path}")
    plt.close(fig)


def _annotate_pvals(ax, results: pd.DataFrame) -> None:
    lr_p = results.loc["Log-rank", "p_value"] \
           if "Log-rank" in results.index else float("nan")
    hc_p = results.loc["Higher Criticism (HC)", "p_value"] \
           if "Higher Criticism (HC)" in results.index else float("nan")
    text = f"Log-rank  p = {lr_p:.3f}\nHCHG      p = {hc_p:.3f}"
    ax.text(0.98, 0.42, text,
            transform=ax.transAxes, ha="right", va="top",
            fontsize=8.5, family="monospace",
            bbox=dict(boxstyle="round", fc="lightyellow", alpha=0.9))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate Application Note Figure 1")
    parser.add_argument("--out", default=str(FIGS_DIR / "figure1.png"),
                        help="Output file path for the 3-domain figure (main.tex)")
    parser.add_argument("--four-domain", action="store_true",
                        help="Also generate the 4-domain figure (adds CSL) "
                             "for the bullet-outline manuscript")
    parser.add_argument("--four-out", default=str(FIGS_DIR / "figure1_4dom.png"),
                        help="Output path for the 4-domain figure")
    args = parser.parse_args()
    build_figure(Path(args.out))
    if args.four_domain:
        build_figure_4domain(Path(args.four_out))


if __name__ == "__main__":
    main()
