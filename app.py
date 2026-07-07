"""
Streamlit web app for Higher Criticism survival analysis.

Run with:
    streamlit run app.py
"""

import io
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st

from lifelines.statistics import logrank_test
from lifelines_hc import (
    higher_criticism_test,
    suspected_deviations,
    KaplanMeierHCIllustrator,
)

# ── Page config ──────────────────────────────────────────────────────────

st.set_page_config(
    page_title="Higher Criticism — Survival Analysis",
    page_icon="📈",
    layout="wide",
)

SAMPLE_CSV = "Data/SCANB_demo.csv"


# ── Helpers ──────────────────────────────────────────────────────────────

@st.cache_data(show_spinner=False)
def load_csv(uploaded_file):
    return pd.read_csv(uploaded_file)


@st.cache_data(show_spinner=False)
def load_sample():
    import os
    path = os.path.join(os.path.dirname(__file__), "..", SAMPLE_CSV)
    if os.path.exists(path):
        return pd.read_csv(path)
    return None


def extract_groups(df, time_col, event_col, group_col):
    """Split a DataFrame into two groups based on a binary column."""
    mask_a = df[group_col] == 0
    mask_b = df[group_col] == 1
    return (
        df.loc[mask_a, time_col].values.astype(float),
        df.loc[mask_b, time_col].values.astype(float),
        df.loc[mask_a, event_col].values.astype(float),
        df.loc[mask_b, event_col].values.astype(float),
    )


def style_suspected(row):
    """Highlight suspected rows in the dataframe."""
    if row.get("suspected", False):
        return ["background-color: rgba(255, 200, 200, 0.4)"] * len(row)
    return [""] * len(row)


# ── Sidebar ──────────────────────────────────────────────────────────────

st.sidebar.title("Parameters")

st.sidebar.subheader("Data source")
uploaded = st.sidebar.file_uploader(
    "Upload CSV", type=["csv"],
    help=(
        "CSV with at least a **time** column, an **event** column "
        "(1 = event, 0 = censored), and one or more binary group columns. "
        "Or use the built-in SCANB demo dataset."
    ),
)

use_sample = False
if uploaded is not None:
    df = load_csv(uploaded)
    st.sidebar.success(f"Active: **{uploaded.name}**")
elif load_sample() is not None:
    df = load_sample()
    use_sample = True
    st.sidebar.info("Active: **SCANB_demo** (built-in)")
else:
    st.info("Upload a CSV file or place the SCANB dataset in `Data/`.")
    st.stop()

all_cols = list(df.columns)

st.sidebar.markdown("---")
st.sidebar.subheader("Column mapping")

time_col = st.sidebar.selectbox(
    "Time column",
    options=[c for c in all_cols if c.lower() in ("time", "duration", "t")],
    index=0,
) if any(c.lower() in ("time", "duration", "t") for c in all_cols) else st.sidebar.selectbox(
    "Time column", options=all_cols,
)

event_col = st.sidebar.selectbox(
    "Event column",
    options=[c for c in all_cols if c.lower() in ("event", "status", "e")],
    index=0,
) if any(c.lower() in ("event", "status", "e") for c in all_cols) else st.sidebar.selectbox(
    "Event column", options=all_cols,
)

candidate_groups = [
    c for c in all_cols
    if c not in (time_col, event_col)
    and df[c].dropna().isin([0, 1, 0.0, 1.0]).all()
]

if not candidate_groups:
    st.error(
        "No binary (0/1) group columns found. "
        "Make sure the CSV contains at least one column with only 0 and 1 values."
    )
    st.stop()

group_col = st.sidebar.selectbox("Group column", options=candidate_groups)
swap_groups = st.sidebar.checkbox(
    "Swap A ↔ B",
    help="By default A = 0, B = 1. Check to swap so A = 1, B = 0. "
         "The test detects increased hazard from A to B.",
)

st.sidebar.markdown("---")
st.sidebar.subheader("HC test settings")

n_intervals = st.sidebar.slider(
    "Intervals to pool (`n_intervals_to_pool`)",
    min_value=10, max_value=2000, value=70, step=10,
    help="Number of equal-width time bins. Higher = finer resolution.",
)
gamma = st.sidebar.slider(
    "Gamma (γ)", min_value=0.05, max_value=0.50, value=0.30, step=0.05,
    help="HC fraction parameter — only top γ·n ordered p-values are considered.",
)
alternative = st.sidebar.selectbox(
    "Alternative", options=["greater", "less", "both"], index=0,
    help=(
        "**greater**: excess hazard in group 1.  "
        "**less**: excess hazard in group 0.  "
        "**both**: maximum of both directions."
    ),
)
stbl = st.sidebar.checkbox("Stabilised HC", value=True)

t_0_option = st.sidebar.checkbox("Restrict time range (t₀)", value=False)
t_0 = -1.0
if t_0_option:
    t_0 = st.sidebar.number_input("t₀", min_value=0.0, value=30.0, step=1.0)

st.sidebar.markdown("---")
st.sidebar.subheader("Permutation test")

n_perms = st.sidebar.number_input(
    "Number of permutations",
    min_value=0, max_value=20000, value=0, step=1000,
    help="Set > 0 to compute a permutation p-value (slower).",
)
seed = st.sidebar.number_input("Random seed", min_value=0, value=42, step=1)

# ── Main panel ───────────────────────────────────────────────────────────

st.title("Higher Criticism for Survival Analysis with Rare, Non-Proportional Hazards Departures")
st.markdown("**Instructions:** Select parameters on the sidebar to the left.")
st.caption("Based on the paper: "
    "Kipnis, A., Galili, B., and Yakhini, Z. "
    "[*Higher criticism for rare and weak non-proportional hazard deviations "
    "in survival analysis.*](https://academic.oup.com/biomet/article/113/1/asaf075/8307530) "
    "Biometrika, Volume 113, Issue 1, 2026."
)

if use_sample:
    st.info(f"Using built-in **SCANB** dataset • grouping by **{group_col}**")

T_A, T_B, event_A, event_B = extract_groups(df, time_col, event_col, group_col)
if swap_groups:
    T_A, T_B = T_B, T_A
    event_A, event_B = event_B, event_A

val_A, val_B = (1, 0) if swap_groups else (0, 1)
label_A = f"A ({group_col}={val_A})"
label_B = f"B ({group_col}={val_B})"

n_A, n_B = len(T_A), len(T_B)
st.markdown(
    f"**{label_A}**: {n_A} subjects &emsp; "
    f"**{label_B}**: {n_B} subjects"
)

if n_A < 2 or n_B < 2:
    st.error("Each group must have at least 2 subjects.")
    st.stop()

# ── Compute ──────────────────────────────────────────────────────────────

with st.spinner("Computing HC test…"):
    hc_result = higher_criticism_test(
        T_A, T_B, event_A, event_B,
        n_intervals_to_pool=n_intervals,
        gamma=gamma, alternative=alternative, stbl=stbl,
        t_0=t_0, n_permutations=int(n_perms), seed=int(seed),
    )

    lr_result = logrank_test(T_A, T_B, event_A, event_B)

    df_dev = suspected_deviations(
        T_A, T_B, event_A, event_B,
        n_intervals_to_pool=n_intervals,
        gamma=gamma, alternative=alternative, stbl=stbl, t_0=t_0,
    )

# ── Summary metrics ──────────────────────────────────────────────────────

col1, col2, col3, col4 = st.columns(4)
col1.metric("HC statistic", f"{hc_result.test_statistic:.3f}")
if np.isfinite(hc_result.p_value):
    col2.metric("HC p-value (perm)", f"{hc_result.p_value:.4f}")
else:
    col2.metric("HC p-value (perm)", "—")
col3.metric("Log-rank statistic", f"{lr_result.test_statistic:.3f}")
col4.metric("Log-rank p-value", f"{lr_result.p_value:.4g}")

perm_stats = getattr(hc_result, "permutation_statistics", None)
if perm_stats is not None:
    fig_perm, ax_perm = plt.subplots(figsize=(4, 1.6))
    KaplanMeierHCIllustrator.plot_permutation_histogram(hc_result, ax=ax_perm)
    fig_perm.tight_layout()
    st.pyplot(fig_perm)
    plt.close(fig_perm)

# ── KM plot ──────────────────────────────────────────────────────────────

if np.isfinite(hc_result.p_value) and hc_result.p_value > 0.05:
    st.warning(
        "The HC statistic is **not significant** "
        f"(p = {hc_result.p_value:.4f}). "
        "Suspected intervals should be disregarded."
    )
elif not np.isfinite(hc_result.p_value) and int(n_perms) == 0:
    st.info(
        "No permutation test was run — significance of the HC statistic is unknown. "
        "Set **Number of permutations > 0** in the sidebar to assess significance."
    )

st.subheader("Kaplan-Meier curves")

ill = KaplanMeierHCIllustrator(
    T_A, T_B, event_A, event_B,
    label_A=label_A,
    label_B=label_B,
)
title = ill.format_title(hc_result)

fig, ax = plt.subplots(figsize=(10, 5))
ill.plot(df_dev, ax=ax, title=title)
fig.tight_layout()
st.pyplot(fig)
plt.close(fig)

# ── Interval table ───────────────────────────────────────────────────────

st.subheader("Per-interval results")

n_suspected = df_dev["suspected"].sum()
st.markdown(
    f"**{n_suspected}** of **{len(df_dev)}** intervals suspected of having greater hazard in {label_B} compared to {label_A}"
    f"(HC threshold = {df_dev['hc_threshold'].iloc[0]:.4f})"
)

display_cols = [
    c for c in (
        "at_risk_A", "at_risk_B", "observed_A", "observed_B",
        "hypergeom_pvalue", "suspected",
    ) if c in df_dev.columns
]
show_df = df_dev[display_cols].copy()
show_df.index.name = "time_interval"

st.dataframe(
    show_df.style.apply(style_suspected, axis=1).format(
        {
            "hypergeom_pvalue": "{:.4f}",
            "at_risk_A": "{:.0f}",
            "at_risk_B": "{:.0f}",
            "observed_A": "{:.0f}",
            "observed_B": "{:.0f}",
        }
    ),
    use_container_width=True,
    height=500,
)

# ── Download ─────────────────────────────────────────────────────────────

csv_buf = io.StringIO()
df_dev.to_csv(csv_buf)
st.download_button(
    "Download full results CSV",
    data=csv_buf.getvalue(),
    file_name=f"hc_results_{group_col}.csv",
    mime="text/csv",
)
