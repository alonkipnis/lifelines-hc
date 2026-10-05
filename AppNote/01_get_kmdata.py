"""
Download and prepare clinical trial data from the kmdata R package.

The kmdata package (github.com/raredd/kmdata) contains reconstructed individual
patient-level data (IPD) from 300+ phase III oncology trials, reverse-engineered
from published Kaplan-Meier curves using the Guyot algorithm.

Trials downloaded by default:
  Checkmate057_1C  - Nivolumab vs docetaxel, NSCLC PFS (Borghaei 2015 NEJM)
                     Crossing PFS curves: LR p=0.35, HC p=0.002
  Checkmate057_1A  - Same trial, OS endpoint (for reference)
  AZURE_2A         - Zoledronic acid vs control, breast cancer DFS
                     (Coleman 2011 NEJM); menopause-dependent effect:
                     LR p=0.30, HC p=0.012

Output:  AppNote/data/<TRIAL>.csv
Columns: time, event, arm  (standardised by kmdata)

Usage:
    python 01_get_kmdata.py [--trials Checkmate057_1C AZURE_2A]

Requires:
    R (>= 4.0) accessible as 'Rscript'
    Internet access (for kmdata package installation from GitHub)
"""

import argparse
import subprocess
import sys
from pathlib import Path

DATA_DIR = Path(__file__).parent / "data"

# Primary trial: CheckMate 057 PFS (crossing curves, ideal HC showcase)
DEFAULT_TRIALS = [
    "Checkmate057_1C",   # PFS, nivolumab vs docetaxel (LR NS, HC sig)
    "Checkmate057_1A",   # OS, nivolumab vs docetaxel (both significant)
    "AZURE_2A",          # DFS, zoledronic acid vs control (LR NS, HC sig)
]

R_EXPORT_SCRIPT = r"""
suppressPackageStartupMessages({
  if (!require("remotes", quietly = TRUE)) {
    install.packages("remotes", repos = "https://cloud.r-project.org", quiet = TRUE)
  }
  if (!require("kmdata", quietly = TRUE)) {
    cat("Installing kmdata from GitHub (raredd/kmdata)...\n")
    remotes::install_github("raredd/kmdata", quiet = TRUE)
  }
  library(kmdata)
})

trials  <- commandArgs(trailingOnly = TRUE)
out_dir <- trials[length(trials)]   # last arg is output dir
trials  <- trials[-length(trials)]  # remaining args are trial names

all_names <- ls("package:kmdata")
cat(sprintf("kmdata contains %d datasets.\n", length(all_names)))
cat(sprintf("Requested: %s\n", paste(trials, collapse = ", ")))

for (trial in trials) {
  # Exact or fuzzy match against package objects
  matched <- trial
  if (!matched %in% all_names) {
    pat <- gsub("_", ".", trial, fixed = TRUE)
    fuzzy <- all_names[grepl(pat, all_names, ignore.case = TRUE)]
    if (length(fuzzy) == 0) {
      cat(sprintf("[SKIP] '%s' not found. Fuzzy search returned nothing.\n", trial))
      next
    }
    matched <- fuzzy[1]
    cat(sprintf("[INFO] Matched '%s' -> '%s'.\n", trial, matched))
  }
  df <- get(matched, envir = asNamespace("kmdata"))
  if (!is.data.frame(df)) {
    cat(sprintf("[SKIP] '%s' is not a data frame.\n", matched))
    next
  }
  fn <- file.path(out_dir, paste0(trial, ".csv"))
  write.csv(df, fn, row.names = FALSE)
  cat(sprintf("[OK]   %s -> %s  (%d rows, cols: %s)\n",
              matched, fn, nrow(df), paste(colnames(df), collapse = ", ")))
}
"""


def run_r_export(trials: list[str], out_dir: Path) -> bool:
    """Call Rscript to install kmdata and export trials to CSV.

    Returns True on success.
    """
    out_dir.mkdir(parents=True, exist_ok=True)

    # Write R script to a temp file (avoids shell quoting issues)
    r_script_path = out_dir / "_export_kmdata.R"
    r_script_path.write_text(R_EXPORT_SCRIPT)

    cmd = ["Rscript", "--no-save", str(r_script_path)] + trials + [str(out_dir)]
    print("Running R export (may take a few minutes on first run)...")
    result = subprocess.run(cmd, text=True, capture_output=True)

    print(result.stdout)
    if result.returncode != 0:
        print("R stderr:", result.stderr[-2000:], file=sys.stderr)
        return False
    return True


def check_output(trials: list[str], out_dir: Path) -> None:
    """Print a summary of exported files."""
    import pandas as pd

    print("\n--- Downloaded files ---")
    for trial in trials:
        fp = out_dir / f"{trial}.csv"
        if fp.exists():
            df = pd.read_csv(fp)
            print(f"  {trial}: {len(df)} patients, columns: {list(df.columns)}")
        else:
            print(f"  {trial}: NOT FOUND")


def main():
    parser = argparse.ArgumentParser(description="Download kmdata trial IPD to CSV")
    parser.add_argument("--trials", nargs="+", default=DEFAULT_TRIALS,
                        help="Trial names (as in kmdata R package)")
    parser.add_argument("--out-dir", default=str(DATA_DIR),
                        help="Output directory for CSV files")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    ok = run_r_export(args.trials, out_dir)
    if ok:
        check_output(args.trials, out_dir)
    else:
        print("\nR export failed. Possible causes:")
        print("  - R not installed or not on PATH")
        print("  - GitHub install blocked (no internet / firewall)")
        print("  - kmdata package renamed or moved")
        print("\nAlternative: install R package manually and re-run:")
        print("  Rscript -e \"remotes::install_github('raredd/kmdata')\"")
        sys.exit(1)


if __name__ == "__main__":
    main()
