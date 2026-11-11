#!/usr/bin/env python3
"""
Analyse BEAST 2 log files: Tracer-style statistics, trace plots,
histograms with 95% HPD, and box/violin comparisons between runs.

Example
-------
    python analyse_beast_logs.py /path/to/logs -o results
    python analyse_beast_logs.py /path/to/logs --pattern "^.+_combined\.log$" --burnin 0.1
"""
import argparse
import re
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

# Parameters compared across files (matched case-insensitively, so the
# "samplingProportionSPEPi" spelling also finds "samplingProportionSPEpi").
COMPARE_PARAMS = [
    "posterior", "likelihood", "prior", "Tree.height", "clockRate",
    "originBDMMPrime", "samplingProportionSPEPi",
]

MAX_LAG = 2000  # same cap Tracer uses for the autocorrelation sum


# --------------------------------------------------------------------------- #
# Reading
# --------------------------------------------------------------------------- #
def read_log(path):
    """Read a BEAST log (lines starting with '#' are comments). Returns a DataFrame."""
    df = pd.read_csv(path, sep="\t", comment="#", skip_blank_lines=True)
    df.columns = [c.strip() for c in df.columns]
    # keep only numeric columns
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.dropna(axis=1, how="all").dropna(axis=0, how="any")
    return df


def apply_burnin(df, burnin_frac):
    n_burn = int(len(df) * burnin_frac)
    return df.iloc[n_burn:].reset_index(drop=True), n_burn


# --------------------------------------------------------------------------- #
# Statistics (Tracer definitions)
# --------------------------------------------------------------------------- #
def hpd_interval(x, mass=0.95):
    """Shortest interval containing `mass` of the samples (as in Tracer)."""
    s = np.sort(np.asarray(x))
    n = len(s)
    k = int(np.floor(mass * n))
    if k < 1 or n < 2:
        return np.nan, np.nan
    widths = s[k:] - s[: n - k]
    i = int(np.argmin(widths))
    return s[i], s[i + k]


def autocorr_time(x, step_size):
    """
    Integrated autocorrelation time and variance estimate, following Tracer's
    IntegerAutocorrelation approach (Geyer-style pairing of lags).
    Returns (ACT in states, var_stat).
    """
    x = np.asarray(x, dtype=float)
    n = len(x)
    xc = x - x.mean()
    max_lag = min(n - 1, MAX_LAG)
    if n < 3 or np.all(xc == 0):
        return np.nan, np.nan

    # autocovariance via FFT, normalised by (n - lag) like Tracer
    nfft = 1 << (2 * n - 1).bit_length()
    f = np.fft.rfft(xc, nfft)
    acov = np.fft.irfft(f * np.conj(f), nfft)[: max_lag + 1]
    gamma = acov / (n - np.arange(max_lag + 1))

    var_stat = gamma[0]
    for lag in range(1, max_lag + 1):
        if lag % 2 == 0:
            pair = gamma[lag - 1] + gamma[lag]
            if pair > 0:
                var_stat += 2.0 * pair
            else:
                break
    act = step_size * var_stat / gamma[0]
    return act, var_stat


def compute_stats(x, step_size):
    x = np.asarray(x, dtype=float)
    n = len(x)
    mean = x.mean()
    var = x.var(ddof=1) if n > 1 else np.nan
    act, var_stat = autocorr_time(x, step_size)
    ess = (step_size * n) / act if act and np.isfinite(act) and act > 0 else np.nan
    stderr = np.sqrt(var_stat / n) if np.isfinite(var_stat) else np.nan
    lo, hi = hpd_interval(x)
    geo = np.exp(np.mean(np.log(x))) if np.all(x > 0) else np.nan  # Tracer: positive data only
    return {
        "mean": mean,
        "stderr_of_mean": stderr,
        "stdev": np.sqrt(var) if np.isfinite(var) else np.nan,
        "variance": var,
        "median": np.median(x),
        "min": x.min(),
        "max": x.max(),
        "range": x.max() - x.min(),
        "geometric_mean": geo,
        "hpd95_lower": lo,
        "hpd95_upper": hi,
        "ACT": act,
        "ESS": ess,
        "n_samples": n,
    }


# --------------------------------------------------------------------------- #
# Plotting
# --------------------------------------------------------------------------- #
def safe_name(s):
    return re.sub(r"[^A-Za-z0-9._-]+", "_", s)


def plot_trace(full_x, full_samples, n_burn, name, run, stats, path):
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(full_samples, full_x, lw=0.6, color="tab:blue")
    if n_burn > 0:
        ax.axvspan(full_samples[0], full_samples[n_burn - 1], color="grey", alpha=0.25,
                   label=f"burn-in ({n_burn} samples)")
    ax.axhline(stats["mean"], color="red", ls="--", lw=1, label=f"mean = {stats['mean']:.4g}")
    ax.axhspan(stats["hpd95_lower"], stats["hpd95_upper"], color="red", alpha=0.08,
               label="95% HPD")
    ax.set_xlabel("State")
    ax.set_ylabel(name)
    ax.set_title(f"{run}: trace of {name} (ESS = {stats['ESS']:.0f})")
    ax.legend(loc="best", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_hist(x, name, run, stats, path):
    fig, ax = plt.subplots(figsize=(6, 4))
    sns.histplot(x, bins=40, kde=len(np.unique(x)) > 1, ax=ax, color="tab:blue")
    ax.axvline(stats["mean"], color="red", ls="--", lw=1.2, label=f"mean = {stats['mean']:.4g}")
    ax.axvline(stats["hpd95_lower"], color="black", ls=":", lw=1.2,
               label=f"95% HPD [{stats['hpd95_lower']:.4g}, {stats['hpd95_upper']:.4g}]")
    ax.axvline(stats["hpd95_upper"], color="black", ls=":", lw=1.2)
    ax.axvspan(stats["hpd95_lower"], stats["hpd95_upper"], color="black", alpha=0.06)
    ax.set_xlabel(name)
    ax.set_title(f"{run}: {name}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def plot_comparison(long_df, param, kind, path):
    n_runs = long_df["run"].nunique()
    fig, ax = plt.subplots(figsize=(max(6, 0.6 * n_runs + 3), 5))
    if kind == "box":
        sns.boxplot(data=long_df, x="run", y="value", ax=ax, color="lightsteelblue",
                    fliersize=1.5)
    else:
        sns.violinplot(data=long_df, x="run", y="value", ax=ax, color="lightsteelblue",
                       inner="quartile", cut=0)
    ax.set_xlabel("Run")
    ax.set_ylabel(param)
    ax.set_title(f"{param} ({'boxplot' if kind == 'box' else 'violin plot'})")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def find_column(columns, wanted):
    """Case-insensitive column lookup."""
    lookup = {c.lower(): c for c in columns}
    return lookup.get(wanted.lower())


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("folder", type=Path, help="folder containing the log files")
    ap.add_argument("-o", "--output", type=Path, default=Path("data/results/posterior_analysis"),
                    help="output folder (default: beast_analysis)")
    ap.add_argument("--pattern", default=r"^.+_combined\.log$",
                    help=r"regex the file NAME must match (default: ^.+_combined\.log$)")
    ap.add_argument("--burnin", type=float, default=0.1,
                    help="burn-in fraction to discard (default 0.1, like Tracer)")
    ap.add_argument("--compare", nargs="*", default=COMPARE_PARAMS,
                    help="parameters for box/violin comparison")
    args = ap.parse_args()

    rx = re.compile(args.pattern)
    files = sorted(p for p in args.folder.iterdir() if p.is_file() and rx.match(p.name))
    if not files:
        sys.exit(f"No files in {args.folder} match pattern {args.pattern!r}")
    print(f"Analysing {len(files)} file(s)")

    out = args.output
    (out / "stats").mkdir(parents=True, exist_ok=True)
    (out / "comparison").mkdir(exist_ok=True)
    sns.set_theme(style="whitegrid")

    all_stats = []
    samples_by_run = {}  # run -> post-burn-in DataFrame

    for f in files:
        run = f.stem
        print(f"  {f.name}")
        df = read_log(f)
        if "Sample" not in df.columns or len(df) < 10:
            print(f"    skipped: no 'Sample' column or too few rows")
            continue
        step = float(np.median(np.diff(df["Sample"]))) if len(df) > 1 else 1.0
        post, n_burn = apply_burnin(df, args.burnin)
        samples_by_run[run] = post

        trace_dir = out / "traces" / run
        hist_dir = out / "histograms" / run
        trace_dir.mkdir(parents=True, exist_ok=True)
        hist_dir.mkdir(parents=True, exist_ok=True)

        for col in df.columns:
            if col == "Sample":
                continue
            st = compute_stats(post[col].values, step)
            all_stats.append({"file": f.name, "parameter": col, **st})
            plot_trace(df[col].values, df["Sample"].values, n_burn, col, run, st,
                       trace_dir / f"{safe_name(col)}.png")
            plot_hist(post[col].values, col, run, st, hist_dir / f"{safe_name(col)}.png")

    stats_df = pd.DataFrame(all_stats)
    stats_df.to_csv(out / "stats" / "all_statistics.csv", index=False)
    # one CSV per file as well
    for fname, sub in stats_df.groupby("file"):
        sub.drop(columns="file").to_csv(out / "stats" / f"{Path(fname).stem}_statistics.csv",
                                        index=False)

    # ---- comparison plots ---------------------------------------------------
    for wanted in args.compare:
        rows = []
        actual = None
        for run, post in samples_by_run.items():
            col = find_column(post.columns, wanted)
            if col is None:
                print(f"  note: '{wanted}' not in {run}; run left out of that comparison")
                continue
            actual = col
            rows.append(pd.DataFrame({"run": run, "value": post[col].values}))
        if not rows:
            print(f"  '{wanted}' not found in any file, skipping")
            continue
        long_df = pd.concat(rows, ignore_index=True)
        for kind in ("box", "violin"):
            plot_comparison(long_df, actual, kind,
                            out / "comparison" / f"{safe_name(actual)}_{kind}.png")

    # focus summary for the key parameter
    sp = stats_df[stats_df["parameter"].str.lower().str.startswith("samplingproportion")]
    if not sp.empty:
        sp.to_csv(out / "stats" / "samplingProportion_summary.csv", index=False)

    print(f"Done. Results in {out.resolve()}")


if __name__ == "__main__":
    main()