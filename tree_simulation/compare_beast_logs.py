#!/usr/bin/env python3
"""
Combine BEAST 2 results across many simulated trees and compare "kinds" of log files
(i.e. different analysis assumptions).

Layout expected
---------------
    root/
        tree_01/  base_combined.log  sections_combined.log  true_values_combined.log
                  whole_period_combined.log  (other files are ignored)
        tree_02/  ...
        ...

* every sub-directory = one simulated tree (the same tree underlies all logs in it)
* every log file     = one run (one per kind and tree)
* "kind"             = the file name without "_combined"
                       (base_combined.log -> "base"); change with --kind-regex/--pattern.
                       Kinds correspond to the different assumptions you want to test.

How results are combined
------------------------
1. Per run: Tracer-style statistics (median/mean, 95% HPD, ESS, ...).
2. Per tree and kind: (if several runs exist they are averaged) -> one value per tree.
3. Per kind: those per-tree values are combined (mean + bootstrap 95% CI over TREES,
   so tree-to-tree variability is properly accounted for), and the posterior samples
   of all runs/trees are pooled with equal weight per run.
4. Kinds are compared with each other: paired by tree (all kinds share every tree),
   with Wilcoxon / Friedman tests, tree-adjusted effects and, if the true simulated
   value is given (--truth), accuracy and 95% HPD coverage.
5. Efficiency: ESS per million MCMC states (ESS / chain length * 1e6) is computed for every
   run and parameter, so kinds can be compared even when the chains have different lengths.
6. Reliability filter: a (run, parameter) estimate whose ESS is below --min-ess (default 200)
   is treated as unreliable and left out of ALL comparisons of the posteriors (tree-level
   values, summaries, tests, pooled posteriors, accuracy/coverage). The ESS plots
   (ESS_overview.png, ESS_per_million_overview.png) and the ESS columns of
   summary_by_kind.csv deliberately use ALL runs, otherwise kinds that mix badly would look
   better than they are. Dropped entries are listed in excluded_low_ESS_runs.csv.

Requires analyse_beast_logs.py (the first script) in the same folder.

Examples
--------
    python summarise_beast_directories.py /path/to/root -o summary \
        --truth samplingProportion=0.05 --reference true_values
"""
import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats as sps

sys.path.insert(0, str(Path(__file__).resolve().parent))
from analyse_beast_logs import apply_burnin, compute_stats, hpd_interval, read_log, safe_name  # noqa: E402

DEFAULT_PARAMS = ["samplingProportion", "posterior", "likelihood", "prior",
                  "Tree.height", "clockRate", "originBDMMPrime"]
DEFAULT_KIND_REGEX = r"^(.*?)(?:_combined)?$"   # base_combined -> base
DEFAULT_PATTERN = r"^(base|sections|true_values|whole_period)_combined\.log$"
RNG = np.random.default_rng(12345)


# --------------------------------------------------------------------------- #
# Collecting
# --------------------------------------------------------------------------- #
def resolve_column(columns, wanted):
    """Case-insensitive exact match, else first column starting with `wanted`."""
    low = {c.lower(): c for c in columns}
    if wanted.lower() in low:
        return low[wanted.lower()]
    cands = [c for c in columns if c.lower().startswith(wanted.lower())]
    return cands[0] if cands else None


def collect(root, args):
    dir_rx = re.compile(args.dir_pattern)
    file_rx = re.compile(args.pattern)
    kind_rx = re.compile(args.kind_regex)

    def kind_of(stem):
        m = kind_rx.search(stem)
        if m:
            return m.group(1) if m.groups() else m.group(0)
        return stem

    dirs = sorted(d for d in root.iterdir() if d.is_dir() and dir_rx.search(d.name))
    if not dirs:
        sys.exit(f"No sub-directories of {root} match {args.dir_pattern!r}")

    rows = []
    pooled = defaultdict(list)  # (kind, parameter) -> list of sub-sampled arrays (one per run)
    for d in dirs:
        n_files = 0
        for f in sorted(d.iterdir()):
            if not (f.is_file() and file_rx.match(f.name)):
                continue
            try:
                df = read_log(f)
            except Exception as e:
                print(f"  {d.name}/{f.name}: could not read ({e}), skipped")
                continue
            if "Sample" not in df.columns or len(df) < 10:
                print(f"  {d.name}/{f.name}: too short, skipped")
                continue
            kind = kind_of(f.stem)
            step = float(np.median(np.diff(df["Sample"])))
            post, _ = apply_burnin(df, args.burnin)
            # chain length in MCMC states: all logged states, or only those after burn-in
            n_states = (len(df) if args.ess_basis == "total" else len(post)) * step
            n_files += 1
            for label in args.params:
                col = resolve_column([c for c in post.columns if c != "Sample"], label)
                if col is None:
                    continue
                x = post[col].values
                st = compute_stats(x, step)
                rows.append({"directory": d.name, "file": f.name, "kind": kind,
                             "parameter": label, **st,
                             "n_states": n_states,
                             "ESS_per_million": st["ESS"] / (n_states / 1e6)})
                idx = np.linspace(0, len(x) - 1, min(len(x), args.pool_samples)).astype(int)
                # keep run identity so the pooled samples can be filtered like the runs
                pooled[(kind, label)].append((d.name, f.name, x[idx]))
        print(f"  {d.name}: {n_files} run(s)")
    if not rows:
        sys.exit("No estimates collected.")
    return pd.DataFrame(rows), pooled


# --------------------------------------------------------------------------- #
# Combining across trees
# --------------------------------------------------------------------------- #
def tree_level(runs, estimate):
    """One value per (tree, kind, parameter): mean of the run-level estimates."""
    return (runs.groupby(["parameter", "directory", "kind"], as_index=False)
                .agg(value=(estimate, "mean"), n_runs=(estimate, "size")))


def bootstrap_ci(values, n_boot=5000, level=0.95):
    v = np.asarray(values, dtype=float)
    v = v[np.isfinite(v)]
    if len(v) < 2:
        return np.nan, np.nan
    means = RNG.choice(v, size=(n_boot, len(v)), replace=True).mean(axis=1)
    a = (1 - level) / 2
    return np.quantile(means, a), np.quantile(means, 1 - a)


def reliable_mask(runs, min_ess):
    """True for (run, parameter) entries with ESS >= min_ess (NaN ESS counts as unreliable).
    min_ess <= 0 switches the filter off."""
    if min_ess <= 0:
        return pd.Series(True, index=runs.index)
    return runs["ESS"] >= min_ess


def select_pooled(pooled, runs):
    """Keep only the pooled samples that belong to runs still present in `runs`."""
    keep = set(zip(runs["kind"], runs["parameter"], runs["directory"], runs["file"]))
    out = defaultdict(list)
    for (kind, par), items in pooled.items():
        for directory, fname, arr in items:
            if (kind, par, directory, fname) in keep:
                out[(kind, par)].append(arr)
    return out


def summarise(runs_all, tl, pooled, truth, args):
    """`runs_all`: every run. Posterior summaries use only the runs with ESS >= --min-ess
    (`tl` and `pooled` must already be filtered the same way); the ESS columns use all runs."""
    out = []
    for (kind, par), g_all in runs_all.groupby(["kind", "parameter"]):
        g = g_all[reliable_mask(g_all, args.min_ess)]
        row = {
            "kind": kind, "parameter": par,
            "n_runs_total": len(g_all), "n_runs": len(g),
            "n_runs_excluded_low_ESS": len(g_all) - len(g),
            "median_ESS": g_all["ESS"].median(),
            "median_ESS_per_million": g_all["ESS_per_million"].median(),
            "frac_runs_ESS_below_min": float((g_all["ESS"] < args.min_ess).mean()),
        }
        if g.empty:   # every run of this kind is unreliable for this parameter
            out.append(row)
            continue
        t = tl[(tl["kind"] == kind) & (tl["parameter"] == par)]["value"]
        lo, hi = bootstrap_ci(t.values)
        p = np.concatenate(pooled[(kind, par)])
        plo, phi = hpd_interval(p)
        width = g["hpd95_upper"] - g["hpd95_lower"]
        row.update({
            "n_trees": t.size,
            "combined_estimate": t.mean(), "combined_CI95_lower": lo, "combined_CI95_upper": hi,
            "between_tree_sd": t.std(ddof=1),
            "pooled_posterior_median": np.median(p),
            "pooled_posterior_HPD95_lower": plo, "pooled_posterior_HPD95_upper": phi,
            "median_run_HPD_width": width.median(),
        })
        if par in truth:
            tv = truth[par]
            err = g[args.estimate] - tv
            row.update({
                "truth": tv, "bias": err.mean(),
                "relative_bias": err.mean() / tv if tv else np.nan,
                "RMSE": float(np.sqrt((err ** 2).mean())),
                "coverage_95HPD": float(((g["hpd95_lower"] <= tv) & (tv <= g["hpd95_upper"])).mean()),
            })
        out.append(row)
    return pd.DataFrame(out)


def comparison_tests(tl, reference):
    """Tree-paired comparisons: every kind vs reference + Friedman over all kinds."""
    pair_rows, fr_rows = [], []
    for par, g in tl.groupby("parameter"):
        wide = g.pivot(index="directory", columns="kind", values="value")
        complete = wide.dropna()
        if wide.shape[1] >= 3 and len(complete) >= 5:
            try:
                with np.errstate(all="ignore"):
                    s = sps.friedmanchisquare(*[complete[c] for c in complete.columns])
                fr_rows.append({"parameter": par, "n_trees": len(complete),
                                "friedman_chi2": s.statistic, "friedman_p": s.pvalue})
            except ValueError:
                pass
        if reference not in wide.columns:
            continue
        for kind in wide.columns:
            if kind == reference:
                continue
            d = (wide[kind] - wide[reference]).dropna()
            p = np.nan
            if len(d) >= 5 and (d != 0).any():
                try:
                    p = sps.wilcoxon(d).pvalue
                except ValueError:
                    pass
            lo, hi = bootstrap_ci(d.values)
            pair_rows.append({"parameter": par, "kind": kind, "reference": reference,
                              "n_trees": len(d), "mean_difference": d.mean(),
                              "difference_CI95_lower": lo, "difference_CI95_upper": hi,
                              "median_difference": d.median(), "wilcoxon_p": p})
    return pd.DataFrame(pair_rows), pd.DataFrame(fr_rows)


# --------------------------------------------------------------------------- #
# Plots
# --------------------------------------------------------------------------- #
def palette_for(kinds):
    return dict(zip(kinds, sns.color_palette("tab10" if len(kinds) <= 10 else "husl", len(kinds))))


def tree_palette(trees):
    """One fixed colour per simulated tree (directory), identical in every plot."""
    n = len(trees)
    pal = sns.color_palette("tab10" if n <= 10 else "tab20" if n <= 20 else "husl", n)
    return dict(zip(trees, pal))


def strip_by_tree(ax, data, x, y, order, tree_pal, size=5):
    """Dots coloured by the underlying simulated tree (same colour = same tree)."""
    sns.stripplot(data=data, x=x, y=y, order=order, hue="directory",
                  hue_order=list(tree_pal), palette=tree_pal, size=size, alpha=0.9,
                  edgecolor="black", linewidth=0.4, legend=False, ax=ax, zorder=3)


def add_tree_legend(ax, tree_pal):
    """Legend mapping colour -> tree (skipped when there are too many trees)."""
    if len(tree_pal) > 40:
        return
    handles = [Line2D([0], [0], marker="o", ls="", color=c, markeredgecolor="black",
                      markeredgewidth=0.4, label=d) for d, c in tree_pal.items()]
    prev = ax.get_legend()
    if prev is not None:
        ax.add_artist(prev)
    ax.legend(handles=handles, title="simulated tree", fontsize=7, title_fontsize=8,
              bbox_to_anchor=(1.01, 1), loc="upper left", ncol=1 if len(tree_pal) <= 20 else 2)


def _finish(fig, path):
    fig.tight_layout()
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _w(n_kinds):
    return max(5.5, 1.0 * n_kinds + 3)


def plot_combined(tl_g, summ, par, kinds, pal, tree_pal, truth, estimate, path):
    """Combined estimate per kind (mean over trees, bootstrap CI over trees) + tree-level points."""
    fig, ax = plt.subplots(figsize=(_w(len(kinds)), 5))
    for i, kind in enumerate(kinds):
        sub = tl_g[tl_g["kind"] == kind]
        ax.scatter(i + RNG.uniform(-0.18, 0.18, len(sub)), sub["value"].values,
                   color=[tree_pal[d] for d in sub["directory"]], edgecolor="black",
                   linewidth=0.4, alpha=0.9, s=28, zorder=2)
        r = summ[summ["kind"] == kind].iloc[0]
        ax.errorbar(i, r["combined_estimate"],
                    yerr=[[r["combined_estimate"] - r["combined_CI95_lower"]],
                          [r["combined_CI95_upper"] - r["combined_estimate"]]],
                    fmt="D", color="black", ms=7, capsize=5, lw=2, zorder=4,
                    label="mean over trees\n(95% bootstrap CI)" if i == 0 else None)
    if truth is not None:
        ax.axhline(truth, color="red", ls="--", lw=1.2, label=f"true value = {truth:g}")
    ax.set_xticks(range(len(kinds)))
    ax.set_xticklabels(kinds, rotation=45, ha="right")
    ax.set_ylabel(f"{par} (per-tree {estimate})")
    ax.set_title(f"{par}: combined over {tl_g['directory'].nunique()} simulated trees")
    ax.legend(fontsize=8, loc="upper left")
    add_tree_legend(ax, tree_pal)
    _finish(fig, path)


def plot_pooled(pooled, par, kinds, pal, truth, path, max_points=30000):
    rows = []
    for kind in kinds:
        p = np.concatenate(pooled[(kind, par)])
        if len(p) > max_points:
            p = RNG.choice(p, max_points, replace=False)
        rows.append(pd.DataFrame({"kind": kind, "value": p}))
    long = pd.concat(rows)
    fig, ax = plt.subplots(figsize=(_w(len(kinds)), 5))
    sns.violinplot(data=long, x="kind", y="value", order=kinds, palette=pal, hue="kind",
                   legend=False, ax=ax, cut=0, inner=None, linewidth=0.8)
    sns.boxplot(data=long, x="kind", y="value", order=kinds, ax=ax, width=0.12,
                showfliers=False, boxprops={"facecolor": "white", "zorder": 3},
                medianprops={"color": "black"})
    if truth is not None:
        ax.axhline(truth, color="red", ls="--", lw=1.2, label=f"true value = {truth:g}")
        ax.legend(fontsize=8)
    ax.set_xlabel("Kind of log file")
    ax.set_ylabel(par)
    ax.set_title(f"{par}: pooled posterior of all runs and trees")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    _finish(fig, path)


def plot_runs(g, par, kinds, pal, tree_pal, truth, estimate, path):
    fig, ax = plt.subplots(figsize=(_w(len(kinds)), 5))
    sns.boxplot(data=g, x="kind", y=estimate, order=kinds, palette=pal, hue="kind",
                legend=False, ax=ax, showfliers=False, boxprops={"alpha": 0.5})
    strip_by_tree(ax, g, "kind", estimate, kinds, tree_pal)
    if truth is not None:
        ax.axhline(truth, color="red", ls="--", lw=1.2, label=f"true value = {truth:g}")
        ax.legend(fontsize=8, loc="upper left")
    ax.set_xlabel("Kind of log file")
    ax.set_ylabel(f"{par} ({estimate} of each run)")
    ax.set_title(f"{par}: every individual run (n = {len(g)}), coloured by simulated tree")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    add_tree_legend(ax, tree_pal)
    _finish(fig, path)


def plot_tree_adjusted(tl_g, par, kinds, pal, tree_pal, path):
    """Per-tree value minus the tree's mean over kinds: removes the tree effect."""
    wide = tl_g.pivot(index="directory", columns="kind", values="value")[kinds].dropna()
    if len(wide) < 2 or len(kinds) < 2:
        return
    adj = wide.sub(wide.mean(axis=1), axis=0).reset_index().melt(
        id_vars="directory", var_name="kind", value_name="adj")
    fig, ax = plt.subplots(figsize=(_w(len(kinds)), 5))
    sns.boxplot(data=adj, x="kind", y="adj", order=kinds, palette=pal, hue="kind",
                legend=False, ax=ax, showfliers=False, boxprops={"alpha": 0.5})
    strip_by_tree(ax, adj, "kind", "adj", kinds, tree_pal)
    ax.axhline(0, color="red", ls="--", lw=1)
    ax.set_xlabel("Kind of log file")
    ax.set_ylabel(f"{par}: estimate − mean of all kinds\n(same tree)")
    ax.set_title(f"{par}: tree-adjusted effect of each assumption (n = {len(wide)} trees)")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    add_tree_legend(ax, tree_pal)
    _finish(fig, path)


def plot_differences(tl_g, par, kinds, reference, pairs, tree_pal, path):
    wide = tl_g.pivot(index="directory", columns="kind", values="value")
    others = [k for k in kinds if k != reference]
    if reference not in wide.columns or not others:
        return
    parts = []
    for k in others:
        if k in wide.columns:
            d = (wide[k] - wide[reference]).dropna()
            parts.append(pd.DataFrame({"kind": k, "diff": d.values, "directory": d.index}))
    long = pd.concat(parts)
    fig, ax = plt.subplots(figsize=(_w(len(others)), 5))
    sns.boxplot(data=long, x="kind", y="diff", order=others, ax=ax, color="lightsteelblue",
                showfliers=False)
    strip_by_tree(ax, long, "kind", "diff", others, tree_pal)
    ax.axhline(0, color="red", ls="--", lw=1)
    t = pairs[pairs["parameter"] == par].set_index("kind") if len(pairs) else None
    if t is not None:
        ymax = ax.get_ylim()[1]
        for i, k in enumerate(others):
            if k in t.index and np.isfinite(t.loc[k, "wilcoxon_p"]):
                ax.text(i, ymax, f"p = {t.loc[k, 'wilcoxon_p']:.2g}", ha="center",
                        va="bottom", fontsize=8)
    ax.set_xlabel("Kind of log file")
    ax.set_ylabel(f"{par}: kind − {reference} (per tree)")
    ax.set_title(f"{par}: paired difference to '{reference}'\n(Wilcoxon signed-rank p on top)",
                 pad=18)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    add_tree_legend(ax, tree_pal)
    _finish(fig, path)


def plot_truth_metrics(g, par, kinds, pal, tree_pal, truth, estimate, path):
    g = g.copy()
    g["rel_error"] = (g[estimate] - truth) / truth
    g["covered"] = (g["hpd95_lower"] <= truth) & (truth <= g["hpd95_upper"])
    fig, axes = plt.subplots(1, 2, figsize=(max(9, 1.4 * len(kinds) + 6), 4.8))
    sns.boxplot(data=g, x="kind", y="rel_error", order=kinds, palette=pal, hue="kind",
                legend=False, ax=axes[0], showfliers=False, boxprops={"alpha": 0.5})
    strip_by_tree(axes[0], g, "kind", "rel_error", kinds, tree_pal)
    axes[0].axhline(0, color="red", ls="--", lw=1)
    axes[0].set_ylabel(f"relative error of {par}\n(estimate − truth) / truth")
    axes[0].set_title("Accuracy (all runs)")
    cov = g.groupby("kind")["covered"].agg(["mean", "count"]).reindex(kinds)
    axes[1].bar(range(len(kinds)), cov["mean"], color=[pal[k] for k in kinds])
    axes[1].axhline(0.95, color="red", ls="--", lw=1, label="nominal 95%")
    for i, (m, n) in enumerate(zip(cov["mean"], cov["count"])):
        axes[1].text(i, m + 0.01, f"{m:.0%}\n(n={n})", ha="center", fontsize=8)
    axes[1].set_xticks(range(len(kinds)))
    axes[1].set_xticklabels(kinds)
    axes[1].set_ylim(0, 1.15)
    axes[1].set_ylabel("fraction of runs with truth in 95% HPD")
    axes[1].set_title("Coverage")
    axes[1].legend(fontsize=8, loc="lower right")
    for a in axes:
        plt.setp(a.get_xticklabels(), rotation=45, ha="right")
        a.set_xlabel("Kind of log file")
    add_tree_legend(axes[1], tree_pal)
    _finish(fig, path)


def plot_hpd_width(g, par, kinds, pal, tree_pal, path):
    g = g.assign(hpd_width=g["hpd95_upper"] - g["hpd95_lower"])
    fig, ax = plt.subplots(figsize=(_w(len(kinds)), 5))
    sns.boxplot(data=g, x="kind", y="hpd_width", order=kinds, palette=pal, hue="kind",
                legend=False, ax=ax, showfliers=False, boxprops={"alpha": 0.5})
    strip_by_tree(ax, g, "kind", "hpd_width", kinds, tree_pal)
    ax.set_ylabel(f"95% HPD width of {par} (per run)")
    ax.set_xlabel("Kind of log file")
    ax.set_title(f"{par}: posterior uncertainty")
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
    add_tree_legend(ax, tree_pal)
    _finish(fig, path)


def plot_ess(runs, kinds, pal, path, column="ESS", ylabel="ESS",
             title="Effective sample sizes (all runs)", threshold=None):
    """Boxplot of an ESS-type column per parameter and kind (log scale).
    `threshold` draws a red dashed line (used for the ESS = 200 cut-off)."""
    runs = runs[np.isfinite(runs[column]) & (runs[column] > 0)]   # log scale needs > 0
    fig, ax = plt.subplots(figsize=(max(8, 1.3 * runs["parameter"].nunique() + 3), 5))
    sns.boxplot(data=runs, x="parameter", y=column, hue="kind", hue_order=kinds, palette=pal,
                ax=ax, showfliers=False)
    ax.set_yscale("log")
    if threshold:
        ax.axhline(threshold, color="red", ls="--", lw=1, label=f"ESS = {threshold:g}")
    ax.set_ylabel(ylabel)
    ax.legend(fontsize=8, bbox_to_anchor=(1.01, 1), loc="upper left")
    ax.set_title(title)
    plt.setp(ax.get_xticklabels(), rotation=30, ha="right")
    _finish(fig, path)


# --------------------------------------------------------------------------- #
# Main
# --------------------------------------------------------------------------- #
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("root", type=Path, help="folder with one sub-directory per simulated tree")
    ap.add_argument("-o", "--output", type=Path, default=Path("beast_summary"))
    ap.add_argument("--dir-pattern", default=r".*",
                    help="regex for sub-directory names to include (default: all)")
    ap.add_argument("--pattern", default=DEFAULT_PATTERN,
                    help="regex for log file NAMES to include (default: the four *_combined.log "
                         "kinds); all other files in the directories are ignored")
    ap.add_argument("--kind-regex", default=DEFAULT_KIND_REGEX,
                    help="regex applied to the file stem; first group = kind name "
                         f"(default {DEFAULT_KIND_REGEX!r}: strips '_combined'). "
                         "Use '^(.*)$' to make every file its own kind.")
    ap.add_argument("--params", nargs="*", default=DEFAULT_PARAMS,
                    help="parameters to summarise (case-insensitive, prefix match); "
                         "the first one is the focus parameter")
    ap.add_argument("--truth", nargs="*", default=[], metavar="PARAM=VALUE",
                    help="true simulated values, e.g. samplingProportion=0.05")
    ap.add_argument("--estimate", choices=["median", "mean"], default="median",
                    help="point estimate per run (default: median)")
    ap.add_argument("--reference", default=None,
                    help="kind used as reference for paired differences (default: first kind)")
    ap.add_argument("--burnin", type=float, default=0.1)
    ap.add_argument("--pool-samples", type=int, default=500,
                    help="samples taken (evenly spaced) from each run for the pooled posterior")
    ap.add_argument("--min-ess", type=float, default=200.0,
                    help="(run, parameter) estimates with ESS below this are unreliable and are "
                         "left out of the comparisons (default 200; 0 = keep everything). "
                         "ESS plots always show all runs; dropped entries are listed in "
                         "excluded_low_ESS_runs.csv")
    ap.add_argument("--ess-basis", choices=["total", "postburnin"], default="total",
                    help="chain length used for ESS per million states: all logged states "
                         "('total', default; burn-in counts as computational cost) or only the "
                         "states kept after burn-in ('postburnin')")
    args = ap.parse_args()

    truth = {}
    for item in args.truth:
        k, _, v = item.partition("=")
        truth[k] = float(v)

    sns.set_theme(style="whitegrid")
    args.output.mkdir(parents=True, exist_ok=True)

    print("Reading logs...")
    runs_all, pooled = collect(args.root, args)
    runs_all.to_csv(args.output / "all_runs.csv", index=False)

    # ---- reliability filter: ESS < --min-ess => (run, parameter) is left out of comparisons
    ok = reliable_mask(runs_all, args.min_ess)
    runs = runs_all[ok]
    dropped = runs_all[~ok]
    dropped[["directory", "file", "kind", "parameter", "ESS", "ESS_per_million"]].to_csv(
        args.output / "excluded_low_ESS_runs.csv", index=False)
    if len(dropped):
        print(f"\nExcluded {len(dropped)} of {len(runs_all)} (run, parameter) estimates with "
              f"ESS < {args.min_ess:g} from the comparisons:")
        print(dropped.groupby(["parameter", "kind"]).size().rename("n_excluded")
              .to_string())
        print(f"  (full list: {args.output / 'excluded_low_ESS_runs.csv'})\n")
    else:
        print(f"All {len(runs_all)} estimates have ESS >= {args.min_ess:g}.")
    if runs.empty:
        sys.exit(f"Every estimate has ESS < {args.min_ess:g}; nothing left to compare. "
                 f"Run longer chains or lower --min-ess.")
    pooled_kept = select_pooled(pooled, runs)

    # kinds / colours come from ALL runs so they are identical in every plot,
    # even if a kind is completely excluded for some parameter
    kinds = sorted(runs_all["kind"].unique())
    reference = args.reference or kinds[0]
    if reference not in kinds:
        sys.exit(f"Reference kind {reference!r} not found. Kinds: {kinds}")
    print(f"Kinds found: {kinds}  (reference: {reference})")
    pal = palette_for(kinds)
    tree_pal = tree_palette(sorted(runs_all['directory'].unique()))

    tl = tree_level(runs, args.estimate)
    tl.to_csv(args.output / "tree_level_estimates.csv", index=False)
    summary = summarise(runs_all, tl, pooled_kept, truth, args)
    summary.to_csv(args.output / "summary_by_kind.csv", index=False)
    pairs, friedman = comparison_tests(tl, reference)
    pairs.to_csv(args.output / "paired_differences_vs_reference.csv", index=False)
    friedman.to_csv(args.output / "friedman_all_kinds.csv", index=False)
    (runs_all.groupby(["directory", "kind"]).size().unstack(fill_value=0)
     .to_csv(args.output / "runs_per_directory_and_kind.csv"))

    for i, par in enumerate(args.params):
        g = runs[runs["parameter"] == par]
        if g.empty:
            why = ("not found in any file" if par not in set(runs_all["parameter"])
                   else f"all runs have ESS < {args.min_ess:g}")
            print(f"  '{par}': {why}, skipped")
            continue
        sub = args.output / ("00_FOCUS_" + safe_name(par) if i == 0 else safe_name(par))
        sub.mkdir(exist_ok=True)
        t = truth.get(par)
        tl_g = tl[tl["parameter"] == par]
        summ_g = summary[summary["parameter"] == par]
        pk = [k for k in kinds if k in set(g["kind"])]
        plot_combined(tl_g, summ_g, par, pk, pal, tree_pal, t, args.estimate, sub / "combined_over_trees.png")
        plot_pooled(pooled_kept, par, pk, pal, t, sub / "pooled_posterior.png")
        plot_runs(g, par, pk, pal, tree_pal, t, args.estimate, sub / "individual_runs.png")
        plot_hpd_width(g, par, pk, pal, tree_pal, sub / "hpd_width.png")
        if len(pk) > 1:
            plot_tree_adjusted(tl_g, par, pk, pal, tree_pal, sub / "tree_adjusted_effects.png")
            plot_differences(tl_g, par, pk, reference, pairs, tree_pal, sub / "paired_difference_vs_reference.png")
        if t is not None:
            plot_truth_metrics(g, par, pk, pal, tree_pal, t, args.estimate, sub / "accuracy_and_coverage.png")

    # ESS plots use ALL runs on purpose (see module docstring)
    plot_ess(runs_all, kinds, pal, args.output / "ESS_overview.png",
             threshold=args.min_ess if args.min_ess > 0 else 200)
    plot_ess(runs_all, kinds, pal, args.output / "ESS_per_million_overview.png",
             column="ESS_per_million", ylabel="ESS per million MCMC states",
             title="Sampling efficiency: ESS per million states (all runs)")

    focus = args.params[0]
    pd.set_option("display.width", 220)
    pd.set_option("display.max_columns", 30)
    fs = summary[summary["parameter"] == focus]
    if not fs.empty:
        print(f"\nFocus parameter '{focus}':")
        print(fs.drop(columns="parameter").round(5).to_string(index=False))
    fp = pairs[pairs["parameter"] == focus]
    if not fp.empty:
        print(f"\nTree-paired differences vs '{reference}':")
        print(fp.drop(columns=["parameter", "reference"]).round(5).to_string(index=False))
    print(f"\nDone. Results in {args.output.resolve()}")


if __name__ == "__main__":
    main()