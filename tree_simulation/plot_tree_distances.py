#!/usr/bin/env python3
"""
Summarize and plot distance matrices from many replicates, using seaborn.

Expected layout (each tree_<k> folder holds the CSVs written by
`tree_distances.py ... --csv PREFIX` for the same set of trees, e.g. 5 tree types):

  input_dir/
    tree_0/  PREFIX_rf.csv  PREFIX_mrca_rmse.csv  ...  PREFIX_trees.csv
    tree_1/  ...
    tree_n/  ...

Trees are matched across replicates by name, so the same tree names must be present in
every folder. For each metric the script computes the mean and the standard deviation of
every pairwise distance across the replicates, and writes to --output_dir:

  heatmaps_mean.png          mean distance per metric, cells annotated "mean ± SD"
  heatmaps_sd.png            standard deviation across replicates per metric
  pairwise_distributions.png box + points of every pair's distance across replicates
  mds_mean.png               2D map (classical MDS) of the mean distances, coloured by tree
  clustered_<m>_mean.png     clustered heatmap of the mean distances (--cluster), with
                             side bars coloured by tree
Colours are consistent across plots: a tree always has the same colour, and a tree pair
has the same colour in every panel of pairwise_distributions.png.
  summary/<metric>_mean.csv, <metric>_sd.csv, pairwise_distances_long.csv

Usage
  python plot_distances.py --input_dir results --output_dir figs
  python plot_distances.py --input_dir results --metrics rf_norm mrca_rmse --cluster rf_norm
  python plot_distances.py --input_dir results --prefix out     # if prefix auto-detect fails

Requires: pip install numpy pandas matplotlib seaborn scipy
"""
import argparse
import itertools
import re
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Patch
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

DISTANCES = ["rf", "rf_norm", "wrf", "bsd", "triplet",
             "mrca_rmse", "mrca_mad", "clade_age_rmse", "clade_age_mad"]
SIMILARITIES = ["n_shared"]  # larger = more similar, so excluded from MDS/clustering


# ------------------------------------------------------------------ loading
def find_replicates(input_dir):
    """tree_0, tree_1, ... sub-folders of input_dir, in numeric order."""
    dirs = [d for d in Path(input_dir).iterdir()
            if d.is_dir() and re.fullmatch(r"tree_\d+", d.name)]
    if len(dirs) < 2:
        raise SystemExit(f"Need at least two tree_<k> folders in {input_dir} "
                         f"(found {len(dirs)}).")
    return sorted(dirs, key=lambda d: int(d.name.split("_")[1]))


def find_prefix(folder, prefix):
    """CSV prefix: given by --prefix, else taken from the single *_trees.csv file."""
    if prefix is not None:
        return prefix
    found = [f.name[: -len("_trees.csv")] for f in folder.glob("*_trees.csv")]
    if len(found) != 1:
        raise SystemExit(f"{folder}: cannot tell the CSV prefix ({len(found)} "
                         f"'*_trees.csv' files found). Pass --prefix.")
    return found[0]


def load_all(dirs, prefix, wanted):
    """Returns (tree names, {metric: array of shape (replicates, n, n)})."""
    metrics = wanted or DISTANCES + SIMILARITIES
    names, data = None, {m: [] for m in metrics}
    for d in dirs:
        p = find_prefix(d, prefix)
        for m in metrics:
            f = d / f"{p}_{m}.csv"
            if not f.exists():
                continue
            df = pd.read_csv(f, index_col=0)
            if names is None:
                names = list(df.index)
            if set(df.index) != set(names) or set(df.columns) != set(names):
                raise SystemExit(f"{f}: tree names differ from those in {dirs[0].name} "
                                 f"({sorted(set(df.index) ^ set(names))}).")
            data[m].append(df.loc[names, names].to_numpy(dtype=float))
    out = {}
    for m, mats in data.items():
        if len(mats) == len(dirs):
            out[m] = np.stack(mats)
        elif mats:
            print(f"Skipping {m}: found in only {len(mats)} of {len(dirs)} folders.")
    if not out:
        raise SystemExit("No metric CSVs found in the tree_<k> folders.")
    return names, out


def summarise(arr):
    """Mean and sample standard deviation across replicates (NaNs ignored)."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        return np.nanmean(arr, axis=0), np.nanstd(arr, axis=0, ddof=1)


def long_table(arrs, names, dirs):
    rows = []
    for m, arr in arrs.items():
        for r, d in enumerate(dirs):
            for i, j in itertools.combinations(range(len(names)), 2):
                rows.append((d.name, names[i], names[j], m, arr[r, i, j]))
    return pd.DataFrame(rows, columns=["replicate", "tree_a", "tree_b", "metric", "value"])


# ------------------------------------------------------------------ helpers
def fill_nan(D):
    """Replace NaN (e.g. no shared clades) by the largest finite distance."""
    D = D.copy()
    D[np.isnan(D)] = np.nanmax(D) if np.isfinite(D).any() else 0.0
    np.fill_diagonal(D, 0.0)
    return D


def classical_mds(D, k=2):
    """Classical (Torgerson) MDS; returns coordinates and variance explained."""
    n = len(D)
    J = np.eye(n) - np.ones((n, n)) / n
    vals, vecs = np.linalg.eigh(-0.5 * J @ (D ** 2) @ J)
    top = np.argsort(vals)[::-1][:k]
    pos_total = np.clip(vals, 0, None).sum()
    coords = vecs[:, top] * np.sqrt(np.clip(vals[top], 0, None))
    return coords, np.clip(vals[top], 0, None) / (pos_total or 1.0)


def palette(keys):
    """One fixed colour per key, reused in every plot."""
    pal = sns.color_palette("tab10" if len(keys) <= 10 else "husl", len(keys))
    return dict(zip(keys, pal))


def legend_outside(fig, colours, title):
    """Legend placed to the right of the figure so it never covers a plot."""
    fig.legend(handles=[Patch(color=c, label=k) for k, c in colours.items()],
               title=title, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False)


def grid(n_plots, ncols, w, h):
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * w, nrows * h), squeeze=False)
    axes = axes.ravel()
    for ax in axes[n_plots:]:
        ax.axis("off")
    return fig, axes


# -------------------------------------------------------------------- plots
def plot_heatmaps(means, sds, names, n_rep, out):
    n = len(names)
    eye = np.eye(n, dtype=bool)
    size = min(max(4.2, 0.9 * n), 9)
    small = n <= 8  # annotate cells only when they stay readable

    # mean, annotated with "mean ± SD"
    fig, axes = grid(len(means), 3, size, size)
    for ax, m in zip(axes, means):
        mean, sd = means[m], sds[m]
        ann = np.array([[f"{a:.2g}\n±{b:.2g}" if np.isfinite(a) else "" for a, b in zip(ra, rb)]
                        for ra, rb in zip(mean, sd)])
        ann[eye] = ""
        df = pd.DataFrame(mean, index=names, columns=names)
        if m in SIMILARITIES:
            df = df.mask(eye)
        sns.heatmap(df, ax=ax, square=True, cmap="mako_r" if m in SIMILARITIES else "mako",
                    annot=ann if small else False, fmt="", annot_kws={"size": 6},
                    cbar_kws={"shrink": 0.8})
        ax.set_title(m)
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Mean ± SD over {n_rep} replicates")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / "heatmaps_mean.png", dpi=200)

    # standard deviation alone
    fig, axes = grid(len(sds), 3, size, size)
    for ax, m in zip(axes, sds):
        df = pd.DataFrame(sds[m], index=names, columns=names).mask(eye)
        sns.heatmap(df, ax=ax, square=True, cmap="flare", annot=small, fmt=".2g",
                    annot_kws={"size": 6}, cbar_kws={"shrink": 0.8})
        ax.set_title(f"{m}: SD")
        ax.tick_params(labelsize=7)
    fig.suptitle(f"Standard deviation over {n_rep} replicates")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / "heatmaps_sd.png", dpi=200)


def plot_distributions(long, names, out):
    pairs = [f"{names[i]} – {names[j]}"
             for i, j in itertools.combinations(range(len(names)), 2)]
    pair_colours = palette(pairs)  # same pair -> same colour in every panel
    long = long.assign(pair=long.tree_a + " – " + long.tree_b)
    metrics = list(dict.fromkeys(long.metric))
    ncols = 2 if len(pairs) <= 15 else 1
    fig, axes = grid(len(metrics), ncols, min(max(6, 0.55 * len(pairs)), 14), 3.8)
    for ax, m in zip(axes, metrics):
        sub = long[long.metric == m]
        sns.boxplot(data=sub, x="pair", y="value", hue="pair", order=pairs,
                    hue_order=pairs, palette=pair_colours, legend=False,
                    fliersize=0, linewidth=0.8, showmeans=True, ax=ax,
                    meanprops=dict(marker="D", markerfacecolor="white",
                                   markeredgecolor="black", markersize=4))
        sns.stripplot(data=sub, x="pair", y="value", order=pairs, size=3, alpha=0.7,
                      color="0.25", ax=ax)
        ax.set_title(m)
        ax.set_xlabel("")
        ax.set_ylabel(m)
        ax.tick_params(axis="x", labelrotation=45, labelsize=7)
        plt.setp(ax.get_xticklabels(), ha="right")
    fig.suptitle("Distance of each tree pair across replicates (diamond = mean)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(out / "pairwise_distributions.png", dpi=200, bbox_inches="tight")


def plot_mds(dist_means, names, colours, n_rep, out):
    fig, axes = grid(len(dist_means), 3, 4.2, 4.2)
    for ax, (m, D) in zip(axes, dist_means.items()):
        xy, var = classical_mds(fill_nan(D))
        sns.scatterplot(x=xy[:, 0], y=xy[:, 1], hue=names, hue_order=names,
                        palette=colours, s=70, legend=False, ax=ax)
        for (x, y), lab in zip(xy, names):
            ax.annotate(lab, (x, y), xytext=(4, 4), textcoords="offset points", fontsize=7)
        ax.set_title(m)
        ax.set_xlabel(f"MDS 1 ({var[0]:.0%})")
        ax.set_ylabel(f"MDS 2 ({var[1]:.0%})")
    fig.suptitle(f"MDS of mean distances ({n_rep} replicates)")
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    legend_outside(fig, colours, "Tree")
    fig.savefig(out / "mds_mean.png", dpi=200, bbox_inches="tight")


def plot_clustermap(D, names, colours, metric, out):
    n = len(names)
    D = fill_nan(D)
    Z = linkage(squareform(D, checks=False), method="average")
    side = pd.Series([colours[k] for k in names], index=names, name="")
    g = sns.clustermap(pd.DataFrame(D, index=names, columns=names),
                       row_linkage=Z, col_linkage=Z, cmap="mako",
                       row_colors=side, col_colors=side,
                       annot=n <= 10, fmt=".2g", annot_kws={"size": 6},
                       figsize=(min(max(6, 0.4 * n + 3), 14),) * 2,
                       cbar_kws={"label": f"mean {metric}"})
    g.ax_heatmap.tick_params(labelsize=7)
    g.ax_row_colors.set_xticks([])  # drop the empty tick marks next to the colour bars
    g.ax_col_colors.set_yticks([])
    g.fig.suptitle(f"mean {metric} (average-linkage clustering)", y=1.01)
    legend_outside(g.fig, colours, "Tree")
    g.savefig(out / f"clustered_{metric}_mean.png", dpi=200, bbox_inches="tight")


# --------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input_dir", required=True,
                    help="folder containing tree_0, tree_1, ... sub-folders with the CSVs")
    ap.add_argument("--output_dir", default=".", help="where to write figures and summaries")
    ap.add_argument("--prefix", help="CSV prefix used with tree_distances.py --csv "
                                     "(default: detected from *_trees.csv)")
    ap.add_argument("--metrics", nargs="+", help="subset of metrics to plot")
    ap.add_argument("--cluster", default="mrca_rmse", help="metric for the clustered heatmap")
    ap.add_argument("--show", action="store_true", help="also open the figures")
    args = ap.parse_args()

    sns.set_theme(style="white", context="paper")
    dirs = find_replicates(args.input_dir)
    names, arrs = load_all(dirs, args.prefix, args.metrics)
    n, n_rep = len(names), len(dirs)
    print(f"{n_rep} replicates, {n} trees: {', '.join(names)}")

    out = Path(args.output_dir)
    (out / "summary").mkdir(parents=True, exist_ok=True)

    means, sds = {}, {}
    for m, arr in arrs.items():
        means[m], sds[m] = summarise(arr)
        pd.DataFrame(means[m], index=names, columns=names).to_csv(out / "summary" / f"{m}_mean.csv")
        pd.DataFrame(sds[m], index=names, columns=names).to_csv(out / "summary" / f"{m}_sd.csv")
    long = long_table(arrs, names, dirs)
    long.to_csv(out / "summary" / "pairwise_distances_long.csv", index=False)

    tree_colours = palette(names)
    plot_heatmaps(means, sds, names, n_rep, out)
    if n * (n - 1) // 2 <= 45:
        plot_distributions(long, names, out)
    else:
        print("Skipping pairwise_distributions.png (more than 45 tree pairs).")

    dist_means = {m: D for m, D in means.items() if m in DISTANCES}
    if n >= 3 and dist_means:
        plot_mds(dist_means, names, tree_colours, n_rep, out)
    if n >= 3 and args.cluster in dist_means:
        plot_clustermap(dist_means[args.cluster], names, tree_colours, args.cluster, out)

    print(f"Saved figures and summaries to {out.resolve()}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()