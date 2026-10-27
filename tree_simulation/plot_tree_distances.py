#!/usr/bin/env python3
"""
Visualize the distance matrices written by tree_distances.py (--csv PREFIX).

Produces, in --outdir:
  heatmaps.png        one heatmap per metric
  mds.png             2D "tree space" map (classical MDS) per distance metric
  clustered_<m>.png   heatmap + dendrogram for one metric (--cluster, default mrca_rmse)

Usage
  python plot_distances.py out                           # reads out_<metric>.csv
  python plot_distances.py out --labels names.txt        # one tree name per line
  python plot_distances.py out --metrics rf_norm bsd --show

Requires: pip install numpy matplotlib scipy
"""
import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.spatial.distance import squareform

DISTANCES = ["rf", "rf_norm", "wrf", "bsd", "triplet",
             "mrca_rmse", "mrca_mad", "clade_age_rmse", "clade_age_mad"]
SIMILARITIES = ["n_shared"]  # larger = more similar, so excluded from MDS/clustering


def load(input_dir, prefix, wanted):
    mats = {}
    for m in wanted or DISTANCES + SIMILARITIES:
        input_path = Path(input_dir)
        f = input_path / f"{prefix}_{m}.csv"
        if f.exists():
            mats[m] = np.atleast_2d(np.loadtxt(f, delimiter=","))
    if not mats:
        raise SystemExit(f"No files like {prefix}_<metric>.csv found.")
    return mats


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
    order = np.argsort(vals)[::-1][:k]
    vals, vecs = vals[order], vecs[:, order]
    coords = vecs * np.sqrt(np.clip(vals, 0, None))
    pos_total = np.clip(np.linalg.eigvalsh(-0.5 * J @ (D ** 2) @ J), 0, None).sum()
    return coords, (np.clip(vals, 0, None) / pos_total if pos_total else vals * 0)


def grid(n_plots, ncols=3, size=4.2):
    nrows = int(np.ceil(n_plots / ncols))
    fig, axes = plt.subplots(nrows, ncols, figsize=(ncols * size, nrows * size),
                             squeeze=False)
    axes = axes.ravel()
    for ax in axes[n_plots:]:
        ax.axis("off")
    return fig, axes


def heatmap(ax, D, labels, title, cmap):
    n = len(D)
    im = ax.imshow(np.ma.masked_invalid(D), cmap=cmap)
    ax.set_title(title)
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    if n <= 30:
        ax.set_xticklabels(labels, rotation=90, fontsize=7)
        ax.set_yticklabels(labels, fontsize=7)
    else:
        ax.set_xticks([])
        ax.set_yticks([])
    if n <= 10:
        for i in range(n):
            for j in range(n):
                if np.isfinite(D[i, j]):
                    ax.text(j, i, f"{D[i, j]:.2g}", ha="center", va="center",
                            fontsize=6, color="white" if D[i, j] > np.nanmean(D) else "black")
    plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("prefix", help="the PREFIX given to tree_distances.py --csv")
    ap.add_argument("--indir", default=".", help="Directory of the csvs.")
    ap.add_argument("--metrics", nargs="+", help="subset of metrics to plot")
    ap.add_argument("--labels", help="text file with one tree name per line")
    ap.add_argument("--cluster", default="mrca_rmse", help="metric for the dendrogram plot")
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--show", action="store_true", help="also open the figures")
    args = ap.parse_args()

    mats = load(args.indir, args.prefix, args.metrics)
    n = len(next(iter(mats.values())))
    if args.labels:
        labels = [l.strip() for l in open(args.labels) if l.strip()]
        if len(labels) != n:
            raise SystemExit(f"{len(labels)} labels but {n} trees.")
    else:
        labels = [str(i + 1) for i in range(n)]
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    # 1. heatmaps
    fig, axes = grid(len(mats))
    for ax, (m, D) in zip(axes, mats.items()):
        if m in SIMILARITIES:  # the diagonal is 0 by construction; hide it so it doesn't skew colours
            D = D.astype(float)
            np.fill_diagonal(D, np.nan)
        heatmap(ax, D, labels, m, "viridis_r" if m in SIMILARITIES else "viridis")
    fig.tight_layout()
    fig.savefig(out / "heatmaps.png", dpi=200)

    # 2. MDS tree-space maps (needs >= 3 trees)
    dist_mats = {m: D for m, D in mats.items() if m in DISTANCES}
    if n >= 3 and dist_mats:
        fig, axes = grid(len(dist_mats))
        for ax, (m, D) in zip(axes, dist_mats.items()):
            xy, var = classical_mds(fill_nan(D))
            ax.scatter(xy[:, 0], xy[:, 1], s=40)
            for (x, y), lab in zip(xy, labels):
                ax.annotate(lab, (x, y), xytext=(3, 3), textcoords="offset points", fontsize=7)
            ax.set_title(m)
            ax.set_xlabel(f"MDS 1 ({var[0]:.0%})")
            ax.set_ylabel(f"MDS 2 ({var[1]:.0%})")
        fig.tight_layout()
        fig.savefig(out / "mds.png", dpi=200)

    # 3. clustered heatmap with dendrogram
    if n >= 3 and args.cluster in dist_mats:
        D = fill_nan(dist_mats[args.cluster])
        Z = linkage(squareform(D, checks=False), method="average")
        fig = plt.figure(figsize=(7, 6.5))
        axd = fig.add_axes([0.05, 0.08, 0.17, 0.78])
        tree = dendrogram(Z, orientation="left", ax=axd, no_labels=True,
                          color_threshold=0, above_threshold_color="gray")
        axd.set_ylim(0, 10 * n)  # leaf spacing is 10, so tips line up with heatmap rows
        axd.axis("off")
        order = tree["leaves"][::-1]
        axh = fig.add_axes([0.25, 0.08, 0.6, 0.78])
        im = axh.imshow(D[np.ix_(order, order)], cmap="viridis")
        ticks = range(n)
        axh.set_xticks(ticks)
        axh.set_yticks(ticks)
        axh.set_xticklabels([labels[i] for i in order], rotation=90, fontsize=7)
        axh.set_yticklabels([labels[i] for i in order], fontsize=7)
        fig.colorbar(im, cax=fig.add_axes([0.88, 0.08, 0.02, 0.78]))
        fig.suptitle(f"{args.cluster} (average-linkage clustering)")
        fig.savefig(out / f"clustered_{args.cluster}.png", dpi=200)

    print(f"Saved figures to {out.resolve()}")
    if args.show:
        plt.show()


if __name__ == "__main__":
    main()
