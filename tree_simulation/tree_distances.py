#!/usr/bin/env python3
"""
Pairwise distances between rooted phylogenies that share the same taxa.

Metrics
  rf        Robinson-Foulds (rooted, clade-based) and normalized RF
  wrf       Weighted RF (sum of absolute branch-length differences per clade)
  bsd       Branch score distance (Kuhner-Felsenstein; Euclidean on branch lengths)
  triplet   Rooted triplet distance (fraction of 3-taxon subsets that differ)

  Time-aware metrics (use node heights/ages derived from branch lengths)
  mrca_rmse, mrca_mad
            Compare, over all taxon pairs, the height/age of the pair's MRCA in the
            two trees (RMSE and mean absolute difference). Captures topology and timing.
  clade_age_rmse, clade_age_mad, n_shared
            Compare node heights/ages only for clades present in BOTH trees, and the
            number of such clades. Isolates timing error from topology error.

  --anchor root (default): heights = distance from the root.
  --anchor tips: ages = (distance of the latest tip from the root) - height, i.e. time
                 before the most recent sample. Use this when root positions differ.

Usage
  python tree_distances.py trees.nwk                # one Newick tree per line
  python tree_distances.py trees.nex                # Nexus (auto-detected from #NEXUS)
  python tree_distances.py trees.tre --format nexus # force a format
  python tree_distances.py a.nwk b.nwk c.nwk        # several files (all trees pooled)
  python tree_distances.py trees.nwk --csv out      # also write out_<metric>.csv

Requires: pip install dendropy numpy
"""
import argparse
import itertools
import math
from pathlib import Path

import dendropy
import numpy as np


# ----------------------------------------------------------------- loading
def detect_schema(path):
    """Return 'nexus' if the file starts with #NEXUS, otherwise 'newick'."""
    with open(path, "r", errors="replace") as fh:
        for line in fh:
            if line.strip():
                return "nexus" if line.strip().upper().startswith("#NEXUS") else "newick"
    raise SystemExit(f"{path} is empty.")


def load_trees(paths, schema="auto"):
    """Read all trees from the given Newick/Nexus files into one shared namespace.
    schema: 'auto' (sniff each file), 'newick' or 'nexus'."""
    tns = dendropy.TaxonNamespace()
    trees = dendropy.TreeList(taxon_namespace=tns)
    for p in paths:
        fmt = detect_schema(p) if schema == "auto" else schema
        trees.read(path=p, schema=fmt, rooting="force-rooted",
                   taxon_namespace=tns, preserve_underscores=True)
    if len(trees) < 2:
        raise SystemExit("Need at least two trees.")
    leafsets = [frozenset(l.taxon.label for l in t.leaf_node_iter()) for t in trees]
    if len(set(leafsets)) != 1:
        raise SystemExit("Trees do not all have the same taxa.")
    return trees


# -------------------------------------------------------- tree summaries
def clade_lengths(tree):
    """Map each non-root clade (frozenset of leaf labels) -> its branch length.
    Includes terminal branches (singleton clades)."""
    clades = {}
    leaves = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            leaves[node] = frozenset([node.taxon.label])
        else:
            leaves[node] = frozenset().union(*(leaves[c] for c in node.child_node_iter()))
        if node is not tree.seed_node:
            clades[leaves[node]] = node.edge_length or 0.0
    return clades


def lca_depth_matrix(tree, labels, use_lengths=False):
    """D[i, j] = depth of the MRCA of taxa i and j, counted in edges from the root,
    or as root-to-node branch-length distance if use_lengths=True."""
    idx = {l: i for i, l in enumerate(labels)}
    n = len(labels)
    D = np.zeros((n, n), dtype=float if use_lengths else int)
    members = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            members[node] = [idx[node.taxon.label]]
            continue
        groups = [members[c] for c in node.child_node_iter()]
        depth = node.distance_from_root() if use_lengths else node.level()
        for g1, g2 in itertools.combinations(groups, 2):
            D[np.ix_(g1, g2)] = depth
            D[np.ix_(g2, g1)] = depth
        members[node] = [i for g in groups for i in g]
    return D


def triplet_topologies(D):
    """For every triple (i<j<k) return 0,1,2 for the resolved pair
    (ij, ik, jk) or -1 if the triple is an unresolved star."""
    n = D.shape[0]
    out = []
    for i, j, k in itertools.combinations(range(n), 3):
        d = (D[i, j], D[i, k], D[j, k])
        m = max(d)
        out.append(d.index(m) if d.count(m) == 1 else -1)
    return np.array(out)


def tree_height(tree):
    """Largest root-to-tip distance."""
    return max(l.distance_from_root() for l in tree.leaf_node_iter())


def clade_node_values(tree, anchor="root"):
    """Map each internal clade (frozenset of leaf labels) -> node height (anchor='root',
    root-to-node distance, root excluded since it is always 0) or node age
    (anchor='tips', tree_height - height, root included)."""
    tmax = tree_height(tree)
    out, leaves = {}, {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            leaves[node] = frozenset([node.taxon.label])
            continue
        leaves[node] = frozenset().union(*(leaves[c] for c in node.child_node_iter()))
        h = node.distance_from_root()
        if anchor == "tips":
            out[leaves[node]] = tmax - h
        elif node is not tree.seed_node:
            out[leaves[node]] = h
    return out


def mrca_values(tree, labels, anchor="root"):
    """Vector over all taxon pairs (i<j) of the MRCA height (or age if anchor='tips')."""
    H = lca_depth_matrix(tree, labels, use_lengths=True)
    if anchor == "tips":
        H = tree_height(tree) - H
    return H[np.triu_indices(len(labels), 1)]


# --------------------------------------------------------------- distances
def rf_distance(c1, c2):
    """Rooted RF on non-trivial clades. Returns (RF, normalized RF)."""
    a = {c for c in c1 if len(c) > 1}
    b = {c for c in c2 if len(c) > 1}
    rf = len(a ^ b)
    denom = len(a) + len(b)
    return rf, (rf / denom if denom else 0.0)


def weighted_rf(c1, c2):
    keys = set(c1) | set(c2)
    return sum(abs(c1.get(k, 0.0) - c2.get(k, 0.0)) for k in keys)


def branch_score(c1, c2):
    keys = set(c1) | set(c2)
    return math.sqrt(sum((c1.get(k, 0.0) - c2.get(k, 0.0)) ** 2 for k in keys))


def triplet_distance(t1, t2):
    """Fraction of triplets with different topology (0 = identical, 1 = all differ)."""
    return float(np.mean(t1 != t2))


def mrca_distance(v1, v2):
    """(RMSE, mean absolute difference) of pairwise MRCA heights/ages."""
    d = v1 - v2
    return float(np.sqrt(np.mean(d ** 2))), float(np.mean(np.abs(d)))


def shared_clade_age_distance(a1, a2):
    """(RMSE, MAD, number of shared clades) of node heights/ages on clades in both trees.
    RMSE and MAD are NaN if the trees share no clades."""
    keys = set(a1) & set(a2)
    if not keys:
        return float("nan"), float("nan"), 0
    d = np.array([a1[k] - a2[k] for k in keys])
    return float(np.sqrt(np.mean(d ** 2))), float(np.mean(np.abs(d))), len(keys)


# -------------------------------------------------------------------- main
def all_pairwise(trees, anchor="root"):
    n = len(trees)
    labels = sorted(l.taxon.label for l in trees[0].leaf_node_iter())
    clades = [clade_lengths(t) for t in trees]
    trip = [triplet_topologies(lca_depth_matrix(t, labels)) for t in trees]
    mrca = [mrca_values(t, labels, anchor) for t in trees]
    ages = [clade_node_values(t, anchor) for t in trees]

    names = ["rf", "rf_norm", "wrf", "bsd", "triplet",
             "mrca_rmse", "mrca_mad", "clade_age_rmse", "clade_age_mad", "n_shared"]
    M = {k: np.zeros((n, n)) for k in names}
    for i, j in itertools.combinations(range(n), 2):
        rf, rfn = rf_distance(clades[i], clades[j])
        m_rmse, m_mad = mrca_distance(mrca[i], mrca[j])
        c_rmse, c_mad, n_sh = shared_clade_age_distance(ages[i], ages[j])
        vals = {"rf": rf, "rf_norm": rfn,
                "mrca_rmse": m_rmse, "mrca_mad": m_mad,
                "clade_age_rmse": c_rmse, "clade_age_mad": c_mad, "n_shared": n_sh,
                "wrf": weighted_rf(clades[i], clades[j]),
                "bsd": branch_score(clades[i], clades[j]),
                "triplet": triplet_distance(trip[i], trip[j])}
        for k, v in vals.items():
            M[k][i, j] = M[k][j, i] = v
    return M


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("files", nargs="+", help="Newick or Nexus file(s)")
    ap.add_argument("--output_dir", type=str, default=".",
                    help="Output directory of the csvs.")
    ap.add_argument("--format", choices=["auto", "newick", "nexus"], default="auto",
                    help="input format (default: auto-detect each file)")
    ap.add_argument("--anchor", choices=["root", "tips"], default="root",
                    help="time-aware metrics: use heights from the root (default) or "
                         "ages before the latest tip")
    ap.add_argument("--csv", metavar="PREFIX", help="write PREFIX_<metric>.csv matrices")
    args = ap.parse_args()

    trees = load_trees(args.files, args.format)
    n = len(trees)
    print(f"{n} trees, {len(trees.taxon_namespace)} taxa\n")
    M = all_pairwise(trees, args.anchor)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    for k, mat in M.items():
        print(f"== {k} ==")
        print(np.array2string(mat, precision=4, suppress_small=True), "\n")
        if args.csv:
            np.savetxt(f"{output_dir / args.csv}_{k}.csv", mat, delimiter=",", fmt="%.6g")


if __name__ == "__main__":
    main()