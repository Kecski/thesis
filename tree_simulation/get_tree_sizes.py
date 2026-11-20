#!/usr/bin/env python3
"""Report the size of every tree in a NEXUS .trees file (no dependencies).

Usage: python tree_sizes.py --trees INPUT.trees --output OUTPUT.csv
"""
import argparse, csv, re

def parse_newick(s):
    """Return (n_tips, n_internal, total_branch_length, root_to_tip_max)."""
    s = re.sub(r"\[[^\]]*\]", "", s).strip().rstrip(";")   # drop [&...] annotations
    pos = 0
    stats = {"tips": 0, "internal": 0, "length": 0.0}

    def node():
        nonlocal pos
        depth_below = 0.0
        if s[pos] == "(":
            pos += 1
            stats["internal"] += 1
            while True:
                depth_below = max(depth_below, node())
                if s[pos] == ",":
                    pos += 1
                    continue
                break
            pos += 1                                         # skip ')'
        else:
            stats["tips"] += 1
        m = re.compile(r"[^:,()]*").match(s, pos)            # label (ignored)
        pos = m.end()
        bl = 0.0
        if pos < len(s) and s[pos] == ":":
            m = re.compile(r":([0-9.eE+-]+)").match(s, pos)
            bl = float(m.group(1))
            pos = m.end()
        stats["length"] += bl
        return depth_below + bl

    height = node()
    return stats["tips"], stats["internal"], stats["length"], height

def main():
    ap = argparse.ArgumentParser(description="Write the size of every tree in a NEXUS file to a CSV.")
    ap.add_argument("--trees", help="input NEXUS trees file")
    ap.add_argument("--output", help="output CSV path")
    args = ap.parse_args()

    with open(args.trees) as f:
        text = f.read()
    trees = re.findall(r"^\s*tree\s+(\S+)\s*=\s*(?:\[&[^\]]*\]\s*)?(\(.*;)\s*$",
                       text, flags=re.M | re.I)
    with open(args.output, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["tree", "tips", "internal_nodes", "total_branch_length", "height"])
        for name, nwk in trees:
            tips, internal, tl, h = parse_newick(nwk)
            w.writerow([name, tips, internal, tl, h])
    print(f"Wrote {len(trees)} trees to {args.output}")

if __name__ == "__main__":
    main()