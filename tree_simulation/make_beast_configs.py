"""
Create BEAST config files (for make_beast_xml.py) from EpiSewer summary csvs.

Input layout:
  input_dir/tree_0/episewer_summary.csv
  input_dir/tree_1/episewer_summary.csv
  ...
Output layout (two configs per tree):
  output_dir/tree_0/output_sections.toml      one Re prior per section; re_change_times
                                              spread evenly over --simulation-length
  output_dir/tree_0/output_whole_period.toml  one constant-Re prior (re_change_times = [])
  output_dir/tree_1/...

The template is edited as text, so comments and formatting are preserved.
Per tree, the template's tree index in the fasta path (e.g. tree0 or tree_0) is
replaced by that tree's index, and the output xml name gets a tree suffix.

Usage:
  python make_configs_from_episewer.py \
      --input_dir  data/wastewater \
      --output_dir configs \
      --config     config_beast_xml.toml \
      --simulation-length 0.3
"""
import argparse
import csv
import math
import os
import re
import sys
from pathlib import Path
from typing import Optional


def read_R_rows(path):
    """Return (list of section rows sorted by start_date, whole_period row)."""
    with open(path, newline="") as f:
        rows = [r for r in csv.DictReader(f) if r["variable"] == "R"]
    if not rows:
        sys.exit(f"No rows with variable == 'R' in {path}")
    whole = [r for r in rows if r["section"] == "whole_period"]
    secs = sorted((r for r in rows if r["section"] != "whole_period"),
                  key=lambda r: r["start_date"])
    if not whole:
        sys.exit(f"No 'whole_period' row for R in {path} "
                 "(re-run episewer.R with the whole-period addition).")
    return secs, whole[0]


def to_prior_params(mean, sd, sigma_mode, mean_in_real_space):
    """Convert posterior mean/sd of Re into the (mean, sigma) the config expects."""
    if sigma_mode == "raw":
        sigma = sd
    else:  # "lognormal": moment-match a LogNormal with this mean and SD on the Re scale
        sigma = math.sqrt(math.log(1.0 + (sd / mean) ** 2))
    m = mean if mean_in_real_space else math.log(mean) - sigma ** 2 / 2.0
    return m, sigma


def fmt(values):
    vals = [f"{v:.4f}".rstrip("0").rstrip(".") for v in values]
    return vals[0] if len(vals) == 1 else "[" + ", ".join(vals) + "]"

def fmt_times(ts):
    return "[" + ", ".join(f"{t:.6f}".rstrip("0").rstrip(".") for t in ts) + "]"


def even_change_times(n_intervals, simulation_length):
    """n intervals of equal length over simulation_length -> n-1 change times (ages)."""
    ts = [simulation_length * k / n_intervals for k in range(1, n_intervals)]
    return "[" + ", ".join(f"{t:.6f}".rstrip("0").rstrip(".") for t in ts) + "]", ts


def replace_key(block, key, new_value, required=True):
    """Replace `key = ...` (keeping any trailing comment)."""
    pat = re.compile(rf"^(\s*{key}\s*=\s*)([^#\n]*?)(\s*(#.*)?)$", re.M)
    if not pat.search(block):
        if required:
            sys.exit(f"Key '{key}' not found in the template.")
        return block
    return pat.sub(lambda m: f"{m.group(1)}{new_value}{m.group(3)}", block, count=1)


def get_string_key(text, key):
    m = re.search(rf'^\s*{key}\s*=\s*"([^"]*)"', text, re.M)
    return m.group(1) if m else None


def split_section(text, header):
    """Return (before, block, after) for the [header] table."""
    m = re.search(rf"^\[{re.escape(header)}\]\s*$", text, re.M)
    if not m:
        sys.exit(f"Section [{header}] not found in the template.")
    nxt = re.search(r"^\[", text[m.end():], re.M)
    end = m.end() + nxt.start() if nxt else len(text)
    return text[:m.start()], text[m.start():end], text[end:]


def relocate(path, template_dir, new_dir):
    """Paths in the config are relative to the config file: re-express a path that was
    relative to the template's folder so it is relative to the new config's folder."""
    if os.path.isabs(path):
        return path
    return os.path.relpath((template_dir / path).resolve(), new_dir.resolve()).replace(os.sep, "/")


def build(template, template_dir, new_dir, tree_idx, tag, means, sigmas,
          change_times, set_initial, relocate_paths):
    # [re_prior]
    before, block, after = split_section(template, "re_prior")
    block = replace_key(block, "mean", fmt(means))
    block = replace_key(block, "sigma", fmt(sigmas))
    if set_initial:
        block = replace_key(block, "initial", fmt(means))
    text = before + block + after
 
    text = replace_key(text, "re_change_times", change_times)
 
    # paths, expressed relative to the template's folder first
    paths = {k: get_string_key(text, k) for k in ("template", "output", "fasta", "dates_csv")}
    if paths["output"] and tree_idx is not None:
        o = Path(paths["output"])
        paths["output"] = (o.parent / f"{o.stem}_tree{tree_idx}_{tag}{o.suffix}").as_posix()
    if paths["fasta"] and tree_idx is not None:
        # tree0 / tree_0 -> tree{idx} / tree_{idx}
        paths["fasta"] = re.sub(r"(tree_?)\d+", rf"\g<1>{tree_idx}", paths["fasta"])
    for k, p in paths.items():
        if p is None:
            continue
        if relocate_paths:
            p = relocate(p, template_dir, new_dir)
        text = replace_key(text, k, f'"{p}"', required=False)
    return text, paths["output"], paths["fasta"]


def find_trees(input_dir, tree_id):
    trees = []
    for d in Path(input_dir).iterdir():
        m = re.fullmatch(r"tree_(\d+)", d.name)
        if m and d.is_dir():
            trees.append((int(m.group(1)), d))
    return sorted(trees)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input_dir", "--input-dir", default=None,
                    help="folder containing tree_<i>/episewer_summary.csv "
                         "(omit when giving --re-means etc. manually)")
    ap.add_argument("--output_dir", "--output-dir", required=True,
                    help="folder in which tree_<i>/ with the new configs is created")
    ap.add_argument("--tree_id", "--tree-id", type=int, default=0,
                    help="ID ([0, N-1]) of the tree if there are N trees in total.")
    ap.add_argument("--config", required=True, help="template config toml")
    ap.add_argument("--re-means", type=float, nargs="+", default=None, metavar="X",
                    help="MANUAL MODE: Re prior means, one per interval")
    ap.add_argument("--re-sigmas", type=float, nargs="+", default=None, metavar="X",
                    help="MANUAL MODE: Re prior sigmas, one per interval (or a single "
                         "value used for all)")
    ap.add_argument("--re-change-times", type=float, nargs="+", default=None, metavar="T",
                    help="MANUAL MODE: Re change times (ages in years), n-1 for n means; "
                         "alternatively use --simulation-length")
    ap.add_argument("--output_name", "--output-name", default="output.toml",
                    help="MANUAL MODE: file name of the config inside output_dir")
    ap.add_argument("--summary-name", default="episewer_summary.csv",
                    help="csv file name inside each tree folder")
    ap.add_argument("--simulation-length", "--simulation_length", type=float, default=None,
                    help="length of the simulation in years; re_change_times in the "
                         "'sections' config are set evenly over it (n sections -> n-1 "
                         "change times, e.g. 3 sections and 0.3 -> [0.1, 0.2]). If omitted, "
                         "the template's re_change_times are kept.")
    ap.add_argument("--sigma-mode", choices=["lognormal", "raw"], default="lognormal",
                    help="'lognormal' (default): convert the posterior SD of Re into the "
                         "log-scale sigma by moment matching. 'raw': write the csv sd as is.")
    ap.add_argument("--keep-initial", action="store_true",
                    help="do not touch the 'initial' MCMC starting values "
                         "(default: set them to the new means)")
    ap.add_argument("--reverse-sections", action="store_true",
                    help="reverse the order of sections (use if index 0 in the config is "
                         "the MOST RECENT interval)")
    ap.add_argument("--keep-paths", action="store_true",
                    help="do not re-express relative paths (template/output/fasta/dates_csv) "
                         "relative to the new config's folder")
    args = ap.parse_args()
    
    if args.re_change_times == [-1]:
        args.re_change_times = None
 
    config_path = Path(args.config)
    template = config_path.read_text()
    template_dir = config_path.resolve().parent

    # prior settings in the template
    _, prior_block, _ = split_section(template, "re_prior")
    mirs = re.search(r"^\s*mean_in_real_space\s*=\s*(\w+)", prior_block, re.M)
    mean_in_real_space = (mirs.group(1).lower() == "true") if mirs else True
    dist = re.search(r'^\s*distribution\s*=\s*"(\w+)"', prior_block, re.M)
    sigma_mode = args.sigma_mode
    if dist and dist.group(1) != "LogNormal":
        print("Warning: template prior is not LogNormal; mean/sigma are written "
              "from the csv without conversion.", file=sys.stderr)
        mean_in_real_space, sigma_mode = True, "raw"
 
    if args.simulation_length is not None and args.simulation_length <= 0:
        sys.exit("--simulation-length must be positive.")
    og = re.search(r"^origin_initial\s*=\s*([\d.eE+-]+)", template, re.M)
    ct = re.search(r"^re_change_times\s*=\s*\[([^\]]*)\]", template, re.M)
    if not ct:
        sys.exit("re_change_times not found in the template.")

    if args.re_means is not None or args.re_sigmas is not None \
            or args.re_change_times is not None:
        # ---------------- manual mode ----------------
        if args.input_dir:
            sys.exit("Use either --input_dir or the manual --re-means/--re-sigmas/"
                     "--re-change-times, not both.")
        if args.re_means is None or args.re_sigmas is None:
            sys.exit("Manual mode needs both --re-means and --re-sigmas.")
        means, sigmas = args.re_means, args.re_sigmas
        n = len(means)
        if len(sigmas) not in (1, n):
            sys.exit(f"--re-sigmas needs 1 or {n} values (got {len(sigmas)}).")
        if args.re_change_times is not None:
            ts = args.re_change_times
            if len(ts) != n - 1:
                sys.exit(f"{n} means need {n - 1} change times (got {len(ts)}).")
            if any(b <= a for a, b in zip(ts, ts[1:])) or (ts and ts[0] <= 0):
                sys.exit("--re-change-times must be positive and strictly increasing.")
            change_times = fmt_times(ts)
        elif args.simulation_length is not None:
            change_times, ts = even_change_times(n, args.simulation_length)
        elif n == 1:
            change_times, ts = "[]", []
        else:
            sys.exit(f"{n} means need --re-change-times ({n - 1} values) or "
                     "--simulation-length.")
        if ts and og and ts[-1] >= float(og.group(1)):
            print(f"Warning: last change time {ts[-1]} >= origin_initial "
                  f"{og.group(1)}; origin_initial must be larger.", file=sys.stderr)
        new_dir = Path(args.output_dir)
        new_dir.mkdir(parents=True, exist_ok=True)
        text, _, _ = build(template, template_dir, new_dir, None, None, means, sigmas,
                           change_times, set_initial=not args.keep_initial,
                           relocate_paths=not args.keep_paths)
        dest = new_dir / args.output_name
        dest.write_text(text)
        print(f"{dest}: mean={fmt(means)} sigma={fmt(sigmas)} change_times={change_times}")
        return


    tree_dir = Path(args.input_dir) / f"tree_{args.tree_id}"

    def params(rows):
        out = [to_prior_params(float(r["mean"]), float(r["sd"]),
                               sigma_mode, mean_in_real_space) for r in rows]
        return [o[0] for o in out], [o[1] for o in out]

    # for idx, tree_dir in trees:
    csv_path = tree_dir / args.summary_name
    # if not csv_path.exists():
    #     print(f"[tree_{idx}] skipped: {csv_path} not found", file=sys.stderr)
    #     continue
    secs, whole = read_R_rows(csv_path)
    if args.reverse_sections:
        secs = secs[::-1]

    if args.simulation_length is not None:
        change_times_sections, ts = even_change_times(len(secs), args.simulation_length)
        if ts and og and ts[-1] >= float(og.group(1)):
            print(f"[tree_{args.tree_id}] Warning: last change time {ts[-1]} >= origin_initial "
                    f"{og.group(1)}; origin_initial must be larger.", file=sys.stderr)
    else:
        n_template = len([x for x in ct.group(1).split(",") if x.strip()]) + 1
        if n_template != len(secs):
            sys.exit(f"[tree_{args.tree_id}] The csv has {len(secs)} sections for R but the "
                        f"template defines {n_template} intervals (re_change_times). "
                        "Adjust one of them or pass --simulation-length.")
        change_times_sections = "[" + ct.group(1).strip() + "]"

    new_dir = Path(args.output_dir) / f"tree_{args.tree_id}"
    new_dir.mkdir(parents=True, exist_ok=True)

    variants = [
        ("sections", params(secs), change_times_sections),
        ("whole_period", params([whole]), "[]"),
    ]
    for tag, (means, sigmas), change_times in variants:
        text, out_xml, fasta = build(template, template_dir, new_dir, args.tree_id, tag,
                                        means, sigmas, change_times,
                                        set_initial=not args.keep_initial,
                                        relocate_paths=not args.keep_paths)
        dest = new_dir / f"config_beast_xml_{tag}.toml"
        dest.write_text(text)
        print(f"{dest}: mean={fmt(means)} sigma={fmt(sigmas)} "
                f"change_times={change_times} fasta={fasta}")

    print(f"Done: {args.tree_id} written to {args.output_dir}")


if __name__ == "__main__":
    main()
