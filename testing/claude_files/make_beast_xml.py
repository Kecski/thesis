#!/usr/bin/env python3
"""
Fill a BEAST2 / BDMM-Prime XML template with sequences, tip dates and parameters.

Usage:
    python make_beast_xml.py config.toml
    python make_beast_xml.py config.toml -o my_analysis.xml
    python make_beast_xml.py config.toml --fasta other.fasta

Needs Python >= 3.11 (uses the built-in tomllib). No other dependencies.

What it does
  1. Reads sequences + tip names from a FASTA file.
  2. Gets a sampling date for every tip, either parsed from the tip name with a
     regular expression or read from a CSV file (name,date).
  3. Builds the <sequence> elements, the date trait and the (untyped) type trait.
  4. Builds one prior per Re interval (different mean/sigma per interval).
  5. Replaces every {{PLACEHOLDER}} in template.xml and writes the result.
"""
import argparse
import csv
import re
import sys
import tomllib
from pathlib import Path
from xml.sax.saxutils import quoteattr

# Package versions written to the <beast required="..."> attribute.
DEFAULT_REQUIRED = "BEAST.base v2.7.8:BDMM-Prime v2.7.2:SA v2.1.1"
IUPAC = set("ACGTUMRWSYKVHDBNacgtumrwsykvhdbn-?")


# ----------------------------------------------------------------------------
# Input readers
# ----------------------------------------------------------------------------
def read_fasta(path):
    """Return a list of (name, sequence). The name is the header up to the first whitespace."""
    records, name, chunks = [], None, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                if name is None:
                    sys.exit(f"ERROR: {path}: sequence data before the first '>' header")
                chunks.append(line)
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def read_dates_csv(path):
    """CSV with a header row; first column = tip name, second column = date (decimal year)."""
    out = {}
    with open(path, newline="", encoding="utf-8") as fh:
        rows = list(csv.reader(fh))
    for row in rows[1:]:
        if len(row) >= 2 and row[0].strip():
            out[row[0].strip()] = row[1].strip()
    return out


# ----------------------------------------------------------------------------
# Validation
# ----------------------------------------------------------------------------
def validate_records(records):
    if not records:
        sys.exit("ERROR: no sequences found in the FASTA file")
    names = [n for n, _ in records]
    dup = {n for n in names if names.count(n) > 1}
    if dup:
        sys.exit(f"ERROR: duplicate tip names: {sorted(dup)[:5]}")
    for n in names:
        if not n or re.search(r"[,=\s\"<>&']", n):
            sys.exit(f"ERROR: tip name {n!r} is empty or contains a character that breaks "
                     "BEAST trait strings (comma, '=', whitespace, quotes, <, >, &)")
    lengths = {len(s) for _, s in records}
    if len(lengths) != 1:
        sys.exit(f"ERROR: sequences are not aligned (lengths differ: {sorted(lengths)[:6]} ...)")
    for n, s in records:
        bad = set(s) - IUPAC
        if bad:
            sys.exit(f"ERROR: sequence {n} has non-nucleotide characters: {sorted(bad)}")


def tip_dates(records, cfg):
    names = [n for n, _ in records]
    dcfg = cfg.get("data", {})
    if dcfg.get("dates_csv"):
        table = read_dates_csv(dcfg["dates_csv"])
        missing = [n for n in names if n not in table]
        if missing:
            sys.exit(f"ERROR: no date in {dcfg['dates_csv']} for {len(missing)} tips, e.g. {missing[:3]}")
        for n in names:
            float(table[n])  # raises if not numeric
        return {n: table[n] for n in names}
    rx = re.compile(dcfg.get("date_regex", r"_(\d{4}(?:\.\d+)?)$"))
    out = {}
    for n in names:
        m = rx.search(n)
        if not m:
            sys.exit(f"ERROR: date_regex {rx.pattern!r} does not match tip name {n!r}")
        out[n] = m.group(1)
    return out


# ----------------------------------------------------------------------------
# Block builders
# ----------------------------------------------------------------------------
def fmt(x):
    """Round-trippable number formatting (1 -> '1.0', 30.42 -> '30.42')."""
    return repr(float(x))


def build_sequences(records, indent=8):
    pad = " " * indent
    lines = []
    for n, s in records:
        lines.append(f'{pad}<sequence id="seq_{n}" spec="Sequence" taxon="{n}" '
                     f'totalcount="4" value="{s}"/>')
    return "\n".join(lines)


def build_re_priors(intervals, indent=16):
    """One <prior> per Re interval, each picking its element of ReEpi with feast's Slice."""
    pad = " " * indent
    blocks = []
    for i, iv in enumerate(intervals, start=1):
        dist = iv.get("distribution", "LogNormal")
        mean, sigma = float(iv["mean"]), float(iv["sigma"])
        if sigma <= 0:
            sys.exit(f"ERROR: interval {i}: sigma must be > 0")
        if dist == "LogNormal":
            if mean <= 0:
                sys.exit(f"ERROR: interval {i}: LogNormal mean (real space) must be > 0")
            body = (
                f'{pad}    <LogNormal id="LogNormalReEpi{i}" meanInRealSpace="true" name="distr">\n'
                f'{pad}        <parameter id="ReEpiPriorM{i}" spec="parameter.RealParameter" '
                f'estimate="false" name="M">{fmt(mean)}</parameter>\n'
                f'{pad}        <parameter id="ReEpiPriorS{i}" spec="parameter.RealParameter" '
                f'estimate="false" lower="0.0" upper="5.0" name="S">{fmt(sigma)}</parameter>\n'
                f'{pad}    </LogNormal>\n')
        elif dist == "Normal":
            body = (
                f'{pad}    <Normal id="NormalReEpi{i}" name="distr">\n'
                f'{pad}        <parameter id="ReEpiPriorMean{i}" spec="parameter.RealParameter" '
                f'estimate="false" name="mean">{fmt(mean)}</parameter>\n'
                f'{pad}        <parameter id="ReEpiPriorSigma{i}" spec="parameter.RealParameter" '
                f'estimate="false" lower="0.0" name="sigma">{fmt(sigma)}</parameter>\n'
                f'{pad}    </Normal>\n')
        else:
            sys.exit(f"ERROR: interval {i}: unknown distribution {dist!r} (use LogNormal or Normal)")
        blocks.append(
            f'{pad}<prior id="ReEpiPrior{i}.t:sequences" name="distribution">\n'
            f'{pad}    <x id="ReEpiSlice{i}.t:sequences" spec="feast.function.Slice" '
            f'arg="@ReEpi.t:sequences" index="{i - 1}" count="1"/>\n'
            f'{body}'
            f'{pad}</prior>')
    # the template already indents the first line, so strip it here
    return "\n".join(blocks).lstrip()


# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="TOML config file")
    ap.add_argument("-o", "--output", help="output XML (overrides config)")
    ap.add_argument("--fasta", help="FASTA file (overrides config)")
    ap.add_argument("--template", help="template XML (overrides config)")
    args = ap.parse_args()

    cfg_path = Path(args.config)
    with open(cfg_path, "rb") as fh:
        cfg = tomllib.load(fh)
    base = cfg_path.parent  # relative paths in the config are relative to the config file

    def path(p):
        p = Path(p)
        return p if p.is_absolute() else base / p

    data, mcmc, model = cfg.get("data", {}), cfg.get("mcmc", {}), cfg.get("model", {})
    subst, clock, samp = cfg.get("substitution", {}), cfg.get("clock", {}), cfg.get("sampling", {})
    intervals = cfg.get("re_interval", [])
    if not intervals:
        sys.exit("ERROR: define at least one [[re_interval]] block in the config")

    # paths given on the command line are relative to the current directory,
    # paths from the config file are relative to the config file
    fasta = Path(args.fasta) if args.fasta else (path(data["fasta"]) if data.get("fasta") else None)
    if not fasta:
        sys.exit("ERROR: no FASTA file given (config [data] fasta, or --fasta)")
    template = Path(args.template) if args.template else path(cfg.get("template", "template.xml"))
    for label, f in (("FASTA", fasta), ("template", template)):
        if not Path(f).is_file():
            sys.exit(f"ERROR: {label} file not found: {f}")
    out = Path(args.output) if args.output else path(cfg.get("output", "analysis.xml"))
    if data.get("dates_csv"):
        data["dates_csv"] = str(path(data["dates_csv"]))

    # ---- sequences and dates
    records = read_fasta(fasta)
    validate_records(records)
    dates = tip_dates(records, cfg)

    # ---- Re skyline
    n_int = len(intervals)
    change_times = [float(t) for t in model.get("re_change_times", [])]
    if len(change_times) != n_int - 1:
        sys.exit(f"ERROR: {n_int} [[re_interval]] blocks need {n_int - 1} re_change_times, "
                 f"got {len(change_times)}")
    if any(b <= a for a, b in zip(change_times, change_times[1:])) or any(t <= 0 for t in change_times):
        sys.exit("ERROR: re_change_times must be positive and strictly increasing")
    re_upper = float(model.get("re_upper", 10.0))
    re_init = [float(iv.get("initial", iv["mean"])) for iv in intervals]
    if any(not 0 < v <= re_upper for v in re_init):
        sys.exit(f"ERROR: initial Re values must lie in (0, {re_upper}]: {re_init}")

    # ---- sanity warnings
    dnum = [float(v) for v in dates.values()]
    span = max(dnum) - min(dnum)
    origin = float(model.get("origin_initial", 1.0))
    if origin <= span:
        print(f"WARNING: origin_initial ({origin}) <= sampling span ({span:.3f}); "
              "BEAST may fail to initialise the tree. Increase origin_initial.", file=sys.stderr)
    if change_times and origin <= change_times[-1]:
        print(f"WARNING: origin_initial ({origin}) <= last Re change time ({change_times[-1]}).",
              file=sys.stderr)

    required = DEFAULT_REQUIRED
    if cfg.get("feast_version"):
        required += f":feast v{cfg['feast_version']}"

    # ---- placeholder values
    names = [n for n, _ in records]
    values = {
        "REQUIRED": required,
        "SEQUENCES": build_sequences(records),
        "DATE_TRAIT_VALUE": ",".join(f"{n}={dates[n]}" for n in names),
        "TYPE_TRAIT_VALUE": ",".join(f"{n}=NOT_SET" for n in names),
        "CHAIN_LENGTH": int(mcmc.get("chain_length", 30_000_000)),
        "STORE_EVERY": int(mcmc.get("store_every", 5000)),
        "LOG_EVERY": int(mcmc.get("log_every", 1000)),
        "CLOCK_RATE_INITIAL": fmt(clock.get("initial", 1.0)),
        "CLOCK_PRIOR_M": fmt(clock.get("prior_mean", 0.001)),
        "CLOCK_PRIOR_S": fmt(clock.get("prior_sigma", 1.25)),
        "ORIGIN_INITIAL": fmt(origin),
        "RE_DIMENSION": n_int,
        "RE_UPPER": fmt(re_upper),
        "RE_INITIAL": " ".join(fmt(v) for v in re_init),
        "RE_CHANGE_TIMES_DIM": len(change_times),
        "RE_CHANGE_TIMES": " ".join(fmt(t) for t in change_times),
        "RE_PRIORS": build_re_priors(intervals),
        "BECOME_UNINFECTIOUS_RATE": fmt(model["become_uninfectious_rate"]),
        "SAMPLING_PROPORTION_INITIAL": fmt(samp.get("initial", 0.1)),
        "SAMPLING_PRIOR_M": fmt(samp.get("prior_m", 1.0)),
        "SAMPLING_PRIOR_S": fmt(samp.get("prior_s", 1.25)),
        "KAPPA": fmt(subst.get("kappa", 3.77)),
        "GAMMA_SHAPE": fmt(subst.get("gamma_shape", 0.5)),
        "GAMMA_CATEGORIES": int(subst.get("gamma_categories", 6)),
    }

    # ---- fill template
    xml = Path(template).read_text(encoding="utf-8")
    wanted = set(re.findall(r"\{\{(\w+)\}\}", xml))
    missing = wanted - set(values)
    if missing:
        sys.exit(f"ERROR: template has placeholders the script does not know: {sorted(missing)}")
    unused = set(values) - wanted
    if unused:
        print(f"WARNING: values not used by the template: {sorted(unused)}", file=sys.stderr)
    for key in wanted:
        xml = xml.replace("{{" + key + "}}", str(values[key]))
    assert "{{" not in xml

    out.write_text(xml, encoding="utf-8")

    # ---- summary
    print(f"Wrote {out}")
    print(f"  sequences : {len(records)} x {len(records[0][1])} sites")
    print(f"  tip dates : {min(dnum):.4f} - {max(dnum):.4f} (span {span:.4f})")
    print(f"  Re skyline: {n_int} intervals, change times (ages) = {change_times}")
    for i, iv in enumerate(intervals, start=1):
        print(f"    interval {i} (index {i - 1}): {iv.get('distribution', 'LogNormal')} "
              f"mean={iv['mean']} sigma={iv['sigma']} initial={re_init[i - 1]}")


if __name__ == "__main__":
    main()
