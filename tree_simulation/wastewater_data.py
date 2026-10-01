#!/usr/bin/env python3
"""
Command-line interface: simulate wastewater measurements from a full,
time-scaled phylogeny using the functions in `wastewater_utils.py`.

Example
-------
python run_simulation.py tree.nwk --start-date 2020-03-01 \
    --branch-length-unit-days 365.25 --sample-days-every 2 \
    --n-replicates 2 --seed 42 --output sim.csv --episewer-prefix sim_ep

Run `python wastewater_data.py -h` for all options.
"""
import argparse
import dataclasses
import math
import sys
import pandas as pd
from pathlib import Path
from io import StringIO
from Bio import Phylo
import dendropy

import wastewater_data_utils as wu

# help texts for the model parameters (defaults are appended automatically)
PARAM_HELP = {
    "total_load_per_infection": "Total gene copies shed per infected person.",
    "shed_mean": "Mean (days) of the gamma shedding profile since infection.",
    "shed_sd": "SD (days) of the gamma shedding profile.",
    "shed_max_days": "Max days after infection that shedding is modelled.",
    "indiv_shape": "Gamma shape of per-person load multiplier (mean 1); "
                   "'inf' = no between-person variability.",
    "decay_fraction": "Fraction of shed load surviving to the sampler (0-1).",
    "transport_mean": "Mean (days) of the sewer transport delay.",
    "transport_sd": "SD (days) of the sewer transport delay.",
    "flow_mean": "Mean wastewater flow (litres/day).",
    "flow_cv": "Coefficient of variation of daily flow (0 = constant).",
    "meas_cv": "CV of multiplicative lognormal lab noise (0 = none).",
    "n_replicates": "Number of technical replicates per sampling day.",
    "lod": "Limit of detection (copies/litre); lower values -> NaN.",
    "sample_days_every": "Sampling interval in days (1 = daily).",
}


def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        description="Simulate wastewater measurements from a full phylogeny.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    io = ap.add_argument_group("input / output")
    io.add_argument("--tree", help="Path to the time-scaled tree file.")
    io.add_argument("--tree-format", default="nexus",
                    help="Tree format understood by Bio.Phylo "
                         "(newick, nexus, phyloxml, ...).")
    io.add_argument("--output_dir", default="episewer",
                    help="Path to output directory.")
    io.add_argument("--episewer-prefix", default="episewer",
                    help="If given, also write <prefix>_measurements.csv and "
                         "<prefix>_flows.csv for EpiSewer.")

    tr = ap.add_argument_group("tree -> infections")
    tr.add_argument("--mode", choices=["transmission_nodes", "tip_parent"],
                    default="transmission_nodes",
                    help="How infections are read from the tree.")
    tr.add_argument("--no-index-case", dest="include_index_case",
                    action="store_false",
                    help="Do not add an index-case infection at the root "
                         "(transmission_nodes mode only).")
    tr.add_argument("--branch-length-unit-days", type=float, default=365.25,
                    help="Days per branch-length unit (365.25 = years, "
                         "30.44 = months, 7 = weeks, 1 = days).")
    tr.add_argument("--root-day", type=float, default=0.0,
                    help="Simulation day assigned to the tree root.")

    tm = ap.add_argument_group("time axis")
    tm.add_argument("--start-date", default="2002-01-01",
                    help="Calendar date of simulation day 0. This is the date of the start of the simulation. As the simulation represents the sequencing samples, the wastewater data might start earlier than this depending on the sequencing_delay parameter.")
    tm.add_argument("--sequencing-delay", type=int, default=0,
                    help="Delay of the sequencing data in days. (Default: 0)")
    tm.add_argument("--n-days", type=int, default=None,
                    help="Length of the simulation in days. Default: tree "
                         "span + shed_max_days + 10, so the shedding tail "
                         "is fully covered.")
    tm.add_argument("--seed", type=int, default=None,
                    help="Random seed for reproducibility.")

    g = ap.add_argument_group("wastewater model parameters")
    for f in dataclasses.fields(wu.WastewaterParams):
        g.add_argument(f"--{f.name.replace('_', '-')}", type=f.type,
                       default=f.default, help=PARAM_HELP.get(f.name, ""))
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    params = wu.WastewaterParams(**{f.name: getattr(args, f.name)
                                    for f in dataclasses.fields(
                                        wu.WastewaterParams)})

    trees = dendropy.TreeList.get(path=args.tree, schema="nexus")
    
    for i, tree in enumerate(trees):
        phylo_tree = Phylo.read(StringIO(tree.as_string(schema="newick")), "newick")
        
        times = wu.infection_times_from_tree(
            phylo_tree,
            root_day=args.root_day,
            mode=args.mode,
            include_index_case=args.include_index_case,
            fmt=args.tree_format,
            branch_length_unit_days=args.branch_length_unit_days,
        )

        n_days = args.n_days or (math.ceil(times.max()) + params.shed_max_days + 10)
        inc = wu.daily_incidence(times, n_days, start_day=0)

        if inc.sum() < len(times):
            print(f"Warning: {int(len(times) - inc.sum())} of {len(times)} "
                f"infections fall outside the {n_days}-day window and were "
                f"dropped (increase --n-days or check --root-day).",
                file=sys.stderr)

        start_date_minus_delay = pd.Timestamp(args.start_date) - pd.Timedelta(days=args.sequencing_delay)
        df = wu.simulate_wastewater(inc,  
                                    start_date=start_date_minus_delay, 
                                    p=params, 
                                    seed=args.seed)
        
        tree_output_dir = output_dir / f"tree_{i}"
        tree_output_dir.mkdir(parents=True, exist_ok=True)
        
        df.to_csv(tree_output_dir / "full_data.csv", index=False)
        print(f"{len(times)} infections over {n_days} days -> "
            f"{len(df)} sampling days written to {args.output_dir}/tree_{i}/full_data.csv")

        wu.export_for_episewer(df, tree_output_dir)
        print(f"Files for wastewater analysis written to {args.output_dir}/tree_{i}")
            
    return 0


if __name__ == "__main__":
    sys.exit(main())