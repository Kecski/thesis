#!/usr/bin/env python3
"""
Command-line interface: simulate wastewater measurements from a full,
time-scaled phylogeny using parameters loaded from a TOML configuration file.

Example
-------
python wastewater_data.py --config config_wastewater.toml
"""

import argparse
import dataclasses
import math
from io import StringIO
from pathlib import Path
import sys
import pandas as pd
from Bio import Phylo
import dendropy
import tomllib
from general_utils import log_and_save_parameters
import wastewater_data_utils as wu


def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate wastewater measurements from a full phylogeny using a TOML configuration file."
    )
    p.add_argument(
        "--config",
        "-c",
        type=str,
        default="config_wastewater.toml",
        help="Path to the TOML configuration file (default: config_wastewater.toml).",
    )
    p.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="data/episewer",
        help="Output directory for the generated wastewater data and logs (default: 'data/episewer').",
    )
    p.add_argument(
        "--tree",
        "-t",
        type=str,
        required=True,
        help="Path to the simulated phylogeny file (Newick or Nexus format)",
    )
    p.add_argument(
        "--tree_format",
        "-f",
        type=str,
        default="nexus",
        help="Format of the input tree file (default: 'nexus'). Options: 'newick', 'nexus'"
    )
    p.add_argument(
        "--tree_id",
        "-i",
        type=int,
        required=True,
        help="ID ([0, N-1]) of the tree in the trees file containing N trees."
    )
    p.add_argument(
        "--sequencing_delay",
        "-s",
        type=int,
        default=0,
        help="Delay in the reporting of the clinically sequenced cases (in days) (default: 0)."
    )
    p.add_argument(
        "--n_days",
        "-n",
        type=int,
        default=0,
        help="Number of days to simulate wastewater measurements (default: 0). If 0, the simulation will run for the full duration of the phylogeny.",
    )
    return p.parse_args()


def load_config(config_path: Path) -> dict:
    if not config_path.exists():
        sys.exit(f"Error: Configuration file not found at '{config_path}'.")

    with open(config_path, "rb") as f:
        raw_cfg = tomllib.load(f)

    # if "tree" not in raw_cfg:
    #     sys.exit(f"Error: Missing required key 'tree' in '{config_path}'.")

    # Target default types from WastewaterParams dataclass
    dataclass_types = {f.name: f.type for f in dataclasses.fields(wu.WastewaterParams)}

    cfg = {
        "mode": str(raw_cfg.get("mode", "transmission_nodes")),
        "include_index_case": bool(raw_cfg.get("include_index_case", True)),
        "branch_length_unit_days": float(raw_cfg.get("branch_length_unit_days", 365.25)),
        "root_day": float(raw_cfg.get("root_day", 0.0)),
        "start_date": str(raw_cfg.get("start_date", "2002-01-01")),
        "seed": raw_cfg.get("seed"),
    }

    # Load and cast model parameters expected by wu.WastewaterParams
    for f in dataclasses.fields(wu.WastewaterParams):
        if f.name in raw_cfg:
            val = raw_cfg[f.name]
            
            # Special case for handling indiv_shape which accepts "inf" or math.inf
            if f.name == "indiv_shape":
                if str(val).lower() == "inf":
                    val = math.inf
                else:
                    val = float(val)
            # Handle float conversions from TOML numbers
            if f.type == float and val is not None:
                val = float(val)
            elif f.type == int and val is not None:
                val = int(val)
            cfg[f.name] = val
        else:
            cfg[f.name] = f.default

    return cfg


def validate_config(cfg: dict):
    errors = []

    if cfg["mode"] not in ["transmission_nodes", "tip_parent"]:
        errors.append("'mode' must be either 'transmission_nodes' or 'tip_parent'.")

    if cfg["branch_length_unit_days"] <= 0:
        errors.append("'branch_length_unit_days' must be > 0.")

    if cfg["sequencing_delay"] < 0:
        errors.append("'sequencing_delay' must be >= 0.")

    if cfg["n_days"] < 0:
        errors.append("'n_days' must be >= 0.")

    if errors:
        sys.exit("Validation Error(s) in configuration:\n  - " + "\n  - ".join(errors))


def main(argv=None) -> int:
    args = parse_args()
    config_path = Path(args.config)

    cfg = load_config(config_path)
    cfg["output_dir"] = args.output_dir
    cfg["tree"] = args.tree
    cfg["tree_format"] = args.tree_format
    cfg["tree_id"] = args.tree_id
    cfg["sequencing_delay"] = args.sequencing_delay
    cfg["n_days"] = args.n_days
    
    validate_config(cfg)
    
    tree_id = cfg["tree_id"]

    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    # Initialize model parameter dataclass directly from parsed configuration
    param_kwargs = {f.name: cfg[f.name] for f in dataclasses.fields(wu.WastewaterParams)}
    params = wu.WastewaterParams(**param_kwargs)

    db_trees = dendropy.TreeList.get(path=cfg["tree"], schema=cfg["tree_format"])

    print(f"Loaded config from: {config_path}")
    log_and_save_parameters(cfg, output_dir, filename="parameters.txt")

    phylo_tree = Phylo.read(StringIO(db_trees[tree_id].as_string(schema="newick")), "newick")

    times = wu.infection_times_from_tree(
        phylo_tree,
        root_day=cfg["root_day"],
        mode=cfg["mode"],
        include_index_case=cfg["include_index_case"],
        fmt=cfg["tree_format"],
        branch_length_unit_days=cfg["branch_length_unit_days"],
    )

    n_days = cfg["n_days"] + cfg["sequencing_delay"]
    inc = wu.daily_incidence(times, n_days, start_day=0)

    if inc.sum() < len(times):
        print(
            f"Warning: {int(len(times) - inc.sum())} of {len(times)} "
            f"infections fall outside the {n_days}-day window and were "
            f"dropped (increase n_days or check root_day).",
            file=sys.stderr,
        )

    start_date_minus_delay = pd.Timestamp(cfg["start_date"]) - pd.Timedelta(
        days=cfg["sequencing_delay"]
    )
    df = wu.simulate_wastewater(
        inc, start_date=start_date_minus_delay, p=params, seed=cfg["seed"]
    )

    tree_output_dir = output_dir / f"tree_{tree_id}"
    tree_output_dir.mkdir(parents=True, exist_ok=True)

    df.to_csv(tree_output_dir / "full_data.csv", index=False)
    print(
        f"{len(times)} infections over {n_days} days -> "
        f"{len(df)} sampling days written to {tree_output_dir / 'full_data.csv'}"
    )

    wu.export_for_episewer(df, tree_output_dir)
    print(f"Files for wastewater analysis written to {tree_output_dir}")

    return 0


if __name__ == "__main__":
    sys.exit(main())