"""Sequence Mutation Simulator along Phylogenetic Trees using TOML configuration.

This module provides a command-line utility to simulate nucleotide sequence evolution
along a specified phylogenetic tree using configuration parameters loaded from a TOML file.

Example:
    $ python sequence_mutation.py --config config_mutation.toml
"""

import argparse
from pathlib import Path
import random
import tomllib
import sys
from Bio import SeqIO
from general_utils import log_and_save_parameters
from sequence_mutation_utils import simulate_sequences_on_tree


def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate sequence mutations along a given phylogeny using a TOML configuration file."
    )
    p.add_argument(
        "--config",
        "-c",
        type=str,
        default="config_mutation.toml",
        help="Path to the TOML configuration file (default: config_mutation.toml).",
    )
    p.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="data/simulated_sequences",
        help="Output directory for the generated sequences and logs (default: 'data/simulated_sequences').",
    )
    p.add_argument(
        "--root_sequence_file",
        "-r",
        type=str,
        required=True,
        help="Path to the root sequence file (default: 'data/root_sequences/test_root_sequence_500.txt').",
    )
    p.add_argument(
        "--tree_file",
        "-t",
        type=str,
        required=True,
        help="Path to the phylogenetic tree file (default: 'data/remaster_trees/logs/sim.sampled.trees').",
    )
    p.add_argument(
        "--tree_id",
        "-i",
        type=int,
        required=True,
        help="ID ([0, N-1]) of the tree in the trees file containing N trees."
    )
    return p.parse_args()


def load_config(config_path: Path) -> dict:
    if not config_path.exists():
        sys.exit(f"Error: Configuration file not found at '{config_path}'.")

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    # required_keys = ["tree_file", "root_sequence_file"]
    # missing = [k for k in required_keys if k not in config]
    # if missing:
    #     sys.exit(f"Error: Missing required keys in '{config_path}': {', '.join(missing)}")

    return {
        "sub_model": str(config.get("sub_model", "HKY")),
        "kappa": float(config.get("kappa", 2.0)),
        "transition_rates": config.get("transition_rates"),
        "frequencies": [float(f) for f in config.get("frequencies", [0.25, 0.25, 0.25, 0.25])],
        "clock_rate": float(config.get("clock_rate", 1.0)),
        "gamma_shape": float(config.get("gamma_shape", 0.5)),
        "gamma_categories": int(config.get("gamma_categories", 4)),
        "root_date": float(config.get("root_date", 2002.0)),
        "seed": config.get("seed"),
    }


def validate_config(cfg: dict):
    errors = []

    if cfg["sub_model"] not in ["HKY", "GTR"]:
        errors.append("'sub_model' must be either 'HKY' or 'GTR'.")

    # Use round() to avoid standard floating-point precision error (e.g. 0.30 + 0.18 + 0.20 + 0.32 == 1.0000000000000002)
    if round(sum(cfg["frequencies"]), 5) != 1.0:
        errors.append(f"'frequencies' must sum to 1.0 (currently sums to {sum(cfg['frequencies'])}).")

    if cfg["transition_rates"] is not None and len(cfg["transition_rates"]) != 6:
        errors.append("'transition_rates' must contain exactly 6 values (AC, AG, AT, CG, CT, GT).")

    if errors:
        sys.exit("Validation Error(s) in configuration:\n  - " + "\n  - ".join(errors))


def main():
    cli_args = parse_args()
    config_path = Path(cli_args.config)

    cfg = load_config(config_path)
    cfg["output_dir"] = cli_args.output_dir
    cfg["root_sequence_file"] = cli_args.root_sequence_file
    cfg["tree_file"] = cli_args.tree_file
    cfg["tree_id"] = cli_args.tree_id
    validate_config(cfg)

    if cfg["seed"] is not None:
        random.seed(cfg["seed"])

    if cfg["transition_rates"]:
        custom_mu = {
            "AC": cfg["transition_rates"][0],
            "AG": cfg["transition_rates"][1],
            "AT": cfg["transition_rates"][2],
            "CG": cfg["transition_rates"][3],
            "CT": cfg["transition_rates"][4],
            "GT": cfg["transition_rates"][5],
        }
    else:
        custom_mu = None

    gamma_dist = {
        "alpha": cfg["gamma_shape"],
        "num_categories": cfg["gamma_categories"],
    }

    root_seq_path = Path(cfg["root_sequence_file"])
    root_file_type = root_seq_path.suffix.lstrip(".").lower()

    if root_file_type == "txt":
        with open(root_seq_path, "r") as f:
            root_sequence = f.read().replace("\n", "").replace("\r", "").strip()
    elif root_file_type in ["fasta", "fa"]:
        root_sequence = next(SeqIO.parse(root_seq_path, "fasta")).seq
        root_sequence = str(root_sequence).replace("N", "").replace("-", "")
    else:
        sys.exit(f"Error: Sequence root file extension must be '.txt' or '.fasta' (got: '.{root_file_type}')")

    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    log_and_save_parameters(cfg, output_dir, filename="parameters.txt")

    simulate_sequences_on_tree(
        tree_file=cfg["tree_file"],
        tree_id=cfg["tree_id"],
        root_sequence=root_sequence,
        sub_model=cfg["sub_model"],
        clock_rate=cfg["clock_rate"],
        frequencies=cfg["frequencies"],
        kappa=cfg["kappa"],
        mutation_rates=custom_mu,
        gamma_dist=gamma_dist,
        output_dir=cfg["output_dir"],
        root_date=cfg["root_date"],
    )


if __name__ == "__main__":
    main()