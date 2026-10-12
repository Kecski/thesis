"""
Command-line entry point for simulating stochastic birth-death trees with
BEAST2 + ReMASTER.
 
This script parses command-line arguments describing a birth-death-sampling
process (rates, optional piecewise rate-change times, population/tree
stopping conditions, etc.), validates them, writes a ReMASTER XML
configuration file via `write_single_tree_xml` (from `xml_utils`), and logs
the chosen parameters to a text file alongside the generated XML.
 
Run directly, e.g.:
    python run_remaster_sim.py --birth_rate 1.4 --death_rate 0.5 --max_time 10
"""

import argparse
import random
import sys
import tomllib
from pathlib import Path
from general_utils import log_and_save_parameters
from create_remaster_xml_utils import write_single_tree_xml, write_single_tree_for_wastewater_xml

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate stochastic birth-death trees with BEAST2 + ReMASTER using a TOML config file."
    )
    p.add_argument(
        "--config",
        "-c",
        type=str,
        default="config.toml",
        help="Path to the TOML configuration file (default: config.toml).",
    )
    p.add_argument(
        "--max_time",
        required=True,
        type=float,
        help="Maximum time for the clinical sequencing during the simulation (in years). The simulation will run for some additional days equal to the sequencing_delay parameter for wastewater simulation purposes."
    )
    p.add_argument(
        "--output_dir",
        "-o",
        type=str,
        default="data/remaster_trees",
        help="Output directory for the generated XML and logs (default: 'data/remaster_trees')."
        )
    p.add_argument(
        "--output_name",
        "-n",
        type=str,
        default="sim",
        help="Base name for the generated XML and logs (default: 'sim')."
    )
    p.add_argument(
        "--n_sims",
        type=int,
        default=1,
        help="Number of independent phylogenies simulated by the xml file with the same parameters (default: 1)."
    )
    p.add_argument(
        "--sequencing_delay",
        type=int,
        default=0,
        help="Delay in the reporting of the clinically sequenced cases (in days) (default: 0)."
    )
    return p.parse_args()

def load_config(config_path: Path) -> dict:
    if not config_path.exists():
        sys.exit(f"Error: Configuration file not found at '{config_path}'.")

    with open(config_path, "rb") as f:
        config = tomllib.load(f)

    # Helper function to convert single values or lists to list of floats
    def to_float_list(val):
        if val is None:
            return None
        if isinstance(val, (int, float)):
            return [float(val)]
        return [float(x) for x in val]

    # Required fields
    required = ["birth_rate", "death_rate"]
    missing = [k for k in required if k not in config]
    if missing:
        sys.exit(f"Error: Missing required keys in '{config_path}': {', '.join(missing)}")

    # Parsed config dictionary with default fallback values
    return {
        "birth_rate": to_float_list(config["birth_rate"]),
        "death_rate": to_float_list(config["death_rate"]),
        "sampling_rate": to_float_list(config.get("sampling_rate", [1.0])),
        "change_times_birth": config.get("change_times_birth"),
        "change_times_death": config.get("change_times_death"),
        "change_times_sampling": config.get("change_times_sampling"),
        "n_sims": int(config.get("n_sims", 1)),
        "init_pop": int(config.get("init_pop", 1)),
        "traj_type": str(config.get("traj_type", "StochasticTrajectory")),
        "max_pop": config.get("max_pop"),
        "min_tips": config.get("min_tips"),
        "beast_cmd": str(config.get("beast_cmd", "beast")),
        "seed": config.get("seed"),
    }

def validate_config(cfg: dict):
    errors = []

    for name in ["birth_rate", "death_rate", "sampling_rate"]:
        if any(r < 0 for r in cfg[name]):
            errors.append(f"'{name}' values must all be >= 0.")

    if cfg["max_time"] <= 0:
        errors.append("'max_time' must be > 0.")
    if cfg["sequencing_delay"] < 0:
        errors.append("'sequencing_delay' must be >= 0.")
    if cfg["n_sims"] < 1:
        errors.append("'n_sims' must be >= 1.")
    if cfg["init_pop"] < 1:
        errors.append("'init_pop' must be >= 1.")
    if cfg["max_pop"] is not None and cfg["max_pop"] < 1:
        errors.append("'max_pop' must be >= 1.")
    if cfg["min_tips"] is not None and cfg["min_tips"] < 1:
        errors.append("'min_tips' must be >= 1.")

    for ct_key in ["change_times_birth", "change_times_death", "change_times_sampling"]:
        cts = cfg[ct_key]
        if cts:
            is_ascending = all(cts[i] < cts[i + 1] for i in range(len(cts) - 1))
            if not is_ascending:
                errors.append(f"'{ct_key}' must be in ascending order.")
            if cts[-1] >= cfg["max_time"]:
                errors.append(f"'{ct_key}' values must be less than 'max_time'.")

    valid_traj_types = ["StochasticTrajectory", "DeterministicTrajectory"]
    if cfg["traj_type"] not in valid_traj_types:
        errors.append(f"'traj_type' must be one of: {', '.join(valid_traj_types)}.")

    if errors:
        sys.exit("Validation Error(s) in configuration:\n  - " + "\n  - ".join(errors))

def main():
    cli_args = parse_args()
    config_path = Path(cli_args.config)
    
    cfg = load_config(config_path)
    cfg["output_dir"] = cli_args.output_dir
    cfg["output_name"] = cli_args.output_name
    cfg["max_time"] = cli_args.max_time
    cfg["sequencing_delay"] = cli_args.sequencing_delay
    cfg["n_sims"] = cli_args.n_sims
    
    validate_config(cfg)

    # Set random seed if specified
    if cfg["seed"] is not None:
        random.seed(cfg["seed"])

    print(f"Generating single-tree XML file using config: {config_path}")
    for key, val in cfg.items():
        print(f"\t{key}: {val}")

    # Output directory setup
    output_dir = Path(cfg["output_dir"])
    output_dir.mkdir(parents=True, exist_ok=True)

    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)

    # Parameter record log
    log_and_save_parameters(cfg, output_dir, filename="parameters.txt")

    # Generate XML file
    xml_file_path = output_dir / f"{cfg['output_name']}.xml"
    write_single_tree_for_wastewater_xml(
        file_path=str(xml_file_path),
        n_sims=cfg["n_sims"],
        birth_rate=cfg["birth_rate"],
        death_rate=cfg["death_rate"],
        change_times_birth=cfg["change_times_birth"],
        change_times_death=cfg["change_times_death"],
        change_times_sampling=cfg["change_times_sampling"],
        max_time=cfg["max_time"],
        sequencing_delay=cfg["sequencing_delay"],
        sampling_rate=cfg["sampling_rate"],
        loggers_base=str(Path("logs") / cfg["output_name"]),
        traj_type=cfg["traj_type"],
        init_pop=cfg["init_pop"],
        ends_when=f"X>={cfg['max_pop']}" if cfg["max_pop"] is not None else None,
        must_have=f"S>={cfg['min_tips']}" if cfg["min_tips"] is not None else None,
    )


if __name__ == "__main__":
    main()