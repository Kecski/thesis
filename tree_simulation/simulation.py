import argparse
import random
import re
import shutil
import subprocess
import sys
from pathlib import Path
from xml_utils import write_single_tree_xml

def parse_args():
    p = argparse.ArgumentParser(
        description="Simulate stochastic birth-death trees with BEAST2 + ReMASTER."
    )
    # Required arguments
    p.add_argument("--birth_rate", type=float, required=True,
                   help="Per-lineage birth rate (lambda) for the birth-death process.")
    p.add_argument("--death_rate", type=float, required=True,
                   help="Per-lineage death rate (mu) for the birth-death process.")
    p.add_argument("--max_time", type=float, required=True,
                   help="Maximum simulation time.")

    # Optional arguments
    p.add_argument("--sampling_rate", type=float, default=1.0,
                    help="Per-lineage sampling rate through time. 0 = no through-time sampling reaction.")
    p.add_argument("--n_sims", type=int, default=1,
                    help="Number of independent replicate simulations to run.")
    p.add_argument("--init_pop", type=int, default=1,
                    help="Initial size of population S.")
    p.add_argument("--traj_type", choices=["StochasticTrajectory", "DeterministicTrajectory"], default="StochasticTrajectory",
                    help="'tree' logs only the sampled tree; 'full' also logs the population trajectory.")
    p.add_argument("--output_dir", type=str, default="./remaster_sims",
                    help="Directory for generated XML and output files.")
    p.add_argument("--output_name", type=str, default="sim",
                    help="Base filename used for XML/tree/trajectory outputs.")
    p.add_argument("--max_pop", type=int, default=None,
                    help="Optional cap on population size; simulation stops early if reached.")
    p.add_argument("--min_tips", type=int, default=None,
                    help="Minimum sampled tips required; re-simulates with a new seed until met.")
    p.add_argument("--beast_cmd", type=str, default="beast",
                    help="Command/path used to invoke BEAST2.")
    p.add_argument("--seed", type=int, default=None,
                    help="Base random seed for reproducibility of the seed sequence.")
    return p.parse_args()

def validate_args(args, parser):
    """Required fields already enforced by argparse; here we check value ranges / allowed choices, 
    raising a clean error (no traceback) via parser.error(), same spirit as stop(call.=FALSE)."""
 
    if args.birth_rate < 0:
        parser.error("--birth_rate must be >= 0.")
    if args.death_rate < 0:
        parser.error("--death_rate must be >= 0.")
    if args.sampling_rate < 0:
        parser.error("--sampling_rate must be >= 0.")
    if args.max_time <= 0:
        parser.error("--max_time must be > 0.")
    if args.n_sims < 1:
        parser.error("--n_sims must be >= 1.")
    if args.init_pop < 1:
        parser.error("--init_pop must be >= 1.")
    if args.max_pop is not None and args.max_pop < 1:
        parser.error("--max_pop must be >= 1.")
    if args.min_tips is not None and args.min_tips < 1:
        parser.error("--min_tips must be >= 1.")
 
    valid_traj_types = ["StochasticTrajectory", "DeterministicTrajectory"]
    if args.traj_type not in valid_traj_types:
        parser.error("--traj_type must be one of: " + ", ".join(valid_traj_types) + ".")
        

def main():
    args = parse_args()
    validate_args(args, parser=argparse.ArgumentParser())
    
    # Set random seed for reproducibility
    if args.seed is not None:
        random.seed(args.seed)
    
    # Print simulation parameters
    print("Generating single-tree XML file for simulation...")
    print("\tbirth_rate:", args.birth_rate)
    print("\tdeath_rate:", args.death_rate)
    print("\tsampling_rate:", args.sampling_rate)
    print("\tn_sims:", args.n_sims)
    print("\tmax_time:", args.max_time)
    print("\tinit_pop:", args.init_pop)
    print("\ttraj_type:", args.traj_type)
    print("\toutput_dir:", args.output_dir)
    print("\toutput_name:", args.output_name)
    print("\tmax_pop:", args.max_pop)
    print("\tmin_tips:", args.min_tips)
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Write parameters to txt file for record-keeping
    param_log = output_dir / f"{args.output_name}_parameters.txt"
    with open(param_log, "w") as f:
        f.write(f"birth_rate: {args.birth_rate}\n")
        f.write(f"death_rate: {args.death_rate}\n")
        f.write(f"sampling_rate: {args.sampling_rate}\n")
        f.write(f"n_sims: {args.n_sims}\n")
        f.write(f"max_time: {args.max_time}\n")
        f.write(f"init_pop: {args.init_pop}\n")
        f.write(f"traj_type: {args.traj_type}\n")
        f.write(f"output_dir: {args.output_dir}\n")
        f.write(f"output_name: {args.output_name}\n")
        f.write(f"max_pop: {args.max_pop}\n")
        f.write(f"min_tips: {args.min_tips}\n")
        if args.seed is not None:
            f.write(f"seed: {args.seed}\n")

    # Generate XML file for simulation
    xml_file_path = output_dir / f"{args.output_name}.xml"
    write_single_tree_xml(
        file_path=str(xml_file_path),
        n_sims=args.n_sims,
        birth_rate=args.birth_rate,
        death_rate=args.death_rate,
        max_time=args.max_time,
        sampling_rate=args.sampling_rate,
        loggers_base=args.output_name,
        traj_type=args.traj_type,
        init_pop=args.init_pop,
        ends_when=f"populationSize >= {args.max_pop}" if args.max_pop is not None else None,
        must_have=f"sampledTips >= {args.min_tips}" if args.min_tips is not None else None
    )
    
    deterministic_file_path = output_dir / f"{args.output_name}_deterministic.xml"
    write_single_tree_xml(
        file_path=str(deterministic_file_path),
        n_sims=1,
        birth_rate=args.birth_rate,
        death_rate=args.death_rate,
        max_time=args.max_time,
        sampling_rate=args.sampling_rate,
        loggers_base=args.output_name,
        traj_type="DeterministicTrajectory",
        init_pop=args.init_pop,
        ends_when=f"populationSize >= {args.max_pop}" if args.max_pop is not None else None,
        must_have=f"sampledTips >= {args.min_tips}" if args.min_tips is not None else None
    )

if __name__ == "__main__":
    main()
 
    