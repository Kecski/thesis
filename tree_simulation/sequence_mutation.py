"""Sequence Mutation Simulator along Phylogenetic Trees.

This module provides a command-line utility to simulate nucleotide sequence evolution
along a specified phylogenetic tree using custom evolutionary parameters (transition rates,
equilibrium base frequencies, rate heterogeneity across sites, and molecular clock rates).

Example:
    Basic execution:
        $ python simulate_mutations.py --tree_file tree.nexus --root_sequence_file root.txt

    Advanced execution:
        $ python simulate_mutations.py \\
            --tree_file tree.nexus \\
            --root_sequence_file root.fasta \\
            --root_file_type fasta \\
            --transition_rates 1 2 1 1 2 1 \\
            --frequencies 0.3 0.2 0.2 0.3 \\
            --clock_rate 0.005 \\
            --output_dir results \\
            --seed 42
"""

import random
import argparse
from pathlib import Path
from sequence_utils import simulate_sequences_on_tree
from Bio import SeqIO

def parse_args():
    """Parse command-line arguments for the mutation simulator.

    Returns:
        argparse.Namespace: Object containing all parsed command-line arguments:
            - tree_file (str): Path to input NEXUS phylogenetic tree.
            - root_sequence_file (str): Path to root sequence file.
            - root_file_type (str): Format of root file ('txt', 'fasta', etc.).
            - transition_rates (list[float]): Relative nucleotide transition rates (6 values).
            - frequencies (list[float]): Equilibrium base frequencies (A, C, G, T).
            - clock_rate (float): Substitutions per site per unit time.
            - gamma_shape (float): Shape parameter (alpha) for gamma site rate variation.
            - gamma_categories (int): Number of rate categories for discrete gamma model.
            - root_date (float): Temporal date (in years) associated with tree root.
            - output_dir (str): Destination directory for outputs and parameter logs.
            - seed (float | None): Seed for pseudorandom number generation.
    """
    p = argparse.ArgumentParser(
        description="Simulate sequence mutations along a given phylogeny"
    )
    # Required arguments
    p.add_argument("--tree_file", type=str, required=True,
                   help="Path to the defined phylogeny/phylogenies in NEXUS format.")
    p.add_argument("--root_sequence_file", type=str, required=True,
                   help="Path to the file containig the sequence of the root of the phylogeny.")
    
    # Optional arguments
    # p.add_argument("--root_file_type", type=str, default="txt",
    #                help="Extension of the root sequence file.")
    p.add_argument("--sub_model", type=str, default="HKY",
                   help="Type of substitution model used for the simulation. The possible options are HKY and GTR.")
    p.add_argument("--kappa", type=float, default=2.0,
                   help="Ti/Tv rate for HKY substitution model.")
    p.add_argument("--transition_rates", nargs=6, type=float, default=None,
                   help="Relative transition probabilities between the nucleotides for GTR model. Provide the values in the following order: AC, AG, AT, CG, CT, GT.")
    p.add_argument("--frequencies", nargs=4, type=float, default=[0.25, 0.25, 0.25, 0.25],
                   help="Equilibrium frequencies of the nucleotides. Provide the values in the following order: A, C, G, T.")
    p.add_argument("--clock_rate", type=float, default=1.0,
                   help="The speed of mutations: substitution per site per unite time.")
    p.add_argument("--gamma_shape", type=float, default=0.5,
                   help="Shape parameter of the gamma distribution for the discrete gamma site model.")
    p.add_argument("--gamma_categories", type=int, default=4,
                   help="Number of categories in the gamma site model.")
    p.add_argument("--root_date", type=float, default=2002.0,
                   help="Date in years of the root of the phylogeny.")
    p.add_argument("--output_dir", type=str, default="sequences",
                   help="Directory to which the output files will be saved.")
    p.add_argument("--seed", type=float, default=None,
                   help="Seed for the simulation")
    
    return p.parse_args()

def validate_args(args, parser):
    """Validate argument values and consistency before simulation execution.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
        parser (argparse.ArgumentParser): Parser instance used to emit CLI errors.

    Raises:
        SystemExit: If nucleotide equilibrium frequencies do not sum to 1.0.
    """
    if args.sub_model not in ["HKY", "GTR"]:
        parser.error("The substitution model provided for --sub_model has to be either HKY or GTR.")
    if sum(args.frequencies) != 1:
        parser.error("The values provided to --frequencies has to add up to 1.")
    

def main():
    """Coordinates simulation workflow: parses arguments, loads sequences,

    writes run configurations, and calls external tree simulation utilities.
    """
    args = parse_args()
    validate_args(args, argparse.ArgumentParser())
    
    if args.seed is not None:
        random.seed(args.seed)
    
    if args.transition_rates:
        custom_mu = {"AC":args.transition_rates[0], 
                    "AG":args.transition_rates[1], 
                    "AT":args.transition_rates[2],
                    "CG":args.transition_rates[3],
                    "CT":args.transition_rates[4],
                    "GT":args.transition_rates[5]}
    else:
        custom_mu = None
        
    gamma_dist = {"alpha": args.gamma_shape, "num_categories": args.gamma_categories}
    
    root_file_type = args.root_sequence_file.split(".")[-1]
    if root_file_type == "txt":
        with open(args.root_sequence_file, "r") as f:
            root_sequence = f.read().replace("\n", "").replace("\r", "").strip()
    elif root_file_type == "fasta":
        root_sequence = next(SeqIO.parse(args.root_sequence_file, args.root_file_type)).seq
        root_sequence = root_sequence.replace("N", "")
        root_sequence = str(root_sequence.replace("-", ""))
    else:
        raise ValueError(f"Sequence root file has to be either txt or fasta. Not {root_file_type}")
    
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print("Simulating sequence mutations along phylogeny:")
    print(f"\tTree file: {args.tree_file}")
    print(f"\tRoot sequence file: {args.root_sequence_file}")
    print(f"\tSubstitution model: {args.sub_model}")
    print(f"\tKappa: {args.kappa}")
    print(f"\tTransition rates: {custom_mu}")
    print(f"\tEquilibrium frequences: {args.frequencies}")
    print(f"\tClock rate: {args.clock_rate}")
    print(f"\tGamma site model: {gamma_dist}")
    print(f"\tRoot date: {args.root_date}")
    print(f"\tOutput directory: {args.output_dir}")
    print(f"\tSeed: {args.seed}")
    
    param_log = output_dir / "parameters.txt"
    with open(param_log, "w") as f:
        f.write(f"Tree file: {args.tree_file}\n")
        f.write(f"Root sequence file: {args.root_sequence_file}\n")
        f.write(f"Substitution model: {args.sub_model}\n")
        f.write(f"Kappa: {args.kappa}")
        f.write(f"Transition rates: {custom_mu}\n")
        f.write(f"Equilibrium frequences: {args.frequencies}\n")
        f.write(f"Clock rate: {args.clock_rate}\n")
        f.write(f"Gamma site model: {gamma_dist}\n")
        f.write(f"Root date: {args.root_date}\n")
        f.write(f"Output directory: {args.output_dir}\n")
        if args.seed is not None:
            f.write(f"Seed: {args.seed}\n")
    
    
    # remove_coma_in_nexus_file("./remaster_sims/sim_logs/sim.full.trees", "./remaster_sims/sim_logs/sim.full_fixed.trees")
    
    # custom_mu = {"AC":1, "AG":1, "AT":1, "CG":1, "CT":1, "GT":1}
    # freqs = [0.1, 0.45, 0.3, 0.15] # A, C, G, T
    # gamma_dist = {"alpha": 0.5, "num_categories": 4}
    # clock_rate = 1 # Substitutions per site per unit time
    # root_sequence = "ACGT"*int(1e3)  # Root sequence of the tree
    # root_date = 2002.0
    
    simulate_sequences_on_tree(tree_file=args.tree_file,
                               root_sequence=root_sequence,
                               sub_model=args.sub_model,
                               clock_rate=args.clock_rate,
                               frequencies=args.frequencies,
                               kappa=args.kappa,
                               mutation_rates=custom_mu,
                               gamma_dist=gamma_dist,
                               output_dir=args.output_dir,
                               root_date=args.root_date)
        
    

if __name__ == "__main__":
    main()
