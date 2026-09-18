import dendropy
import pyvolve
from Bio import SeqIO

# def remove_coma_in_nexus_file(input_file: str, output_file: str) -> None:
#     with open(input_file, "r") as f:
#         content = f.read()

#     # Remove trailing comma before a semicolon (handles whitespace/newlines between them)
#     content_fixed = re.sub(r',(\s*;)', r'\1', content)

#     with open(output_file, "w") as f:
#         f.write(content_fixed)
        
# def convert_nexus_to_newick(nexus_file: str, newick_file: str) -> None:
#     trees = dendropy.TreeList.get(path=nexus_file, schema="nexus")
#     trees.write(path=newick_file, schema="newick")

def add_dates(tree, input_fasta: str, root_date: float=2002.0) -> None:
    """
    Calculates the tip dates based on root-to-tip distances and appends them to FASTA headers.
    
    Args:
        tree (dendropy.Tree): The phylogenetic tree with branch lengths.
        input_fasta (str): Path to the FASTA file containing sequences for the tree tips.
        root_date (float): The base chronological date of the root node. Defaults to 2002.0.
        
    Returns:
        None. The input FASTA file is overwritten in place with updated sequence IDs.
    """
    tip_times = {}
    for leaf in tree.leaf_node_iter():
        dist = leaf.distance_from_root()
        label = "_".join(leaf.taxon.label.split(" "))
        tip_times[label] = dist
        
    sampling_dates = {name: root_date + t for name, t in tip_times.items()}
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    for rec in records:
        date = sampling_dates[rec.id]
        rec.id = f"{rec.id}_{date:.4f}"
        rec.description = ""

    # output_fasta = f"{input_fasta}"
    SeqIO.write(records, input_fasta, "fasta")

def simulate_sequences_on_tree(tree_file: str,
                               root_sequence: str,
                               clock_rate: float=1.,
                               frequencies: list[float]=[0.25, 0.25, 0.25, 0.25],
                               mutation_rates: dict={"AC":1.0, "AG":1.0, "AT":1.0, "CG":1.0, "CT":1.0, "GT":1.0},
                               gamma_dist: dict={"alpha": 0.5, "num_categories": 4},
                               output_dir: str="./",
                               root_date: float=2002.0
                               ) -> None:
    """
    Simulates nucleotide sequence evolution along trees in a Nexus file and saves the output.
    
    Args:
        tree_file (str): Path to the input Nexus tree file.
        root_sequence (str): The ancestral nucleotide sequence to evolve.
        clock_rate (float): Branch length multiplier to scale the evolutionary rate. Defaults to 1.0.
        frequencies (list[float]): Equilibrium frequencies of nucleotides [A, C, G, T].
        mutation_rates (dict): Relative substitution rates for nucleotide pairs.
        gamma_dist (dict): Parameters for Gamma rate heterogeneity (alpha and number of categories).
        output_dir (str): Destination directory for the generated fasta, info, and rate files.
        root_date (float): Chronological date of the root, passed to add_dates(). Defaults to 2002.0.
        
    Returns:
        None. Generates sequence, info, and rate files in the output directory.
    """
    model = pyvolve.Model("nucleotide", 
                          {"mu": mutation_rates, "state_freqs": frequencies},
                          alpha=gamma_dist["alpha"],
                          num_categories=gamma_dist["num_categories"])
    partition = pyvolve.Partition(models=model, root_sequence=root_sequence)
        
    trees = dendropy.TreeList.get(path=tree_file, schema="nexus")
    for i, dp_tree in enumerate(trees):
        tree = pyvolve.read_tree(tree=dp_tree.as_string(schema="newick"), scale_tree=clock_rate)
        evolver = pyvolve.Evolver(partitions=partition, tree=tree)
        
        # Output files for each tree
        fasta_file = f"{output_dir}/sequences_{i}.fasta"
        info_file = f"{output_dir}/sequences_{i}_info.txt"
        rates_file = f"{output_dir}/sequences_{i}_rates.txt"

        evolver(seqfile=fasta_file, seqfmt="fasta", 
                infofile=info_file, 
                ratefile=rates_file)
        
        add_dates(dp_tree, fasta_file, root_date)