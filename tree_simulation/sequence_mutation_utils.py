import dendropy
import pyvolve
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from pathlib import Path
from typing import Optional
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

def _tip_sampling_dates(tree, root_date: float=2002.0) -> dict:
    """
    Returns {tip name (spaces replaced by underscores): sampling date} for every tip,
    where date = root_date + root-to-tip distance.
    """
    sampling_dates = {}
    for leaf in tree.leaf_node_iter():
        label = "_".join(leaf.taxon.label.split(" "))
        sampling_dates[label] = root_date + leaf.distance_from_root()
    return sampling_dates


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
    sampling_dates = _tip_sampling_dates(tree, root_date)
    
    records = list(SeqIO.parse(input_fasta, "fasta"))
    for rec in records:
        date = sampling_dates[rec.id]
        rec.id = f"{rec.id}_{date:.4f}"
        rec.description = ""

    # output_fasta = f"{input_fasta}"
    SeqIO.write(records, input_fasta, "fasta")


def save_tree_with_dates(tree, output_file: str, root_date: float=2002.0,
                         schema: str="newick") -> None:
    """
    Saves a copy of the tree whose tip names carry the sampling date, using exactly the
    same naming as add_dates() (e.g. "tipA_2005.1234"), so tree tips match FASTA headers.
    
    Args:
        tree (dendropy.Tree): The phylogenetic tree with branch lengths.
        output_file (str): Path of the tree file to write.
        root_date (float): The base chronological date of the root node. Defaults to 2002.0.
        schema (str): Output format, e.g. "newick" or "nexus". Defaults to "newick".
        
    Returns:
        None. The original tree object is left unchanged.
    """
    sampling_dates = _tip_sampling_dates(tree, root_date)
    
    # depth=2 also copies the taxa, so renaming them doesn't touch the original tree
    dated_tree = tree.clone(depth=2)
    for leaf in dated_tree.leaf_node_iter():
        name = "_".join(leaf.taxon.label.split(" "))
        leaf.taxon.label = f"{name}_{sampling_dates[name]:.4f}"
    
    Path(output_file).parent.mkdir(parents=True, exist_ok=True)
    # unquoted_underscores keeps labels identical to the FASTA IDs (no added quotes)
    dated_tree.write(path=output_file, schema=schema, unquoted_underscores=True)

def simulate_sequences_on_tree(tree_file: str,
                               tree_id: int,
                               root_sequence: str,
                               sub_model: str="HKY",
                               clock_rate: float=1.0,
                               frequencies: list[float]=[0.25, 0.25, 0.25, 0.25],
                               kappa: Optional[float]=2.0,
                               mutation_rates: Optional[dict]=None,
                               gamma_dist: dict={"alpha": 0.5, "num_categories": 4},
                               output_dir: str="./",
                               root_date: float=2002.0,
                               output_tree_file: Optional[str]=None,
                               tree_schema: str="newick"
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
        output_tree_file (str, optional): Path where the tree with dated tip names is saved.
            If None, it is saved as "tree_with_dates.nwk" (or ".nex" for nexus) inside the
            tree_<id> output directory.
        tree_schema (str): Format of the saved tree, "newick" or "nexus". Defaults to "newick".
        
    Returns:
        None. Generates sequence, info, and rate files, plus the dated tree, in the output directory.
    """
    n_sites = len(root_sequence)
    
    if sub_model == "HKY":
        model = pyvolve.Model("nucleotide",
                              {"kappa": kappa, "state_freqs": frequencies},
                              alpha=gamma_dist["alpha"],
                              num_categories=gamma_dist["num_categories"])
    elif sub_model == "GTR":
        model = pyvolve.Model("nucleotide", 
                            {"mu": mutation_rates, "state_freqs": frequencies},
                            alpha=gamma_dist["alpha"],
                            num_categories=gamma_dist["num_categories"])
    else:
        raise ValueError("The sub_model parameter has to be either HKY or GTR")
    
    rate_factors = model.rate_factors
    rate_probs = model.rate_probs
    
    rng = np.random.default_rng(seed=42)
    site_categories = rng.choice(gamma_dist["num_categories"], size=n_sites, p=rate_probs)
    
    partitions, site_order = [], []
    for k in range(gamma_dist["num_categories"]):
        idx = np.nonzero(site_categories == k)[0]
        
        if len(idx) == 0:
            continue
        
        sub_seq = "".join(root_sequence[i] for i in idx)
        sub_model = pyvolve.Model("nucleotide", rate_factors=[rate_factors[k]])
        
        partitions.append(pyvolve.Partition(models=sub_model, root_sequence=sub_seq))
        site_order.extend(idx)   # remember where these columns belong in the original sequence
        
    # partition = pyvolve.Partition(models=model, root_sequence=root_sequence)
        
    dp_trees = dendropy.TreeList.get(path=tree_file, schema="nexus")
    
    tree = pyvolve.read_tree(tree=dp_trees[tree_id].as_string(schema="newick"), scale_tree=clock_rate)
    evolver = pyvolve.Evolver(partitions=partitions, tree=tree)
    
    # Output files for each tree
    tree_dir = Path(output_dir) / f"tree_{tree_id}"
    tree_dir.mkdir(parents=True, exist_ok=True)
    
    fasta_file = tree_dir / "sequences.fasta"
    info_file = tree_dir / "sequences_info.txt"
    rates_file = tree_dir / "sequences_rates.txt"

    evolver(seqfile=None, infofile=None, ratefile=None)
    raw_seqs = evolver.get_sequences()
    
    final_seqs = {}
    for taxon, seq in raw_seqs.items():
        ordered = [None] * n_sites
        for col_pos, orig_pos in enumerate(site_order):
            ordered[orig_pos] = seq[col_pos]
        final_seqs[taxon] = "".join(ordered)
        
    records = [SeqRecord(Seq(seq), id=taxon, description="") 
            for taxon, seq in final_seqs.items()]

    SeqIO.write(records, fasta_file, "fasta")
    add_dates(dp_trees[tree_id], fasta_file, root_date)

    if output_tree_file is None:
        ext = "nex" if tree_schema == "nexus" else "nwk"
        output_tree_file = tree_dir / f"tree_with_dates.{ext}"
    save_tree_with_dates(dp_trees[tree_id], output_tree_file, root_date, tree_schema)

    with open(rates_file, "w") as f:
        f.write("Site_Index\tPartition_Index\tRate_Category\n")
        for site_idx in range(n_sites):
            category = int(site_categories[site_idx]) + 1  # pyvolve uses 1-based indexing
            f.write(f"{site_idx + 1}\t1\t{category}\n")

    # site_rates_info.txt  — one row per rate category
    with open(info_file, "w") as f:
        f.write("Partition_Index\tModel_Name\tRate_Category\tRate_Probability\tRate_Factor\n")
        for k in range(gamma_dist["num_categories"]):
            f.write(f"1\tnucleotide\t{k + 1}\t{rate_probs[k]}\t{rate_factors[k]}\n")