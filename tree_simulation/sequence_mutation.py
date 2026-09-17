import re
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

    output_fasta = f"{input_fasta.split('.')[0]}_with_dates.fasta"
    SeqIO.write(records, output_fasta, "fasta")

def simulate_sequences_on_tree(tree_file: str,
                               root_sequence: str,
                               clock_rate: float=1.,
                               frequencies: list[float]=[0.25, 0.25, 0.25, 0.25],
                               mutation_rates: dict={"AC":1.0, "AG":1.0, "AT":1.0, "CG":1.0, "CT":1.0, "GT":1.0},
                               gamma_dist: dict={"alpha": 0.5, "num_categories": 4},
                               output_dir: str="./",
                               root_date: float=2002.0
                               ) -> None:
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
        fasta_file = f"{output_dir}/simulated_sequences_{i}.fasta"
        info_file = f"{output_dir}/simulated_sequences_info_{i}.txt"
        rates_file = f"{output_dir}/simulated_sequences_rates_{i}.txt"

        evolver(seqfile=fasta_file, seqfmt="fasta", 
                infofile=info_file, 
                ratefile=rates_file)
        
        add_dates(dp_tree, fasta_file, root_date)
    

def main():
    # remove_coma_in_nexus_file("./remaster_sims/sim_logs/sim.full.trees", "./remaster_sims/sim_logs/sim.full_fixed.trees")
    
    custom_mu = {"AC":1, "AG":2, "AT":1, "CG":2, "CT":3, "GT":1}
    freqs = [0.1, 0.45, 0.3, 0.15] # A, C, G, T
    gamma_dist = {"alpha": 0.5, "num_categories": 4}
    clock_rate = 1e-2 # Substitutions per site per unit time
    root_sequence = "A"*int(1e4)  # Root sequence of the tree
    root_date = 2002.0
    
    simulate_sequences_on_tree(tree_file="./remaster_sims/sim_logs/sim.sampled.trees",
                               root_sequence=root_sequence,
                               clock_rate=clock_rate,
                               frequencies=freqs,
                               mutation_rates=custom_mu,
                               gamma_dist=gamma_dist,
                               output_dir="test_data",
                               root_date=root_date)
        
    

if __name__ == "__main__":
    main()
