from Bio import AlignIO, Phylo
from Bio.Phylo.Consensus import bootstrap_trees, get_support
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt

# 1. Load the alignment (format can be "fasta", "phylip", "clustal", "nexus")
alignment = AlignIO.read("data_cluster/first_multiple_chain_run/data/simulated_sequences/tree_3/sequences.fasta", "fasta")
print(f"Loaded {len(alignment)} sequences, {alignment.get_alignment_length()} columns")

# 2. Set up distance calculation and NJ tree construction
calculator = DistanceCalculator("identity")
constructor = DistanceTreeConstructor(calculator, "nj")

# 3. Build the reference NJ tree
nj_tree = constructor.build_tree(alignment)

# 4. Bootstrap (100-1000 replicates is typical)
boot_trees = list(bootstrap_trees(alignment, 100, constructor))

# 5. Add bootstrap support (in %) to the reference tree
support_tree = get_support(nj_tree, boot_trees)
support_tree.ladderize()

# 6. Save the tree
Phylo.write(support_tree, "testing/nj_tree_bootstrap.nwk", "newick")

# 7. Plot, showing only terminal names and bootstrap values on internal nodes
def tip_label(clade):
    return clade.name if clade.is_terminal() else None

def support_label(clade):
    if clade.is_terminal() or clade.confidence is None:
        return None
    return f"{clade.confidence:.0f}"

fig, ax = plt.subplots(figsize=(10, 8))
Phylo.draw(support_tree, axes=ax, do_show=False,
           label_func=tip_label, branch_labels=support_label)
plt.savefig("testing/nj_tree_bootstrap.png", dpi=300, bbox_inches="tight")
plt.show()