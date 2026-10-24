"""
Exercise 1 — Building a Neighbor-Joining tree

Instructions (to follow in the lab):
1. Reuse sequences from previous labs (FASTA from Lab 2 or FASTQ→FASTA from Lab 3).
2. If you only have single-sequence FASTA files, combine at least 3 into a multi-FASTA file.
3. Save the multi-FASTA file in: data/work/<handle>/lab04/your_sequences.fasta
4. Complete the steps below:
   - load the multi-FASTA file,
   - compute the distance matrix,
   - build the NJ tree,
   - save the result in Newick format (.nwk).
"""

from pathlib import Path
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

if __name__ == "__main__":
    # TODO 1: Load your multi-FASTA file
    handle = "AlexZvinca"
   
    fasta = Path(f"data/work/{handle}/lab04/aligned_sequences.fasta")
    output = Path(f"/workspaces/BDHB-lab-Zvinca/labs/04_phylogenetics/submissions/{handle}/tree_AlexZvinca.nwk")

    # Load aligned sequences
    alignment = AlignIO.read(fasta, "fasta")

    # TODO 2: Compute the distance matrix
    calculator = DistanceCalculator("identity")
    distance_matrix = calculator.get_distance(alignment)
        
    print("Distance matrix:")
    print(distance_matrix)

    # TODO 3: Build the NJ tree
    constructor = DistanceTreeConstructor()
    tree = constructor.nj(distance_matrix)
    
    # TODO 4: Save the tree in Newick format
    Phylo.write(tree, output, "newick")

    # TODO 5: Visualize the tree
    Phylo.draw_ascii(tree)
    print(f"Tree saved to: {output}")
