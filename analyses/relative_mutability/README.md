## Relative mutability ##

```relative_mutability.py``` takes each node in the maximum-clade-credibility of the B cell trees, records its observed mutability and produces a distribution of mutability values obtained by randomizing the original codon sequence while keeping the amino acid sequence constant. Codons are sampled during the randomization according to their usage in human genes.

As an example, to perform the analysis for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

``python relative_mutability.py CH103_con_run1a```

The script will analyze the MCC tree in the file ```CH103_con_run1a_MCC_tree.tree``` (in the corresponding BEAST results folder) and output two csv files to the ```results/relative_mutability/observed_lineages/CH103_constant/``` folder:  ```CH103_con_run1a_observed_mutability.csv``` contains results for the original observed and inferred sequences from the MCC tree, whereas ```CH103_con_run1a_randomized_mutability.csv``` contains results from the randomized sequences for each node.
