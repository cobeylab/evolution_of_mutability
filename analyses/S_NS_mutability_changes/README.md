## Synonymous / Non-synonymous mutability changes ##

We implemented ```S_NS_mutability_changes.py``` to partition changes in S5F mutability on all branches of maximum-clade-credibility trees into a component attributable to synonymous subsitutions and a component attributable to non-synonymous substitutions. In addition, the script simulates sequences descending from each node in the tree under motif-based on motif-independent models of mutation rate variation across sites in B cell receptor sequences, and partitions simulated mutability changes into synonymous and non-synonymous components.

As an example, to perform the analysis for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python S_NS_mutability_changes.py CH103_con_run1a```

The script will analyze the MCC tree in the file ```CH103_con_run1a_MCC_tree.tree``` (in the corresponding BEAST results folder) and output two csv files to the ```results/S_NS_mutability_changes/observed_lineages/CH103_constant/``` folder:  ```CH103_con_run1a_observed_MCC.csv``` contains results for the original observed and inferred sequences from the MCC tree, whereas ```CH103_con_run1a_simulated_MCC.csv``` contains results from the simulated descendant sequences.
