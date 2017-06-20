## Divergence from ancestors ##

```divergence_from_ancestors.py``` takes a maximum-clade-credibility tree of a B cell lineage and computes the amino acid and nucleotide percent divergence between each node in the tree and the root node (representing the unmutated common ancestor of the B cell lineage).

For example, to perform the analysis for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python divergence_from_ancestors.py CH103_con_run1a```

The script will analyze the MCC tree in the file ```CH103_con_run1a_MCC_tree.tree``` (in the corresponding BEAST results folder) and output results to ```results/divergence_from_ancestors/CH103_con_run1a_divergence_from_ancestors.csv```.
