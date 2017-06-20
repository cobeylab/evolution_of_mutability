## Tree annotation ##

```annotate_MCC.py``` takes the NEXUS file containing the maximum-clade-credibility (MCC) tree from an MCMC chain and outputs a NEXUS file with nodes annotated with mutability metrics based on their sequences.

As an example, to annotate the MCC tree of replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python annotate_MCC.py CH103_con_run1a```

The script will analyze the MCC tree in the file ```CH103_con_run1a_MCC_tree.tree``` (in the corresponding BEAST results folder) and output the annotated NEXUS file, ```CH103_con_run1a_annotated_MCC.tree```, to the ```results/tree_annotation``` directory.

