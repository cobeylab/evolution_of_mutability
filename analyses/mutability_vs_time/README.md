## Mutability vs time ##

```mutability_vs_time.py``` takes each tree in the subsample of trees from the posterior distribution obtained by BEAST and records the motif-based sequence mutability and the time (since the tree's root) for each node.

As an example, to perform the analysis for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python mutability_vs_time.py CH103_con_run1a```

The script will analyze the trees in the file ```CH103_con_run1a_tree_sample.trees``` (in the corresponding BEAST results folder) and output two csv files to the ```results/mutability_vs_time/observed_lineages/CH103_constant/``` folder: ```CH103_con_run1a_mutability_vs_time_points.csv``` lists all nodes for all trees with their corresponding values for sequence mutability and time from the root, whereas ```CH103_con_run1a_mutability_vs_time_correlations.csv``` lists, for each tree, results from a linear regression analysis between sequence mutability (response variable) and time (predictor variable) across nodes in that tree.
