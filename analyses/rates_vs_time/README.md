## Rates vs. time ##

```rates_vs_time.py``` takes each tree in the subsample of trees from the posterior distribution obtained by BEAST and records substitution rates inferred for each branch and the time (since the tree's root) at each branch's parent node.

As an example, to perform the analysis for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python rates_vs_time.py CH103_con_run1a```

The script will analyze the trees in the file ```CH103_con_run1a_tree_sample.trees``` (in the corresponding BEAST results folder) and output two csv files to the ```results/rates_vs_time/observed_lineages/CH103_constant/``` folder: ```CH103_con_run1a_rates_vs_time_points.csv``` lists all branches for all trees with their corresponding values for substitution rates and time from the root, whereas ```CH103_con_run1a_rates_vs_time_correlations.csv``` lists, for each tree, results from a linear regression analysis between substitution rates on a branch (response variables) and time (predictor variable) across branches in that tree.
