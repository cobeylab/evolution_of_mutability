## Contrasts ##

```contrasts.py``` quantifies the magnitude, frequency and distribution of mutability changes across branches of the B cell trees in the downsampled set of trees from an MCMC chain. The script can be executed by running:

```python contrasts.py <chain_id>```

where ```<chain_id``` is an MCMC chain identifier. For example, to analyze contrasts for replicate chain "a" of CH103 heavy-chain under a constant growth coalescent prior, run:

```python contrasts.py CH103_con_run1a```

The script will analyze trees from ```CH103_con_run1a_tree_sample.trees``` (in the corresponding BEAST results folder) and output a csv file with the results called ```CH103_con_run1a_mutability_contrasts.csv``` to the ```results/contrasts/observed_lineages/CH103_constant``` folder.
