## Selection ##

We analyze B cell receptor sequence alignments with [BASELINe](http://selection.med.yale.edu/baseline/). BASELINe compares the observed ratio of non-synonymous to synonymous substitutions between each sequence and its unmutated common ancestor model to the expected ratio under neutral mutations following the S5F mutability model.

To prepare fasta files for each dataset by including the unmutated common ancestors inferred from BEAST to the dataset's sequence alignment, run:

```python prepare_fasta_for_baseline.py```

Resulting fasta files are output into the ```baseline_fasta_files``` folder. 

To analyze an individual alignment using baseline, for instance the CH103 heavy chain dataset, run:

```RScript Baseline_Main_Version1.3.R 1 1 5 5 0 0 1:11:19:36:43:81:96 ./baseline_fasta_files/CH103_con_run1a_baseline.fasta ../../results/selection/CH103_con_run1a _baseline```

For a description of the first six arguments after the ```RScript``` command, see the comments in ```Baseline_Main_Version1.3.R```. The argument ```1:11:19:36:43:81:96``` indicates partition points delimiting FRs and CDRs (for heavy chain CH103 specifically, in this example). The first number indicates the first amino acid position of FR1, and the subsequent numbers indicate the last amino acid positions of FRs and CDRs 1-3. A Python script with partition points for each dataset (in nucleotide positions), ```partition_points.py``` can be found in the ```evolution_of_mutability/mutability``` directory. 

A ```.txt`` file with the results is output to the ```results/selection/``` folder, but subsequent analyses require the output files to have the ```.tsv``` extension:

```mv ../../results/selection/CH103_con_run1a_baseline.txt ../../results/selection/CH103_con_run1a_baseline.tsv```

The shell script ```run_baseline.sh``` executes the two commands above for all datasets.
