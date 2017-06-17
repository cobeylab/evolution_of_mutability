## Clonal assignment ##
We used [Partis](http://github.com/psathyrella/partis/) to partition VRC01 heavy sequences into clones and to identify their likely germline V and J genes.

Instructions for installing Partis can be found [here](http://github.com/psathyrella/partis/blob/master/manual.md). 

To partition VRC01 heavy chain sequences into clones, go to the Partis directory and run (as a single command):

```./bin/partis partition --n-procs <number of processors> --infname /path/to/evolution_of_mutability/Data/Wu_2015_VRC01/Wu_2015_VRC01_heavy_curated.fasta --outfname /path/to/evolution_of_mutability/results/clonal_assignment/VRC01/VRC01_partition.csv```

To run Partis to identify germline genes, run (as a single command):

```./bin/partis run-viterbi --n-procs <number of processors> --infname path/to/evolution_of_mutability/Data/Wu_2015_VRC01/Wu_2015_VRC01_heavy_curated.fasta --outfname path/to/evolution_of_mutability/results/clonal_assignment/VRC01/VRC01_annotation.csv```

To process Partis output and get sequences for each clone into separate fasta files (with germline gene information in corresponding csv files), go to the ```analyses/clonal_assignment``` directory and run ```partition_into_clones.py```.

The resulting fasta and csv files for each clone are provided in the ```results/clonal_assignment``` directory. Because of the stochasticity in the estimation procedure used by Partis, different runs may provide somewhat different results.



