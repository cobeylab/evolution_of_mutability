## BEAST analysis ##

We used [BEAST (v.1.8.2)](http://beast.bio.ed.ac.uk/) to fit Bayesian phylogenetic models to the sequence alignments for each of the seven heavy and light chain datasets (which we refer to as "clones" in the description of file nomenclature below): CH103, CH103L, VRC26, VRC26L, VRC01_13, VRC01_01, VRC01_19.

The analysis of each clone involves producing an XML file with data and model specification and passing it to BEAST. We used three different population dynamic models to compute prior tree probabilities: constant population size, logistic growth, and exponential growth. Those priors are identified in the XML files by their names ("constant", "exponential", "logistic") or abbreviations ("con", "exp", "log"). We ran four replicate MCMC chains under a constant population size prior, and two chains each under the exponential and logistic priors. Each replicate chain is identified by a letter (a-d) (note that a single XML file is used for all replicate chains).

To run a single MCMC chain for the CH103 dataset using a constant population size prior, for example, run:

```beast CH103_con_run1.xml```

Subsequent analysis assume BEAST results will be in the ```results``` under the following hierarchy:

```evolution_of_mutability/results/BEAST/observed_lineages/<CLONE>_<prior>/<CLONE>_<prior abbreviation>_run1<replicate chain letter>```. 

For example, to have results from chain ```CH103_con_run1a``` in the correct location, go to ```results/BEAST/observed_lineages/CH103_constant/CH103_con_run1a``` and run:

```../../../../../analyses/BEAST/observed_lineages/CH103_constant/CH103_con_run1.xml```

We describe the steps for setting up the XML files below, but also provide the files in the corresponding clone + prior folders in the ```analyses/BEAST``` directory.

To process the output from BEAST, we implemented a shell script (```process_output.sh```) that downsamples the ```.trees``` file to retain a sample of 1,000 trees from the posterior distribution, and uses the "loganalyzer" and "treeannotator" programs provided with BEAST to summarize beast ```.log``` files and reconstruct a maximum-clade credibility tree, respectively. The shell script requires a file called ```burnins_observed.csv``` to be located in the ```results/BEAST/observed_lineages``` folder. The ```burnins_observed.csv``` file must contain four columns, specifying the chain id (in the ```<CLONE>_<prior abbreviation>_run1<replicate chain letter>``` format), the number of steps in the MCMC chain, the number of steps in the burn-in, and the number of lines before the first line containing a tree in the BEAST ```.trees``` files. To execute the bash script, run:

```./process_output.sh <CLONE> <prior>```.

Note that the BEAST commands in ```process_output.sh``` must be modified to specify the path to the BEAST installation in the user's computer.

## Setting up XML files  ##

1 – Import ```<CLONE>_final_aligment.nex``` into BEAUti v.1.8.2.

2 – On the *Taxa* tab, set up a taxon set containing all sequences except the GERMLINE, name it *ingroup* and set it to be monophyletic.

3 – On the *Tips* tab, use "Guess Dates" to extract sampling time from sequence IDs (defined by the underscore).

4 – On the *Sites* tab, set the substitution model to GTR, with *empirical* base frequencies (to reduce the number of estimated parameters), turn-on partition into 3 codon positions and uncheck the first box (again, to reduce the number of parameters).

5 – On the *Clocks* tab, set the model to *Random local clock*.

6 – On the *Trees* tab, set the tree prior to coalescent, logistic or exponential. Leave the parameterization at the default option.

7 – On the *States* tab, check *Reconstruct states at all ancestors* and *Reconstruct synonymous/non-synonymous change counts*.

8 – On the *Priors* tab, leave default priors, except where specified below:

| Parameter name | Description                        | Prior distribution | Initial value |
|----------------|------------------------------------|--------------------|---------------|
| CP1.mu         | Relative rate for codon position 1 | Gamma [0.05, 10]   | 1             |
| CP2.mu         | Relative rate for codon position 2 | Gamma [0.05, 10]   | 1             |
| CP3.mu         | Relative rate for codon position 3 | Gamma [0.05, 10]   | 1             |
| clock. rate    | Substitution rate                  | Uniform [0,10]      | 0.0001
| logistic.popSize    | Population size under logistic growth  | Default| 100        
| logistic.growthRate    | Logistic growth rate  | Default| 0.5


9 – On the *MCMC* tab, set chain length to 500,000,000 and leave sampling frequency at the default 1 every 1,000 value. Set the file name stem to <CLONE>\_<con/log/exp>\_run1

10 – Save XML file to the ```...analyses/BEAST/observed_alignments/<CLONE>_<constant/logistic/exponential>``` folder.

11 – Edit XML file to replace ```<logTree id="treeFileLog" logEvery="1000"``` by ```<logTree id="treeFileLog" logEvery="10000"```. This will make trees be recorded every 10,000 steps even though the parameters are recorded every 1,000 steps.

