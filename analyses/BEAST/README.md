#BEAST analyses#

##Overview##

We performed the same analysis on each observed B cell clone, under three alternative choices for the coalescent prior (constant population size, logistic growth, exponential growth).

The analysis of each clone (observed or simulated) involves producing an XML file and passing it to BEAST v.1.8.


##Setting up XML files for observed alignments##
Because the germline-as-outgroup constraint is not obviously easy to enforce for observed alignments using BEASTGen, XML files for observed alignments were set up manually in BEAUti v.1.8 following the steps below (which produce the same configuration used for XML files of simulated alignments, except for a longer chain length of 500 million).

1 – Import ```*_final_aligment.nex``` into BEAUti v.1.8.2.

2 – On the *Taxa* tab, set up a taxon set containing all sequences except the GERMLINE, name it *ingroup* and set it to be monophyletic.

3 – On the *Tips* tab, use "Guess Dates" to extract sampling time from sequence IDs (defined by the underscore).

4 – On the *Sites* tab, set the substitution model to GTR, with *empirical* base frequencies (to reduce the number of estimated parameters), turn-on partition into 3 codon positions and uncheck the first box (again, to reduce the number of parameters).

5 – On the *Clocks* tab, set the model to *Random local clock*.

6 – On the *Trees* tab, set the tree prior to coalescent, logistic or exponential. Leave the parameterization at the default option.

7 – On the *States* tab, check *Reconstruct states at all ancestors* and *Reconstruct synonymous/non-synonymous change counts*.

8 – On the *Priors* tab, leave default priors, except where specified below:

| Parameter name | Description                        | Prior distribution | Initial value |   |   |
|----------------|------------------------------------|--------------------|---------------|---|---|
| CP1.mu         | Relative rate for codon position 1 | Gamma [0.05, 10]   | 1             |   |   |
| CP2.mu         | Relative rate for codon position 2 | Gamma [0.05, 10]   | 1             |   |   |
| CP3.mu         | Relative rate for codon position 3 | Gamma [0.05, 10]   | 1             |   |   |
| clock. rate    | Substitution rate                  | Uniform[0,10]      | 0.0001
| logistic.popSize    | Population size under logistic growth  | Default| 100        
| logistic.growthRate    | Logistic growth rate  | Default| 0.5
|   |   |

9 – On the *MCMC* tab, set chain length to 500,000,000 and leave sampling frequency at the default 1 every 1,000 value. Set the file name stem to <CLONE>\_<con/log/exp>\_run1

10 – Save XML file to the ```...analyses/BEAST/observed_alignments/<CLONE>_<constant/logistic/exponential>``` folder.

11 – Edit XML file to replace ```<logTree id="treeFileLog" logEvery="1000"``` by ```<logTree id="treeFileLog" logEvery="10000"```. This will make trees be recorded every 10,000 steps even though the parameters are recorded every 1,000 steps.

##Running the chains##
For each clone (treating heavy and light chains as independent clones) and coalescent prior, we ran replicated MCMC chains using BEAST (4 chains per clone under a constant population size coalescent prior, 2 chains each for logistic and exponential priors).

Replicate chains are identified by letters (*a*, *b*, *c*, *d*).

To start a chain, run ```<CLONE>_<con/log/exp>_run1<a/b/c/d>.sbatch```.


##Processing BEAST output##
We processed the output from the BEAST chains in order to summarize the posterior distributions for the parameters (using the ```logananyzer``` program associated with BEAST) and to extract a sub-sample of 1000 trees from the complete posterior sample of trees. For the sub-sample, we took trees at regular intervals between the first state after the burn-in and the last state of the chain. For each observed-alignment chain, the burn-in was determined visually and by looking at the effective sample sizes of parameters in Tracer v.1.6.0 and recorded in ```results/BEAST/observed_lineages/burnins_observed.csv```. Chain lengths and burn-ins for chains corresponding to simulated alignments are listed in ```results/BEAST/simulated_alignments/burnins_simulated.csv```

The sbatch files for processing beast output are ```<CLONE>_<con/log/exp>_run<1/2/3...><a/b/c/d>_proc_output.sbatch```.
Those invoke the ```process_output.sh``` file, located in the ```analyses/BEAST/observed_lineages``` folder.

