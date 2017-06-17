# Evolution of mutability
This is the repository for *Mutability loss is intrinsic to the evolution of long-lived B cell lineages* (Vieira et al. *in prep*). This document outlines the overall pipeline of the project. Detailed instructions for replicating each analysis can be found in the README files in the corresponding directories in the ```Analyses``` folder.

The project's analyses are divided into two main stages. The first stage is to align B cell receptor sequences and fit phylogenetic models to the resulting alignments using [BEAST (v.1.8.2)](http://http://beast.bio.ed.ac.uk/). Before they can be aligned, heavy chain sequences from the VRC01 dataset must be partitioned and assigned to different clones using [Partis](http://https://github.com/psathyrella/partis).

The second stage is to do a series of analyses on the phylogenetic trees inferred by BEAST:

**Relative mutability**: comparison between the observed mutability of B cell receptor sequences and their expected mutability under codon randomization.

**Mutability vs. time**: changes over time in the average mutability of nodes in the B cell trees.

**Contrasts**: quantification of the magnitude, frequency and distribution of mutability changes across branches of the B cell trees.

**Selection**: testing for positive and purifying selection on B cell receptor sequences using BASELINe.

**Divergence from ancestors**: measuring amino acid divergence between B cell receptor sequences and their unmutated ancestors.

**Synonymous / Non-Synonymous mutability changes**: quantifying mutability changes caused by synonymous and non-synonymous substitutions separately; comparing observed synonymous changes to simulations under different models.

In addition to the code for conducting the analyses, this repository contains a mostly empty results directory with the hierarchy required for analysis scripts to access results from other analyses and output their own. Sequence alignment results (and fasta files with the partitioning of VRC01 sequences) are included in the repository.
