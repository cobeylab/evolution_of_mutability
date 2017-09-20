#!/bin/bash
#SBATCH --job-name=scn4a_rep2
#SBATCH --output=scn4a_rep2.out
#SBATCH --error=scn4a_rep2.err
#SBATCH --time=1500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8
cd /project/cobey/mvieira/evolution_of_mutability/results/BEAST/simulated_alignments/scenario4a_rep2/
beast /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments/scenario4a_rep2/scenario4a_rep2.xml
