#!/bin/bash
#SBATCH --job-name=scn2d_rep4
#SBATCH --output=scn2d_rep4.out
#SBATCH --error=scn2d_rep4.err
#SBATCH --time=1500:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

module load beast/1.8
cd /project/cobey/mvieira/evolution_of_mutability/results/BEAST/simulated_alignments/scenario2d_rep4/
beast /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments/scenario2d_rep4/scenario2d_rep4.xml
