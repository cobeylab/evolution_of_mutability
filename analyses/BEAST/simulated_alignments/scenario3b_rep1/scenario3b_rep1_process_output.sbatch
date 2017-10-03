#!/bin/bash
#SBATCH --job-name=scenario3b_rep1_process_output
#SBATCH --output=scenario3b_rep1_process_output.out
#SBATCH --error=scenario3b_rep1_process_output.err
#SBATCH --time=3:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load java
cd /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments
./process_output.sh 3b 1
