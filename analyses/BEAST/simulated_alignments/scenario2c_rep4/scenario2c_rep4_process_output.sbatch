#!/bin/bash
#SBATCH --job-name=scenario2c_rep4_process_output
#SBATCH --output=scenario2c_rep4_process_output.out
#SBATCH --error=scenario2c_rep4_process_output.err
#SBATCH --time=3:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=10000

module load java
cd /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/simulated_alignments
./process_output.sh 2c 4
