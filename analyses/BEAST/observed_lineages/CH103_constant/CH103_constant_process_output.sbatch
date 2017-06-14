#!/bin/bash
#SBATCH --job-name=CH103_constant_process_output
#SBATCH --output=CH103_constant_process_output.out
#SBATCH --error=CH103_constant_process_output.err
#SBATCH --time=2:00:00
#SBATCH --partition=cobey
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=16000

module load java
cd /project/cobey/mvieira/evolution_of_mutability/analyses/BEAST/observed_lineages
./process_output.sh CH103 constant
